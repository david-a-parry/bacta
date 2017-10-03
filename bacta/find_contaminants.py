import sys
import os
from subprocess import PIPE, Popen, call
import re
import pysam
import tempfile
import gzip
import shutil
import logging
import multiprocessing as mp
from itertools import repeat
from collections import defaultdict
from .cigar_scorer import CigarScorer
from .reverse_complement import reverse_complement
from .alignment_file_utils import get_align_output,get_bamfile
from .alignment_file_utils import seek_back_til_reads, get_tmp_bam

#CIGAR STRING OPERATORS and BAM CODES
#    M   BAM_CMATCH  0
#    I   BAM_CINS    1
#    D   BAM_CDEL    2
#    N   BAM_CREF_SKIP   3
#    S   BAM_CSOFT_CLIP  4
#    H   BAM_CHARD_CLIP  5
#    P   BAM_CPAD    6
#    =   BAM_CEQUAL  7
#    X   BAM_CDIFF   8
#    B   BAM_CBACK   9

_qname_re = re.compile(r"""(\S+)(#\d+)?(/\[1-2])?""")  
__version__ = '0.0.1'

def _get_read_coord(read, no_commas=False):
    com = ','
    if no_commas:
        com = ''
    if read.reference_id > -1:
        coord = "{}:{:{n_format}}".format(read.reference_name, 
                                          read.reference_start + 1, 
                                          n_format=com)
    else:
        coord = "*/*"
    return coord

def _overlaps_loci(read, target_loci):
    start,end = _infer_bounds(read)
    for locus in target_loci:
        #contig is guaranteed to be same for read and target_loci
        if start < locus and end > locus:
            return True
        if read.is_paired:
            #want to avoid fetching mate (slow)
            if read.has_tag('MC'):
                # this is a standard tag, should not need to check type
                mctups = cigar_scorer.cigarstring_to_tuples(read.get_tag('MC'))
                s_off, e_off = cigar_scorer.get_clipped_offset(mctups)
                start = read.next_reference_start - s_off
                end = (cigar_scorer.get_aligned_length(mctups) + start + e_off)
                if start <= locus and end > locus:
                    return True
            else:
                # in absence of MC tag we'll have to make a guess and fix 
                # any unpaired overlapping reads using _mop_up_region_pairs
                start = read.next_reference_start
                end = start + read.infer_read_length()
                if start < locus and end > locus:
                    return True
    return False

def _infer_bounds(read):
    if read.cigartuples is None: #unmapped
        return (0, 0)
    s_off, e_off = cigar_scorer.get_clipped_offset(read.cigartuples)
    start = read.reference_start - s_off
    end = read.reference_end + e_off
    return (start, end)

def check_read_clipping(read, cigar_scorer, min_frac, min_clip=None, 
                        mate_cigar=None):
    ''' 
        Returns True if read is clipped greater than 
        min_fraction_clipped or min_bases_clipped. If cigarstring 
        for mate is provided, returns True if mate is clipped 
        greater than min_fraction_clipped or min_bases_clipped.

        Args:
            read:   A pysam.AlignedSegment object representing a 
                    single read. The cigar string will be read if 
                    available to determine the extent of soft/hard
                    clipping.

            cigar_scorer:
                    CigarScore object.

            min_frac:
                    Minimum proportion of the read that must be clipped 
                    or mismatch the reference. Function will return True
                    if read clipping/mismatch reaches this threshold.

            min_clip:
                    Minimum actual number of bases clipped. Method will
                    return True if the number of clipped/mismatching 
                    bases reaches this threshold. By default this value
                    is None and is not used. If either min_clip or 
                    min_frac conditions are met, this method will return
                    True.

            mate_cigar:
                    Optional cigar string for mate (e.g. as provided
                    by the standard 'MC' tag. If provided, this 
                    function will return True if either the read or
                    mate cigar indicate clipping above the threhold.
            
    '''
    if read.cigartuples is not None: #check if mapped instead?
        md = None
        if read.has_tag('MD'):
            md = read.get_tag('MD')
        if check_cigar_tuple_clipping(read.cigartuples, cigar_scorer, min_frac,
                                      min_bases_clipped=min_clip, md_tag=md, ):
            return True
    if mate_cigar is not None:
        return check_cigarstring_clipping(mate_cigar, cigar_scorer, min_frac,
                                          min_clip=min_clip)
    return False

def check_cigarstring_clipping(cigar, cigar_scorer, min_frac, min_clip=None,
                               md_tag=None):
    ''' 
        Performs the same function as check_read_clipping but 
        operates directly on the cigarstring.
    '''
    cts = cigar_scorer.cigarstring_to_tuples(cigar)
    return check_cigar_tuple_clipping(cts, cigar_scorer, min_frac, 
                                      min_bases_clipped=min_clip, 
                                      md_tag=md_tag)
    
def check_cigar_tuple_clipping(cts, cigar_scorer, min_fraction_clipped,
                                min_bases_clipped=None, md_tag=None):
    ''' 
        Returns True if read is clipped greater than 
        min_fraction_clipped or min_bases_clipped. If md_tag is 
        provided, the number of single nucleotide mismatches are 
        also added to the number of clipped bases.
    '''
    clipping = 0
    length = 0
    for c in cts:
        if 4 <= c[0] <= 5:
            #SOFT or HARD clip
            clipping += c[1]
            length += c[1]
        elif c[0] < 2 or 7 <= c[0] <= 8:
            #MATCH, INS, EQUAL or DIFF
            length += c[1]
    if md_tag is not None:
        clipping += cigar_scorer.md_mismatches(md_tag)
    if clipping:
        if (min_bases_clipped is not None and 
            clipping >= min_bases_clipped):
            return True
        else:
            if float(clipping)/length >= min_fraction_clipped:
                return True
    return False


def should_store_is_clipped(read, cigar_scorer, min_frac, min_clip=None,
                            decoys=None, get_unaligned=False):
    ''' 
        See if we have enough information to determine whether read
        should be stored or not. Returns two booleans - the first 
        indicates whether the read should be stored prior to 
        encountering its mate, the second indicates whether the read
        is clipped beyond the clipping threshold.

        If the MC (mate cigar) tag is present, we can determine 
        whether the mate is not clipped and thus can be ignored if
        this read is also not clipped.
    '''
    store = True
    is_clipped = False
    if read.reference_id == -1: #unmapped
        if get_unaligned:
            is_clipped = True
    elif decoys and read.reference_name in decoys:
        is_clipped = True
    elif read.mate_is_unmapped:
        if check_read_clipping(read, cigar_scorer, min_frac, 
                               min_clip=min_clip):
            is_clipped = True
        else:
            store = False
    elif read.has_tag('MC'):
        if check_read_clipping(read, cigar_scorer, min_frac, min_clip=min_clip,
                               mate_cigar=read.get_tag('MC')):
            is_clipped = True
        else:
            store = False
    elif check_read_clipping(read, cigar_scorer, min_frac, min_clip=min_clip):
        is_clipped = True
    return (store, is_clipped)

def _merge_intervals(intervals, flanks=0):
    if not intervals:
        return [] 
    merged = []
    intervals = sorted(intervals)
    prev_interval = intervals[0]
    for i in range(1, len(intervals)):
        if prev_interval.contig != intervals[i].contig:
            merged.append(prev_interval)
            prev_interval = intervals[i]
        elif intervals[i].start - flanks <= prev_interval.end + flanks:
            if intervals[i].targets:
                prev_interval.targets.update(intervals[i].targets)
            if prev_interval.end < intervals[i].end:
                prev_interval.end = intervals[i].end
        else:
            merged.append(prev_interval)
            prev_interval = intervals[i]
    merged.append(prev_interval)
    return merged
        
def _mop_up_region_pairs(bamfile, pair_dict, clipped_names, fq_writer, logger):
    ''' 
        Fetch unencountered mates in an efficient manner and score/
        output.

        Args:
            pair_dict:  
                dict of read ID to the previously encountered read.

            clipped_names:
                read IDs for any reads that have already exceeded 
                clipping threshold and should be written to output.
     '''
    to_fetch = []
    logger.info("{:,} unpaired reads after parsing region"
                     .format((len(pair_dict[p]) + len(pair_dict[p]))))
    for name,read in (list(pair_dict[1].items()) + 
                      list(pair_dict[2].items())):
        if name in clipped_names:
            to_fetch.append(GenomicInterval(read.reference_name, 
                                            read.reference_start,
                                            read.reference_end or 
                                                read.reference_start, 
                                            name))
        elif read.has_tag('MC'):
            # already got mate cigar string - check before fetching
            # this is a standard tag, should not need to check type
            mcigar = read.get_tag('MC')
            if check_cigarstring_clipping(mcigar):
                clipped_names.add(name)
                to_fetch.append(GenomicInterval(read.reference_name, 
                                                read.reference_start,
                                                read.reference_end or 
                                                    read.reference_start, 
                                                name))
        else:
            to_fetch.append(GenomicInterval(read.reference_name, 
                                            read.reference_start,
                                            read.reference_end or 
                                                read.reference_start, 
                                            name))
    logger.info("{:,} unpaired reads (potentially) over clipping "
                     .format(len(to_fetch)) + "threshold to fetch.")
    to_fetch = _merge_intervals(to_fetch, 200)
    got = 0
    clp = 0
    for interval in to_fetch:
        logger.debug("Fetching reads overlapping " + str(interval))
        for read in bamfile.fetch(region=str(interval)):
            read_name = read.query_name 
            p = 1 #read no. of pair
            if read.is_read1:
                p = 2
            if read_name in interval.targets:
                got += 1
                interval.targets.remove(read_name)
                if read_name in clipped_names:
                    fq_writer.output_pair(read, pair_dict[p][read_name])
                    clipped_names.remove(read_name)
                    clp += 1
                elif check_read_clipping(read):
                    fq_writer.output_pair(read, pair_dict[p][read_name])
                    clp += 1
    logger.info("Found {} previously unpaired mates, {} over clipping"
                     .format(got, clp) + " threshold.")

def _process_runner(tup):
    kwargs1, kwargs2 = tup
    kwargs2.update(kwargs1)
    return process_reads(**kwargs2)

def _initialize_mp_logger(logger, loglevel, logfile=None):
    logger.setLevel(loglevel)
    ch = logging.StreamHandler()  
    ch.setLevel(logging.INFO) 
    formatter = logging.Formatter(
           '[%(asctime)s] BACTA-%(processName)s - %(levelname)s - %(message)s') 
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    if logfile is not None:
        fh = logging.FileHandler(logfile)
        fh.setLevel(logger.level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

def process_reads(bam, tmp_fq1, tmp_fq2, tmp_bam, min_frac, min_clip=None,
                   ignore_dups=False, no_caching=False, region=None, 
                   contig=None, decoys=None, get_unaligned=False, paired=None,
                   unmapped_only=False, logfile=None, loglevel=logging.INFO):
    ''' 
        Read BAM and output reads reaching clip/mismatch threshold
        to self.fq1 and self.fq2. Returns no. reads which have been
        cached to self.bam_cache for processing of pairs mapped to
        different contigs.
    '''
    candidate_qnames = set()
    pair_tracker = defaultdict(dict)
    n = 0
    cached = 0
    cigar_scorer = CigarScorer(loglevel)
    logger = mp.get_logger()
    if logger.level == logging.NOTSET:
        _initialize_mp_logger(logger, loglevel, logfile)
    bf = get_bamfile(bam) #input bamfile
    cache_bamfile = get_align_output(tmp_bam, template=bf)
    fq_writer = Read2Fastq(tmp_fq1, tmp_fq2, decoys)
    kwargs = {}
    target_loci = None
    prev_chrom = None
    if unmapped_only:
        seek_back_til_reads(bf)
        kwargs['until_eof'] = True
    elif region:
        kwargs['region'] = str(region)
        if region.targets:
            target_loci = list(sorted(region.targets))
    elif contig:
        kwargs['contig'] = contig
    else:
        kwargs['until_eof'] = True
    for read in bf.fetch(**kwargs):
        n += 1
        if not n % 1000000:
            coord = _get_read_coord(read)
            logger.info("Reading input: {:,} records read. At pos {}"
                             .format(n, coord))
            logger.debug("Tracking {:,} pairs in RAM, {:,} reads cached "
                              .format((len(pair_tracker[1]) + 
                                      len(pair_tracker[2])), cached) + 
                              "to disk.")
        if read.is_secondary or read.is_supplementary:
            continue
        #read_name = self.parse_read_name(read.query_name) 
        read_name = read.query_name
        p = 2 #read no. of pair
        r = 1 #read no. of this read
        if read.is_read2:
            p = 1
            r = 2
        # do any modern aligners still keep the pair/tag on read ID?
        if ignore_dups and read.is_duplicate: 
            if read_name in pair_tracker[p]:
                # if mate is unmapped, mate won't be flagged as dup
                del pair_tracker[p][read_name]
                # read might be in candidate_qnames due to mate cigar
                if read_name in candidate_qnames:
                    candidate_qnames.remove(read_name)
            continue
        if target_loci is not None:
            if not _overlaps_loci(read, target_loci):
                continue
        if read.reference_id != prev_chrom and not no_caching: 
            # clearing pair_tracker will free up RAM in the case of 
            # malformed BAM files such as those processed with some 
            # versions of bcbio containing duplicated unmapped reads with 
            # mapped mates. For well formed BAMs pair_tracker and 
            # candidate_qnames should already be empty.
            if pair_tracker[1] or pair_tracker[2]:
                contig = '*'
                try:
                    contig = bf.get_reference_name(prev_chrom)
                except ValueError:
                    pass
                logger.warn("Clearing {} unpaired reads at end of "  
                            .format((len(pair_tracker[1]) + 
                                     len(pair_tracker[2]))) + 
                            "contig " + contig)
                pair_tracker.clear()
                if candidate_qnames:
                    logger.warn("Clearing {} unmatched clipped reads "
                                     .format(len(candidate_qnames)) + 
                                    "at end of contig " + contig)
                    candidate_qnames.clear()
            # use refernce_id not reference_name to prevent ValueError for 
            # unmapped reads at end of file
            prev_chrom = read.reference_id 
        if read.is_paired:
            if paired is None:
                paired = True
            elif not paired:
                raise RuntimeError('Mixed paired/unpaired reads can not ' +
                                   'be handled by this program yet. ' + 
                                   'Exiting.')
            if (read_name in candidate_qnames and 
                read_name in pair_tracker[p]):
                # checking pair_tracker[p] allows for problematic 
                # bams where unmapped reads are duplicated
                fq_writer.output_pair(read, pair_tracker[p][read_name])
                candidate_qnames.remove(read_name)
                del pair_tracker[p][read_name]
            else:
                store, is_clipped = should_store_is_clipped(read, cigar_scorer,
                                                   min_frac, min_clip=min_clip, 
                                                   decoys=decoys, 
                                                   get_unaligned=get_unaligned)
                if (not no_caching and 
                   read.reference_id != read.next_reference_id):
                    if store:
                        cache_bamfile.write(read)
                        cached += 1
                elif is_clipped:
                    if read_name in pair_tracker[p]:
                        fq_writer.output_pair(read, pair_tracker[p][read_name])
                        del pair_tracker[p][read_name]
                    else:
                        candidate_qnames.add(read_name)
                        pair_tracker[r][read_name] = read
                elif read_name in pair_tracker[p]:
                    del pair_tracker[p][read_name]
                else:#first encountered of pair
                    pair_tracker[r][read_name] = read
        else: #single-end reads
            if paired is None:
                paired = False
            elif paired:
                coord = _get_read_coord(read)
                logger.warn("Skipping unpaired read '{}' at {}"
                                 .format(read_name, coord))
                continue
            if check_read_clipping(read):
                fq_writer.output_single(read)
    if paired is None:
        paired = False
    f_msg = "Finished reading {:,} reads from input BAM".format(n)
    if unmapped_only:
        f_msg += " for unmapped reads"
    elif region:
        f_msg += " for region{}".format(region) 
    elif contig is not None:
        f_msg += " for contig {}".format(contig) 
    logger.info(f_msg)
    if region is not None and pair_tracker:
        _mop_up_region_pairs(bf, pair_tracker, candidate_qnames, fq_writer, 
                             logger)
        if candidate_qnames:
            logger.warning("{} unpaired candidate reads remaining "
                                .format(len(candidate_qnames)) + "after " +
                                "parsing region {}".format(region))
    bf.close()
    cache_bamfile.close()
    fq_writer.close()
    return {"n_cached": cached, "bam_cache": tmp_bam, "fq1": fq_writer.fq1, 
            "fq2": fq_writer.fq2, "paired": paired}


class BamAnalyzer(object):

    def __init__(self, bam, ref, output=None, contaminants=None, 
                bwa=None, samtools=None, min_fraction_clipped=0.15, 
                min_bases_clipped=None, min_expect_diff=1000, 
                min_aligned_score=50, fastqs=None, tmp=None, paired=None, 
                regions=[], vcf=None, flanks=500, ignore_dups=False, 
                threads=1, quiet=False, debug=False, log_file=None,
                no_caching=False, unaligned=False, decoy_contigs=[]):
        '''
            Read and identify potentially contaminating reads in a BAM 
            file according to read clipping and a reference fasta file.

            Args:
                bam:    Input BAM file.

                ref:    Reference fasta file for potentially 
                        contaminating sequences.

                output: Cleaned BAM output filename. 

                contaminants:
                        Prefix for contaminant output files. Defaults to 
                        the basename of the input + "_contaminants".
                
                regions: 
                        List of regions (in format chr1:1-1000) to scan
                        rather than processing the whole BAM file.

                vcf:    VCF file of candidate variants introduced by 
                        contamination. If provided the BAM will be 
                        scanned for reads that overlap these reads 
                        + flanks rather than processing the whole BAM
                        file.

                flanks: Amount of bp up and downstream of candidate 
                        variants (from vcf argument) to scan for reads.
                        Default=500.

                bwa:    Location of bwa executable. Only required if not 
                        in your PATH.

                samtools:
                        Location of samtools executable. Only required 
                        if not in your PATH.

                min_fraction_clipped:
                        Minimum proportion of a read that is hard or 
                        soft-clipped or single nucleotide mismatch in the BAM 
                        file for a read to be analyzed as a potential 
                        contaminant. Default=0.15.

                min_bases_clipped: 
                        Minimum number of hard or soft-clipped bases for 
                        a read to be analyzed as a potential 
                        contaminant. This overrides the 
                        --min_fraction_clipped argument.

                min_expect_diff:
                        TODO
                        Default=1000

                min_aligned_score:
                        TODO
                        Default=50 

                fastqs: Prefix for fastq files created from reads that
                        are tested for contamination. By default 
                        temporary files are created and deleted after
                        the program exits. If this argument is given,
                        these fastq files will be named 
                        <prefix>_r1.fastq.gz and <prefix>_r2.fastq.gz.

                tmp:    Directory to use for temporary files. If not 
                        specified, the system default temporary 
                        directory will be used.

                paired: Expect paired end reads if True, single end 
                        reads if False. If None, then will work it out 
                        on the fly. 

                no_caching:
                        Do not cache any reads to disk and hold all as
                        yet unpaired reads in memory. By default, reads
                        with a mate mapped to another chomosome are
                        cached to disk and processed after reading the
                        input file. Set this option to True to disable
                        this caching at the expense of more RAM usage.

                unaligned:
                        Also test unmapped reads with unmapped mates (or 
                        simply unmapped reads if input is from single-end 
                        reads).

                decoy_contigs:
                        Names of decoy contigs from input BAM. If specified, 
                        any reads mapped to these contigs will be tested 
                        regardless of clipping. For scoring purposes, these 
                        reads will be treated as if unmapped.
        
                threads:
                        Number of threads to use with BWA alignment step.
                        Default=1.

        '''
        self.commandline = str.join(" ", sys.argv)
        self.bam = bam
        self.bamfile = get_bamfile(bam)
        self.bam_ref_length = sum(self.bamfile.lengths)
        if not os.path.isfile(ref):
            raise RuntimeError('--ref argument "{}" does not ' .format(ref) + 
                               'exist or is not a file')
        self.ref = ref
        self.ref_length = self._length_from_fai(ref)
        if contaminants is None:
            contaminants = (os.path.splitext(bam)[0] + "_contaminants")
        self.output = output
        self.c_bam = contaminants + '.bam'
        self.c_sum = contaminants + '_summary.txt'
        if (self.output is not None):
            if not self.output.endswith(('.bam', '.BAM', '.cram', '.CRAM',
                                         '.sam', '.SAM')):
                self.output += '.bam'
            if self.c_bam == self.output:
                sys.exit("ERROR: --contaminants bam output and cleaned bam "
                         "--output have the same name ('{}'). Please ensure "
                         .format(self.c_bam) + "output file names are unique.")
        self.decoys = set(decoy_contigs)
        self.get_unaligned = unaligned
        self.min_fraction_clipped = min_fraction_clipped
        self.min_bases_clipped = min_bases_clipped
        self.min_expect_diff = min_expect_diff
        self.min_aligned_score  = min_aligned_score 
        self.tmp = tmp
        self.tmp_files = []
        self.no_caching = no_caching
        self.fastqs = fastqs
        self.bam_cache = None
        self.paired = paired
        self.ignore_dups = ignore_dups 
        if fastqs is None:
            tmpdir =  tempfile.mkdtemp(prefix="bacta_fastq", dir=self.tmp)
            self.fq1 = os.path.join(tmpdir, 'r1.fq.gz')
            self.fq2 = os.path.join(tmpdir, 'r2.fq.gz')
        else:
            self.fq1 = self.fastqs + '_r1.fastq.gz'
            self.fq2 = self.fastqs + '_r2.fastq.gz'
        for f in [self.output, self.c_bam, self.c_sum, self.fq1, self.fq2]:
            if f is not None and not self._is_writable(f):
                raise RuntimeError('path {} is not writeable'.format(f))
        if bwa is None:
            bwa = "bwa"
        if samtools is None:
            samtools = "samtools"
        for prog in (bwa, samtools):
            if self._which(prog) is None:
                raise RuntimeError('Could not find {} executable'.format(prog)+
                                   ' - please ensure it is installed and on ' + 
                                   'your PATH or otherwise provide the ' + 
                                   'location of a valid executable using the' +
                                   'appropriate argument.')
        self.bwa = bwa
        self.samtools = samtools
        self.threads = threads
        self.log_file = log_file
        self.debug = debug
        self.quiet = quiet
        self._set_logger()
        self.cigar_scorer = CigarScorer(self.logger.level)
        self.targets = []
        self.flanks = flanks
        if vcf:
            self._parse_vcf_regions(vcf)
        if regions:
            self._parse_regions(regions)
        if self.targets:
            self._check_targets()

    def _length_from_fai(self, fasta):
        l = 0
        fai = fasta + '.fai'
        if not os.path.exists(fai):
            sys.exit("ERROR: could not find fasta index ('{}') for fasta "
                     .format(fai) + "reference. Please index ('samtools faidx"+
                     " {}') before running this program." .format(fasta))
        with open(fai, 'r') as fh:
            for line in fh:
                try:
                    l += int(line.split()[1])
                except ValueError:
                    sys.exit("ERROR: Fasta index '{}' appears to be malformed."
                             .format(fai) + " Could not determine sequence " + 
                             "length for line:\n{}".format(line))
        return l
        
            
    def _check_targets(self):
        ldict = dict(zip( self.bamfile.references, self.bamfile.lengths))
        for gi in self.targets:
            if gi.contig not in ldict:
                raise SystemExit("Contig {} not in BAM index" 
                                 .format(gi.contig))
            if ldict[gi.contig] < gi.start:
                raise SystemExit("ERROR: Start of region {} is ".format(gi) + 
                                 "greater than the length of contig {}."
                                 .format(gi.contig))
            if ldict[gi.contig] < gi.end:
                self.logger.warn("ERROR: End of region {} is ".format(gi) + 
                                 "greater than the length of contig {}. "
                                 .format(gi.contig) + "Reducing region " + 
                                 "length to {}:{}-{}".format(gi.contig, gi.start, 
                                                   ldict[gi.contig]))
                gi.end = ldict[gi.contig]
            

    def _parse_vcf_regions(self, vcf):
        self.logger.info("Parsing VCF input and identifying regions...")
        vfile = pysam.VariantFile(vcf)
        seen_contigs = set()
        region = None
        prev_var = None
        var_count = 0
        for rec in vfile.fetch():
            # create intervals based on variant pos plus self.flanks, merging
            # any overlapping intervals
            var_count += 1
            start = rec.pos - self.flanks
            end = rec.pos + self.flanks
            if region is None:
                region = GenomicInterval(rec.contig, start, end)
            elif region.contig != rec.contig:
                if rec.contig in seen_contigs:
                    raise SystemExit("VCF input must be sorted by chromosome" + 
                                     " and position. Contig '{}' encountered " 
                                     .format(record.contig) + "before and " + 
                                     "after contig '{}'." .format(region[0]))
                self.targets.append(region)
                region = GenomicInterval(rec.contig, start, end, rec.pos)
                seen_contigs.add(rec.contig)
            else:
                if rec.pos < prev_var.pos:
                    raise SystemExit("VCF input must be sorted by chromosome" + 
                                     " and position. Position {}:{} "
                                     .format(rec.contig, rec.pos) + 
                                     "encountered after position {}:{}."
                                     .format(prev_var.contig, prev_var.pos))
                if start <= region.end: 
                    region.end = end #flanks are constand and VCF is sorted
                                     #so end must be > than region.end
                else: #no overlap
                    self.targets.append(region)
                    region = GenomicInterval(rec.contig, start, end)
            prev_var = rec
        self.targets.append(region)
        self.logger.info("Identified {:,} non-overlapping ({:,} bp) regions " 
                         .format(len(self.targets), self._get_targets_length()) 
                         +"from {:,} variants in VCF input.".format(var_count))

    def _get_targets_length(self):
        ''' 
            Calculates length of target intervals, which should be 
            non-overlapping.
        '''
        length = 0;
        for t in self.targets:
            length+= 1 + t.end - t.start
        return length

    def _parse_regions(self, regions):
        self.targets.extend(self._region_strings_to_lists(regions))
        self.targets = self._merge_intervals(self.targets)
        self.logger.info("Processing {:,} total non-overlapping regions"
                         .format(len(self.targets)) + " ({:,} bp)." 
                         .format(self._get_targets_length()))

    def _region_strings_to_lists(self, regions):
        reg_splitter = re.compile(r'[:\-]')
        gis = []
        for reg in regions:
            split = reg_splitter.split(reg)
            if len(split) != 3:
                raise SystemExit("Invalid region specified: '{}'" .format(reg))
            try:
                start = int(split[1].replace(',', ''))
                end = int(split[2].replace(',', ''))
                gis.append(GenomicInterval(split[0], start, end))
            except ValueError as oops:
                raise SystemExit("ERROR: Invalid region specified: '{}'"
                                 .format(reg) + "\n" + str(oops))
        return gis

    def _set_logger(self):
        self.logger = logging.getLogger("BACTA")
        if self.debug:
            self.logger.setLevel(logging.DEBUG)
        elif self.quiet:
            self.logger.setLevel(logging.WARNING)
        else:
            self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter(
                        '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
        ch = logging.StreamHandler()
        ch.setLevel(self.logger.level)
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)
        if self.log_file is not None:
            fh = logging.FileHandler(self.log_file)
            fh.setLevel(self.logger.level)
            fh.setFormatter(formatter)
            self.logger.addHandler(fh)

    def _is_writable(self, file_path):
        if not os.path.dirname(file_path):
            file_path = os.path.join(".", file_path)
        if os.access(os.path.dirname(file_path), os.W_OK):
            return True
        return False
    
    def _which(self, prog):
        path, name = os.path.split(prog)
        if path:
            if os.path.isfile(prog) and os.access(prog, os.X_OK):
                return prog
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, prog)
                if os.path.isfile(exe_file) and os.access(exe_file, os.X_OK):
                    return exe_file
        return None

    def cleanup(self):
        for f in self.tmp_files:
            if os.path.exists(f):
                os.remove(f)
        if self.fastqs is None: #cleanup tmp fastqs
           shutil.rmtree(os.path.split(self.fq1)[0])
        if self.bam_cache:
            fn = self.bam_cache.filename
            collated = os.path.splitext(fn)[0] + "_collated" + '.bam'
            for f in (fn, collated):
                if os.path.exists(f):
                    os.remove(f)

    def clean_bam(self):
        contam_reads = set()
        # get read IDs for contaminants from contaminant summary file created 
        # with  self.align_candidates
        with open(self.c_sum, 'rt') as cfile:
            for line in cfile:
                if line[0] == '#':
                    continue
                contam_reads.add(line.split()[0])
        if not contam_reads:
            self.logger.warn("No candidate reads for cleaning - skipping bam" +
                             " cleaning.")
            return
        if self.bamfile.is_open(): 
            #want to make sure we reset pointer to beggining of file
            self.bamfile.close()
        self.bamfile = get_bamfile(self.bam)
        header = self.bamfile.header
        pgid = self._get_pg_id(header)
        header['PG'].append({'ID' : pgid, 'PN' : 'bacta.py', 
                             'CL' : self.commandline, 'VN' : __version__})
        outbam = get_align_output(self.output, header=header)
        n = 0
        f = 0
        for read in self.bamfile.fetch(until_eof=True):
            n += 1
            if not n % 100000:
                coord = _get_read_coord(read)
                self.logger.info("Cleaning bam: {:,} records read, {:,} "
                                 .format(n, f) + "records filtered. At pos " +
                                 coord)
            if read.query_name in contam_reads:
                f += 1
            else:
                outbam.write(read)
        self.bamfile.close()
        outbam.close()
        self.logger.info("Finished cleaning bam. {:,} records read, {:,} "
                         .format(n, f) + "filtered")
    
    def _get_pg_id(self, header):
        ''' Ensure @PG ID is unique '''
        pg = "BACTA"
        prog_ids = set(x['ID'] for x in header['PG'])
        while 1:
            if pg in prog_ids:
                if pg == "BACTA":
                    pg = "BACTA.1"
                else:
                    b,n = pg.split(".")
                    n = int(n) + 1
                    pg = "BACTA." + str(n)
            else:
                return pg

    def align_candidates(self):
        ''' use bwa to align candidate reads (after running read_bam).'''
        #bam_out = open as pipe to BAM self.c_bam 
        contam_out = open(self.c_sum, 'w')
        bam_out = open(self.c_bam, 'wb')
        args = [self.bwa, 'mem', '-C', self.ref, self.fq1]
        if self.paired:
            args.append(self.fq2)
            pairs = dict()
        if self.threads > 1:
            args += ["-t", str(self.threads)]
        prev_cigar_re = re.compile(r'ZC:Z:(\S+)')
        prev_pos_re = re.compile(r'ZP:Z:(\S+)')
        prev_md_re = re.compile(r'ZM:Z:(\S+)')
        prev_mapq_re = re.compile(r'ZQ:Z:(\S+)')
        md_re = re.compile(r'MD:Z:(\S+)')
        self.logger.info("Attempting alignment of clipped reads against {} " 
                         .format(self.ref) + "using bwa.")
        self.logger.info("Command is: " + str.join(" ", args))
        self.logger.info("Piping output to " + self.c_bam + " via samtools.")
        contam_out.write(str.join("\t", ("#ID", "SCORE", "OLDSCORE", "EXPECT",
                                         "OLDEXPECT", "CIGAR", "OLDCIGAR", 
                                         "OLDPOS", "CONTIG", "POS", 
                                         "MATE_CONTIG", "MATE_POS", "MAPQ",
                                         "OLD_MAPQ")) + "\n")
        dash = '-' #subprocess doesn't seem to like '-' passed as a string
        samwrite = Popen([self.samtools, 'view', '-Sbh', dash], stdout=bam_out, 
                        stdin=PIPE, bufsize=-1)
        bwamem = Popen(args, bufsize=-1, stdout=PIPE)
        nparsed = 0
        for alignment in bwamem.stdout:
            samwrite.stdin.write(alignment)
            alignment = alignment.decode(sys.stdout.encoding)
            if alignment[0] == '@': #header
                continue
            nparsed += 1
            split = alignment.split()
            flag = int(split[1])
            if flag & 256 or flag & 2048: #not primary/supplementary 
                continue
            old_cigar = ''
            old_pos = ''
            old_md = ''
            old_mapq = ''
            md = ''
            for tag in split[:10:-1]:
                match = md_re.match(tag)
                if match:
                    md = match.group(1)
                    break #MD should be before all custom tags - bail
                match = prev_cigar_re.match(tag)
                if match:
                    old_cigar = match.group(1)
                    continue
                match = prev_pos_re.match(tag)
                if match:
                    old_pos = match.group(1)
                    continue
                match = prev_mapq_re.match(tag)
                if match:
                    old_mapq = match.group(1)
                    continue
                match = prev_md_re.match(tag)
                if match:
                    old_md = match.group(1)
                    continue
            score = self.cigar_scorer.score_cigarstring(split[5])
            old_score = self.cigar_scorer.score_cigarstring(old_cigar)
            if old_md and md:
                score -= self.cigar_scorer.md_mismatches(md)
                old_score -= self.cigar_scorer.md_mismatches(old_md)
            read_length = len(split[9])
            e = self._calc_expect(score, read_length, self.ref_length)
            old_e = self._calc_expect(old_score, read_length, 
                                      self.bam_ref_length)
            mapq = split[4]
            if self.paired:
                if split[0] in pairs:
                    (p_split, p_score, p_old_score, p_e, p_old_e, p_old_cigar, 
                     p_old_pos, p_old_mapq) = pairs[split[0]]
                    if (score + p_score >= self.min_aligned_score *2 and
                        (old_e * p_old_e)/(e * p_e) >= self.min_expect_diff):
                        self.write_contam_summary(contam_out, score, 
                                                  old_score, e, old_e,
                                                  split, old_cigar, old_pos, 
                                                  old_mapq)
                        self.write_contam_summary(contam_out, p_score, 
                                                  p_old_score, p_e, p_old_e, 
                                                  p_split, p_old_cigar, 
                                                  p_old_pos, p_old_mapq) 
                    del pairs[split[0]]
                else:
                    pairs[split[0]] = (split, score, old_score, e, old_e, 
                                       old_cigar, old_pos, old_mapq)
            else:
                if (score >= self.min_aligned_score and
                    old_e/e >= self.min_expect_diff):
                    self.write_contam_summary(contam_out, score, old_score, e, 
                                              old_e, split, old_cigar, old_pos,
                                              old_mapq)
        contam_out.close()
        samwrite.stdin.close()
        bwamem.stdout.close()
        bexit = bwamem.wait()#timeout not compatible with python2
        sexit = samwrite.wait()
        if bexit > 0:
            sys.exit("ERROR: bwa command failed with exit code {}"
                     .format(bexit))
        if sexit > 0:
            sys.exit("ERROR: samtools command failed with exit code {}"
                     .format(sexit))
        self.logger.info("Finished bwa mem run - {} reads parsed. "
                         .format(nparsed)) 
        self.logger.info("Alignments of all reads reaching the clip/mismatch" + 
                         " threshold are in {}.".format(self.c_bam))
        self.logger.info("Putative contaminant read information is in {}."
                         .format(self.c_sum))

    def _calc_expect(self, score, read_length, ref_length):
        ''' 
            Crude approximation of E-value which should be acceptable
            when comparing identical sequences and long references.
        '''
        return read_length * ref_length * 2**-score

    def write_contam_summary(self, fh, score, old_score, e, old_e, record, 
                             old_cigar, old_pos, old_mapq):
        ''' Write a summary of a contaminant read to fh. '''
        if int(record[1]) & 2: #pair_mapped
            mate_coord = (record[6] if record[6] != '=' else record[2], 
                          record[7])
        else:
            mate_coord = ('.', '.')
        #format is: ReadID, Score, OldScore, CIGAR, OldCigar, 
        #           OldPos, Chrom, Pos MateChrom, MatePos, MAPQ, OldMAPQ
        fh.write(str.join("\t", (record[0], str(score), str(old_score), 
                                 str.format("{:g}", e), 
                                 str.format("{:g}", old_e),
                                 record[5], old_cigar, old_pos, record[2], 
                                 record[3], mate_coord[0], mate_coord[1],
                                 record[4], old_mapq))
                         + "\n")

    def read_bam(self):
        ''' 
            Read a BAM/SAM/CRAM file and identify reads with excessive
            clipping.
        '''
        results = []
        kwargs = {"bam": self.bam, 
                  "min_frac": self.min_fraction_clipped, 
                  "min_clip": self.min_bases_clipped, 
                  "ignore_dups": self.ignore_dups,
                  "no_caching": self.no_caching, 
                  "decoys": self.decoys,
                  "get_unaligned": self.get_unaligned}
        
        if self.threads > 1:
            with mp.Pool(self.threads) as p:
                if self.targets:
                    rargs = ({'region': x} for x in self.targets)
                    rargs = self._add_process_args(rargs)
                    results = p.map(_process_runner, 
                                    zip(rargs, repeat(kwargs)), 1)
                else:
                    rargs = list({'contig': x} for x in self.bamfile.references)
                    if self.get_unaligned:
                        rargs.append({"unmapped_only": True})
                    rargs = self._add_process_args(rargs)
                    results = p.map(_process_runner, 
                                    zip(rargs, repeat(kwargs)), 1)
        else:
            (f, tmp_bam) = tempfile.mkstemp(suffix='.bam', dir=self.tmp)
            os.close(f)
            kwargs['tmp_bam'] = tmp_bam
            self.tmp_files.append(tmp_bam)
            (f, tmp_fq1) = tempfile.mkstemp(suffix='.r1.fq.gz', dir=self.tmp)
            os.close(f)
            kwargs['tmp_fq1'] = tmp_fq1
            self.tmp_files.append(tmp_fq1)
            (f, tmp_fq2) = tempfile.mkstemp(suffix='.r2.fq.gz', dir=self.tmp)
            os.close(f)
            kwargs['tmp_fq2'] = tmp_fq2
            self.tmp_files.append(tmp_fq2)
            if self.targets:
                for reg in self.targets:
                    results.append(process_reads(region=reg, **kwargs))
            else:
                results.append(process_reads(**kwargs))
        self.bamfile.close()
        #return {"n_cached": cached, "bam_cache": cache_fn, "fq1": fq_writer.fq1, 
        #        "fq2": fq_writer.fq2, "paired": paired}
        #self.bam_cache.close()
        n_cached = sum(x['n_cached'] for x in results)
        paired = sum(x['paired'] for x in results)
        fq1s = list(x['fq1'] for x in results)
        fq2s = list(x['fq2'] for x in results)
        cached_bams = list(x['bam_cache'] for x in results)
        if n_cached:
            cf1, cf2 = self._process_bam_cache(cached_bams) #deletes bam_cache after collating
            fq1s.append(cf1)
            fq2s.append(cf2)
        else:
            for f in cached_bams:
                if os.path.exists(f):
                    os.remove(f)
        self._concat_fqs(fq1s, fq2s, paired) #deletes all tmp fastqs after cat

    def _add_process_args(self, args):
        for a in args:
            (tmpfh1, tmp_bam) = tempfile.mkstemp(suffix='.bam', dir=self.tmp)
            os.close(tmpfh1)
            a['tmp_bam'] = tmp_bam
            self.tmp_files.append(tmp_bam)
            (tmpfh1, tmp_fq1) = tempfile.mkstemp(suffix='.r1.fq.gz', 
                                                 dir=self.tmp)
            os.close(tmpfh1)
            a['tmp_fq1'] = tmp_fq1
            self.tmp_files.append(tmp_fq1)
            (tmpfh2, tmp_fq2) = tempfile.mkstemp(suffix='.r2.fq.gz', 
                                                 dir=self.tmp)
            os.close(tmpfh2)
            a['tmp_fq2'] = tmp_fq2
            self.tmp_files.append(tmp_fq2)
        return args

    def _concat_fqs(self, r1s, r2s, paired):
        if paired:
            if len(r1s) != len(r2s):
                sys.exit("ERROR: Expected the same number of read 1 and " + 
                         "read 2 fastq files after reading input. Instad got"+
                         " {} and {} respectively".format(len(r1s), len(r2s)))
        self.logger.info("Concatanating read 1 fastqs")
        if len(r1s) == 1:
            os.rename(r1s[0], self.fq1)
            if paired:
                os.rename(r2s[0], self.fq2)
            return
        call("cat " + str.join(" ", r1s) + " > " + self.fq1, shell=True)
        self.logger.info("Removing temporary read 1 fastqs")
        for f in r1s:
            os.remove(f)
        if paired:
            self.logger.info("Concatanating read 2 fastqs")
            call("cat " + str.join(" ", r2s) + "> " + self.fq2, shell=True)
            self.logger.info("Removing temporary read 2 fastqs")
            for f in r2s:
                os.remove(f)

    def collate_single_bam(self, bam):
        collated = os.path.splitext(bam)[0] + "_collated"
        self.logger.info("Collating cached reads by read ID...")
        args = ['--output-fmt', 'BAM', bam, collated]
        if self.threads > 1:
            args += ['-@', str(self.threads - 1)]
        pysam.collate(*args)
        self.logger.info("Finished collating cached reads.")
        os.remove(bam)
        return pysam.AlignmentFile(collated + '.bam', 'rb')
 
    def collate_multiple_bams(self, bams):
        self.logger.info("Concatanating cached reads and collating by read " + 
                         " ID using samtools")
        dash = '-'
        (f, collated) = tempfile.mkstemp(suffix='.collated', dir=self.tmp)
        os.close(f)
        os.remove(collated) #samtools will create & add .bam extension
        cat_args = [self.samtools, 'cat'] + bams
        col_args = [self.samtools, 'collate', '-@', str(self.threads - 1), 
                    dash, collated]
        self.logger.debug("Concatanating command: " + str.join(" ", cat_args))
        self.logger.debug("Collating command: " + str.join(" ", col_args))
        bamcat = Popen(cat_args, stdout=PIPE, bufsize=-1)
        bamcollate = Popen(col_args, stdin=PIPE)
        for buff in bamcat.stdout:
            bamcollate.stdin.write(buff)
        bamcat.stdout.close()
        bamcollate.stdin.close()
        cat_exit = bamcat.wait()
        col_exit = bamcollate.wait()
        if cat_exit > 0:
            sys.exit("ERROR: samtools cat command failed with exit code {}"
                     .format(cat_exit))
        if col_exit > 0:
            sys.exit("ERROR: samtools collate command failed with exit code {}"
                     .format(col_exit))
        for f in bams:
            if os.path.exists(f):
                os.remove(f)
        return pysam.AlignmentFile(collated + '.bam', 'rb')

    def _process_bam_cache(self, cache_bams):
        ''' 
            Identify reads with excessive clipping in self.bam_cache.
            self.bam_cache should only contain paired reads where each 
            mate is on a different contig. Have already skipped 
            supplementary/secondary reads.
        '''
        if len(cache_bams) == 1:
            tmp_bam = self.collate_single_bam(cache_bams[0])
        else:
            tmp_bam = self.collate_multiple_bams(cache_bams)
        candidate_qnames = set()
        pair_tracker = defaultdict(dict)
        n = 0
        self.logger.info("Processing cached reads...")
        (f, tmp_fq1) = tempfile.mkstemp(suffix='.r1.fq.gz', dir=self.tmp)
        os.close(f)
        (f, tmp_fq2) = tempfile.mkstemp(suffix='.r2.fq.gz', dir=self.tmp)
        os.close(f)
        fq_writer = Read2Fastq(tmp_fq1, tmp_fq2, self.decoys)
        for read in tmp_bam.fetch(until_eof=True):
            n += 1
            if not n % 100000:
                coord = _get_read_coord(read)
                self.logger.info("Reading cache: {:,} records read. At pos {}"
                                 .format(n, coord))
                self.logger.debug("Tracking {:,} cached pairs in RAM"
                                  .format((len(pair_tracker[1]) + 
                                           len(pair_tracker[2]))))
            read_name = read.query_name
            read_name = read.query_name
            p = 2 #read no. of pair
            r = 1 #read no. of this read
            if read.is_read2:
                p = 1
                r = 2
            if self.ignore_dups and read.is_duplicate: 
                if read_name in pair_tracker[p]:
                    # if mate is unmapped, mate won't be flagged as dup
                    del pair_tracker[p][read_name]
                    # read might be in candidate_qnames due to mate cigar
                    if read_name in candidate_qnames:
                        candidate_qnames.remove(read_name)
                continue
            if (read_name in candidate_qnames and 
                read_name in pair_tracker[p]):
                fq_writer.output_pair(read, pair_tracker[p][read_name])
                candidate_qnames.remove(read_name)
                del pair_tracker[p][read_name]
            else:
                store, is_clipped = should_store_is_clipped(read, 
                                              self.cigar_scorer,
                                              self.min_fraction_clipped,
                                              min_clip=self.min_bases_clipped, 
                                              decoys=self.decoys, 
                                              get_unaligned=self.get_unaligned)
                if is_clipped:
                    if read_name in pair_tracker[p]:
                        fq_writer.output_pair(read, pair_tracker[p][read_name])
                        del pair_tracker[p][read_name]
                    else:
                        candidate_qnames.add(read_name)
                        pair_tracker[r][read_name] = read
                elif read_name in pair_tracker[p]:
                    del pair_tracker[p][read_name]
                else:#first encountered of pair
                    pair_tracker[r][read_name] = read
        tmp_bam.close()
        fq_writer.close()
        os.remove(tmp_bam.filename)
        return(tmp_fq1, tmp_fq2)
   
    def score_read(self, read):
        score = 0
        if read.cigartuples is not None: #check if mapped instead?
            score = cigar_scorer.score_cigartuples(read.cigartuples)
            self.logger.debug("Read {}: cigar = {}, score = {}" 
                             .format(read_name, read.cigarstring, score)) 
        return score

class Read2Fastq(object):
    ''' Write alignments to FASTQ.'''

    def __init__(self, fq1, fq2, decoys):
        self.fq1 = fq1
        self.fq2 = fq2
        self.r1_fh = gzip.open(self.fq1, mode='wt')
        self.r2_fh = gzip.open(self.fq2, mode='wt')
        self.decoys = decoys

    def close(self):
        self.r1_fh.close()
        self.r2_fh.close()

    def read_to_fastq(self, read):
        ''' Get FASTQ string for a given BAM read.'''
        read_name = self.parse_read_name(read.query_name)
        coord = _get_read_coord(read, True)
        mapq = read.mapping_quality
        if read.reference_id < 0:
            cigar = '*'
        elif self.decoys and read.reference_name in self.decoys:
            cigar = '*'
        else:
            cigar = read.cigarstring or '*'
        header = ('@{} ZC:Z:{}\tZP:Z:{}\tZQ:Z:{}'.format(read_name, cigar, 
                                                         coord, mapq))
        if read.has_tag('MD'):
            header += "\tZM:Z:" + read.get_tag('MD')
        if read.is_reverse:
            seq = reverse_complement(read.seq)
            qual = read.qual[::-1]
        else:
            seq = read.seq
            qual = read.qual
        return str.join("\n", (header, seq, "+", qual)) + "\n"

    def output_single(self, read):
        rfq = self.read_to_fastq(read)
        self.r1_fh.write(rfq)

    def output_pair(self, read1, read2):
        r1fq = self.read_to_fastq(read1) 
        r2fq = self.read_to_fastq(read2) 
        if read1.is_read1:
            self.r1_fh.write(r1fq)
            self.r2_fh.write(r2fq)
        else:
            self.r1_fh.write(r2fq)
            self.r2_fh.write(r1fq)

    def parse_read_name(self, query_name):
        match = _qname_re.match(query_name)
        if match:
            return match.group(1)
        return query_name


class GenomicInterval(object):
    ''' Simple class for representing a genomic interval. '''

    __slots__ = ['contig', 'start', 'end', 'targets']

    def __init__(self, contig, start, end, target=None):
        ''' 
            Args:
                contig: Name of contig/chromosome (string)

                start:  Start coordinate (int).

                end:    End coordinate (int). A ValueError will be 
                        raised if this is smaller than start.

                target: One or multiple target positions within the 
                        start and end coordinates. These are stored as
                        a set in the targets property.

        '''
        self.contig = contig
        self.start = start
        self.end = end
        self.targets = set()
        if target is not None:
            try:
                self.targets.add(target)
            except TypeError:
                self.targets.update(target)
        if start > end:
            raise ValueError("GenomicInterval start can not be greater than " + 
                            "end (for interval {}:{}-{})"
                            .format(contig, start, end))
    
    def __str__(self):
        return self.contig + ':' + str(self.start) + '-' + str(self.end)

    def __eq__(self, other):
        return (self.contig == other.contig and self.start == other.start and 
               self.end == other.end)

    def __ne__(self, other):
        return (self.contig != other.contig or self.start != other.start or
               self.end != other.end)

    def __lt__(self, other):
        return (self.contig < other.contig or 
               (self.contig == other.contig and 
               (self.start < other.start or self.start == other.start and 
                self.end < other.end)))

    def __le__(self, other):
        return (self.contig < other.contig or 
               (self.contig == other.contig and self.end <= other.end))

    def __gt__(self, other):
        return (self.contig > other.contig or 
               (self.contig == other.contig and 
               (self.start > other.start or self.start == other.start and 
                self.end > other.end)))

    def __ge__(self, other):
        return (self.contig > other.contig or 
               (self.contig == other.contig and self.start >= other.start))



