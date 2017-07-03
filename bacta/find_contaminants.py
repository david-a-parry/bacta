import sys
import os
from subprocess import PIPE, Popen
import re
import pysam
import tempfile
import gzip
import shutil
import logging
from .cigar_scorer import CigarScorer
from .reverse_complement import reverse_complement

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

class BamAnalyzer(object):

    def __init__(self, bam, ref, output=None, contaminants=None, 
                bwa=None, samtools=None, min_fraction_clipped=0.2, 
                min_bases_clipped=None, min_score_diff=30, min_aligned_score=50,
                fastqs=None, tmp=None, max_pair_distance=1000000, paired=None, 
                regions=[], vcf=None, flanks=500, quiet=False, debug=False):
        '''
            Read and identify potentially contaminating reads in a BAM 
            file according to read clipping and a reference fasta file.

            Args:
                bam:    Input BAM file.

                ref:    Reference fasta file for potentially 
                        contaminating sequences.

                output: Output filename. Default=input + _bacta.bam.

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
                        soft-clipped in the BAM file for a read to be 
                        analyzed as a potential contaminant. Default=0.2.

                min_bases_clipped: 
                        Minimum number of hard or soft-clipped bases for 
                        a read to be analyzed as a potential 
                        contaminant. This overrides the 
                        --min_fraction_clipped argument.

                fastqs: Prefix for fastq files created from reads that
                        are tested for contamination. By default 
                        temporary files are created and deleted after
                        the program exits. If this argument is given,
                        these fastq files will be named 
                        <prefix>_r1.fastq.gz and <prefix>_r2.fastq.gz.

                tmp:    Directory to use for temporary files. If not 
                        specified, the system default temporary 
                        directory will be used.

                max_pair_distance:
                        For speed, reads are stored and analyzed with 
                        their mate unless creater than this distance 
                        apart from each other in the genome. Increase 
                        this value to favour speed at the expense of 
                        memory and decrease this value to favour memory
                        conservation over speed. Default=1000000.
                
                paired: Expect paired end reads if True, single end 
                        reads if False. If None, then will work it out 
                        on the fly.

        '''
        self.bam = bam
        self.bmode = 'rb'
        if bam.endswith(('.sam', '.SAM')):
            self.bmode = 'r'
        elif bam.endswith(('.cram', '.CRAM')):
            self.bmode = 'rc'
        self.bamfile = pysam.AlignmentFile(bam, self.bmode)
        if not os.path.isfile(ref):
            raise RuntimeError('--ref argument "{}" does not ' .format(ref) + 
                               'exist or is not a file')
        self.ref = ref
        if output is None:
            output = os.path.splitext(bam)[0] + "_bacta.bam"
        if contaminants is None:
            contaminants = (os.path.splitext(bam)[0] + "_contaminants")
        self.output = output
        self.c_bam = contaminants + '.bam'
        self.c_sum = contaminants + '_summary.txt'
        self.min_fraction_clipped = min_fraction_clipped
        self.min_bases_clipped = min_bases_clipped
        self.min_score_diff = min_score_diff
        self.min_aligned_score  = min_aligned_score 
        self.tmp = tmp
        self.fastqs = fastqs
        self.max_pair_distance = max_pair_distance
        self.paired = paired
        if fastqs is None:
            tmpdir =  tempfile.mkdtemp(prefix="bacta_fastq", dir=self.tmp)
            self.fq1 = os.path.join(tmpdir, 'r1.fq.gz')
            self.fq2 = os.path.join(tmpdir, 'r2.fq.gz')
        else:
            self.fq1 = self.fastqs + '_r1.fastq.gz'
            self.fq2 = self.fastqs + '_r2.fastq.gz'
        for f in [self.output, self.c_bam, self.c_sum, self.fq1, self.fq2]:
            if not self._is_writable(f):
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
        self.debug = debug
        self.quiet = quiet
        self._set_logger()
        self.cigar_scorer = CigarScorer(self.logger.level)
        self.r1_fq = gzip.open(self.fq1, mode='wt')
        self.r2_fq = gzip.open(self.fq2, mode='wt')
        self.targets = []
        self.flanks = flanks
        if vcf:
            self._parse_vcf_regions(vcf)
        if regions:
            self._parse_regions(regions)
        if self.targets:
            self._check_targets()

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
                region = GenomicInterval(rec.contig, start, end)
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
        self._merge_targets()
        self.logger.info("Processing {:,} total non-overlapping regions"
                         .format(len(self.targets)) + " ({:,} bp)." 
                         .format(self._get_targets_length()))

    def _merge_targets(self):
        merged = []
        self.targets = sorted(self.targets)
        prev_target = self.targets[0]
        for i in range(1, len(self.targets)):
            if prev_target.contig != self.targets[i].contig:
                merged.append(prev_target)
            elif self.targets[i].start <= prev_target.end:
                if prev_target.end < self.targets[i].end:
                    prev_target.end = self.targets[i].end
            else:
                merged.append(prev_target)
                prev_target = self.targets[i]
        merged.append(prev_target)
            
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
        if self.fastqs is None: #cleanup tmp fastqs
           shutil.rmtree(os.path.split(self.fq1)[0])

    def align_candidates(self):
        ''' use bwa to align candidate reads (after running read_bam).'''
        #bam_out = open as pipe to BAM self.c_bam 
        contam_out = open(self.c_sum, 'w')
        bam_out = open(self.c_bam, 'wb')
        args = [self.bwa, 'mem', '-C', self.ref, self.fq1]
        if self.paired:
            args.append(self.fq2)
            pairs = dict()
        prev_cigar_re = re.compile(r'ZC:Z:(\w+)')
        prev_pos_re = re.compile(r'ZP:Z:(\w+:\d+)')
        self.logger.info("Attempting alignment of clipped reads against {} " 
                         .format(self.ref) + "using bwa.")
        self.logger.info("Command is: " + str.join(" ", args))
        contam_out.write(str.join("\t", ("#ID", "SCORE", "OLDSCORE", "CIGAR",
                                         "OLDCIGAR", "OLDPOS", "CONTIG", "POS", 
                                         "MATE_CONTIG", "MATE_POS")) + "\n")
        dash = '-' #subprocess doesn't seem to like '-' passed as a string
        samwrite = Popen([self.samtools, 'view', '-Sbh', dash], stdout=bam_out, 
                        stdin=PIPE, bufsize=-1)
        with Popen(args, bufsize=-1, stdout=PIPE) as bwamem:
            for alignment in bwamem.stdout:
                samwrite.stdin.write(alignment)
                alignment = alignment.decode(sys.stdout.encoding)
                if alignment[0] == '@': #header
                    continue
                split = alignment.split()
                flag = int(split[1])
                if flag & 256 or flag & 2048: #not primary/supplementary 
                    continue
                old_cigar = ''
                old_pos = ''
                for tag in split[:10:-1]:
                    match = prev_cigar_re.match(tag)
                    if match:
                        old_cigar = match.group(1)
                        #tags were added with old cigar first - bail
                        break
                    match = prev_pos_re.match(tag)
                    if match:
                        old_pos = match.group(1)
                score = self.cigar_scorer.score_cigarstring(split[5])
                old_score = self.cigar_scorer.score_cigarstring(old_cigar)
                if self.paired:
                    if split[0] in pairs:
                        (p_split, p_score, p_old_score, p_old_cigar, 
                         p_old_pos) = pairs[split[0]]
                        if (score + p_score >= self.min_aligned_score *2 and
                            score + p_score > old_score + p_old_score + 
                            self.min_score_diff):
                            self.write_contam_summary(contam_out, score, 
                                                      old_score, split, 
                                                      old_cigar, old_pos) 
                            self.write_contam_summary(contam_out, p_score, 
                                                      p_old_score, p_split, 
                                                      p_old_cigar, p_old_pos) 
                        del pairs[split[0]]
                    else:
                        pairs[split[0]] = (split, score, old_score, old_cigar, 
                                           old_pos)
                else:
                    self.write_contam_summary(contam_out, score, old_score, 
                                              split, old_cigar, old_pos)
        samwrite.stdin.close()

    def write_contam_summary(self, fh, score, old_score, record, old_cigar, 
                             old_pos):
        ''' Write a summary of a contaminant read to fh. '''
        if int(record[1]) & 2: #pair_mapped
            mate_coord = (record[6] if record[6] != '=' else record[2], 
                          record[7])
        else:
            mate_coord = ('.', '.')
        #format is: ReadID, Score, OldScore, CIGAR, OldCigar, 
        #           OldPos, Chrom, Pos MateChrom, MatePos
        fh.write(str.join("\t", (record[0], str(score), str(old_score), 
                                 record[5], old_cigar, old_pos, record[2], 
                                 record[3], mate_coord[0], mate_coord[1]))
                         + "\n")

    def read_bam(self):
        ''' 
            Read a BAM/SAM/CRAM file and identify reads with excessive
            clipping.
        '''
        self._mate_fetched = set() #prevent fetching reciprocal mates
                                   # when using regions
        if self.targets:
            for region in self.targets:
                self.logger.info("Parsing region {}".format(region))
                self._process_reads(region)
        else:
            self._process_reads()
        del self._mate_fetched
        self.bamfile.close()
        self.r1_fq.close()
        self.r2_fq.close()

    def _process_reads(self, region=None):
        candidate_qnames = set()
        pair_tracker = dict()
        n = 0
        kwargs = {}
        if region is None:
            kwargs['until_eof'] = True
        else:
            kwargs['region'] = str(region)
        for read in self.bamfile.fetch(**kwargs):
            n += 1
            if not n % 10000:
                self.logger.info("{} records read. At pos {}:{}" .format(n, 
                                 read.reference_name, read.reference_start))
            if read.is_secondary or read.is_supplementary or read.is_duplicate: 
                continue
            #read_name = self.parse_read_name(read.query_name) 
            read_name = read.query_name
            # do any modern aligners still keep the pair/tag on read ID?
            if read.is_paired:
                if self.paired is None:
                    self.paired = True
                elif not self.paired:
                    raise RuntimeError('Mixed paired/unpaired reads can not ' +
                                       'be handled by this program yet. ' + 
                                       'Exiting.')
                if self._should_fetch_mate(read, read_name, region):
                    # we need to fetch mates if very far away to prevent
                    # storing too many reads while waiting to encounter their
                    # mates (but this is too slow to do for all reads) 
                    if read.has_tag('MC'):
                        # already got mate cigar string - check before fetching
                        # this is a standard tag, should not need to check type
                        mcigar = read.get_tag('MC')
                        if self.check_read_clipping(read, mate_cigar=mcigar):
                            try:
                                read2 = self.bamfile.mate(read)
                                self.output_pair(read, read2)
                                if region is not None:
                                    self._mate_fetched.add(read_name)
                            except ValueError:
                                self.logger.warn("Mate not found for {}"
                                                            .format(read_name))
                    else:
                        try:
                            read2 = self.bamfile.mate(read)
                            if region is not None:
                                self._mate_fetched.add(read_name)
                            if (self.check_read_clipping(read) or 
                                self.check_read_clipping(read2)):
                                #one of pair is clipped
                                self.output_pair(read, read2)
                        except ValueError:
                            self.logger.warn("Mate not found for {}"
                                                            .format(read_name))
                else:
                    if read_name in candidate_qnames:
                        self.output_pair(read, pair_tracker[read_name])
                        candidate_qnames.remove(read_name)
                        del pair_tracker[read_name]
                    elif self.check_read_clipping(read):
                        if read_name in pair_tracker:
                            self.output_pair(read, pair_tracker[read_name])
                            del pair_tracker[read_name]
                        else:
                            candidate_qnames.add(read_name)
                            pair_tracker[read_name] = read
                    elif read_name in pair_tracker:
                        del pair_tracker[read_name]
                    else:
                        pair_tracker[read_name] = read
                self.logger.debug("Tracking {} pairs."
                                  .format(len(pair_tracker)))
            else: #single-end reads
                if self.paired is None:
                    self.paired = False
                elif self.paired:
                    raise RuntimeError('Mixed paired/unpaired reads can not ' +
                                       'be handled by this program yet. ' + 
                                       'Exiting.')
                if self.check_read_clipping(read):
                    self.read_to_fastq(read, self.r1_fq)
        self.logger.info("Finished reading {} reads from input BAM" .format(n))
        if candidate_qnames:
            self.logger.warning("{} unpaired candidate reads remaining after "
                                .format(len(candidate_qnames)) + "finished " +
                                "parsing input BAM file.")

    def _should_fetch_mate(self, read, read_name, region):
        if region is not None and read_name in self._mate_fetched:
            return False
        if read.next_reference_id != read.reference_id:
            # if parsing til eof only fetch mates of first in pair so as not to
            # waste memory on self._mate_fetched
            return region is not None or read.is_read1
        elif (abs(read.next_reference_start - read.reference_start) > 
              self.max_pair_distance):
            #if parsing til eof only fetch mates of first in pair
            return region is not None or read.is_read1
        elif region is not None:
            if read.next_reference_start > region.end:
                return True
            #problem - we can't guarantee mate read length, so determining
            #if it overlaps region.start is impossible if it starts before it
            if read.next_reference_start < region.start:
                return True
        return False

    def score_read(self, read):
        score = 0
        if read.cigartuples is not None: #check if mapped instead?
            score = self.cigar_scorer.score_cigartuples(read.cigartuples)
            self.logger.debug("Read {}: cigar = {}, score = {}" 
                             .format(read_name, read.cigarstring, score)) 
        return score

    def check_read_clipping(self, read, mate_cigar=None):
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

                mate_cigar:
                        Optional cigar string for mate (e.g. as provided
                        by the standard 'MC' tag. If provided, this 
                        function will return True if either the read or
                        mate cigar indicate clipping above the threhold.
                
        '''
        if read.cigartuples is not None: #check if mapped instead?
            if self._check_cigar_tuple_clipping(read.cigartuples):
                return True
        if mate_cigar is not None:
            mcts = self.cigar_scorer.cigarstring_to_tuples(mate_cigar)
            return self._check_cigar_tuple_clipping(mcts)
        return False

    def _check_cigar_tuple_clipping(self, cts):
        ''' 
            Returns True if read is clipped greater than 
            min_fraction_clipped or min_bases_clipped.
        '''
        clipping = 0
        length = 0
        for c in cts:
            if c[0] >= 4 and c[0] <= 5:
                #SOFT or HARD clip
                clipping += c[1]
                length += c[1]
            elif c[0] < 2 or c[0] >= 7 and c[0] <= 8:
                #MATCH, INS, EQUAL or DIFF
                length += c[1]
        if clipping:
            if (self.min_bases_clipped is not None and 
                clipping >= self.min_bases_clipped):
                return True
            else:
                if float(clipping)/length >= self.min_fraction_clipped:
                    return True
        return False

    def read_to_fastq(self, read, fh):
        read_name = self.parse_read_name(read.query_name)
        header = ('@{} ZC:Z:{}\tZP:Z:{}:{}'.format(read_name, 
                                               (read.cigarstring or '.'),
                                               read.reference_name, 
                                               read.reference_start))
        if read.is_reverse:
            seq = reverse_complement(read.seq)
            qual = read.qual
        else:
            seq = read.seq
            qual = read.qual
        fh.write(str.join("\n", (header, seq, "+", qual)) + "\n")

    def output_pair(self, read1, read2):
        if read1.is_read1:
            self.read_to_fastq(read1, self.r1_fq)
            self.read_to_fastq(read2, self.r2_fq)
        else:
            self.read_to_fastq(read1, self.r2_fq)
            self.read_to_fastq(read2, self.r1_fq)

    def parse_read_name(self, query_name):
        match = _qname_re.match(query_name)
        if match:
            return match.group(1)
        return query_name


class GenomicInterval(object):
    ''' Simple class for representing a genomic interval. '''

    __slots__ = ['contig', 'start', 'end']

    def __init__(self, contig, start, end):
        ''' 
            Args:
                contig: Name of contig/chromosome (string)

                start:  Start coordinate (int).

                end:    End coordinate (int). A ValueError will be 
                        raised if this is smaller than start.
        '''
        self.contig = contig
        self.start = start
        self.end = end
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



