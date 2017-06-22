import sys
import os
import re
import pysam
import tempfile
import gzip
import shutil
import logging
from collections import defaultdict

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
_cigar_score = {
    0: lambda x: x, 
    1: lambda x: -6 - x, 
    2: lambda x: -6 - x,
    3: lambda x: 0, 
    4: lambda x: x * -1, 
    5: lambda x: x * -1, 
    6: lambda x: 0, 
    7: lambda x: x,
    8: lambda x: x * -1,
    9: lambda x: 0,
}

_qname_re = re.compile(r"""(\S+)(#\d+)?(/\[1-2])?""")  

class BamAnalyzer(object):

    def __init__(self, bam=None, ref=None, output=None, contaminants=None, 
                bwa=None, min_fraction_clipped=0.2, min_bases_clipped=None, 
                fastqs=None, tmp=None, max_pair_distance=1000000, 
                quiet=False, debug=False):
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
                
                bwa:    Location of bwa executable. Only required if not 
                        in your PATH.

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
        if output is None:
            output = os.path.splitext(bam)[0] + "_bacta.bam"
        if contaminants is None:
            contaminants = (os.path.splitext(bam)[0] + "_contaminants")
        self.output = output
        self.c_bam = contaminants + '.bam'
        self.c_sum = contaminants + '_summary.txt'
        self.min_fraction_clipped = min_fraction_clipped
        self.min_bases_clipped = min_bases_clipped
        self.tmp = tmp
        self.fastqs = fastqs
        self.max_pair_distance = max_pair_distance
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
        if self._which(bwa) is None:
            raise RuntimeError('Could not find bwa executable - please ' + 
                               'ensure bwa is installed and on your PATH or ' + 
                               'otherwise provide the location of a valid ' + 
                               'bwa executable using the --bwa argument')
        self.bwa = bwa
        self.debug = debug
        self.quiet = quiet
        self._set_logger()

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

    def read_bam(self):
        ''' 
            Read a BAM/SAM/CRAM file and identify reads with excessive
            clipping.
        '''
        paired = False
        self.r1_fq = gzip.open(self.fq1, mode='wt')
        self.r2_fq = gzip.open(self.fq2, mode='wt')
        self.cig_warned = set()
        candidate_qnames = set()
        pair_tracker = dict()
        n = 0
        for read in self.bamfile.fetch():
            n += 1
            if not n % 10000:
                self.logger.info("{} records read. At pos {}:{}" .format(n, 
                                 read.reference_name, read.reference_start))
            if read.is_secondary or read.is_supplementary or read.is_duplicate: 
                continue
            read_name = self.parse_read_name(read.query_name)
            if read.is_paired:
                paired = True
                if (read.next_reference_id != read.reference_id or 
                    abs(read.next_reference_start - read.reference_start) > 
                    self.max_pair_distance):
                    # we need to fetch mates if very far away to prevent
                    # storing too many reads while waiting to encounter their
                    # mates (but this is too slow to do for all reads) 
                    if read.is_read1: 
                        try:
                            read2 = self.bamfile.mate(read)
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
                if paired:
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
        self.bamfile.close()
            

    def score_read(self, read):
        score = 0
        if read.cigartuples is not None: #check if mapped instead?
            for c in read.cigartuples:
                try:
                    score += _cigar_score[c[0]](c[1])
                except KeyError:
                    if c[0] not in self.cig_warned:
                        self.logger.warn("Unrecognized operator code found " + 
                                         "in CIGAR STRING: '{}'".format(c[0]))
                        self.cig_warned.add(c[0])
            self.logger.debug("Read {}: cigar = {}, score = {}" 
                             .format(read_name, read.cigarstring, score)) 
        return score

    def check_read_clipping(self, read):
        ''' 
            Returns True if read is clipped greater than 
            min_fraction_clipped or min_bases_clipped.
        '''
        clipping = 0
        if read.cigartuples is not None: #check if mapped instead?
            for c in read.cigartuples:
                if c[0] == 4 or c[0] == 5:
                    clipping += c[1]
            if clipping:
                if (self.min_bases_clipped is not None and 
                    clipping >= self.min_bases_clipped):
                    return True
                else:
                    l = read.infer_read_length()
                    self.logger.debug("Read {}: cigar = {}, length = {},"
                                      .format(read.query_name, 
                                              read.cigarstring, 
                                      l) + " clip = {}" .format(clipping))
                    if float(clipping)/l >= self.min_fraction_clipped:
                        return True
        return False

    def read_to_fastq(self, read, fh):
        read_name = self.parse_read_name(read.query_name)
        header = '@' + read_name + " ZC:Z:" + (read.cigarstring or '.')
        fh.write(str.join("\n", (header, read.seq, "+", 
                          read.qual)) + "\n")

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
