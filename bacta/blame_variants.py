import sys
import os
import re
import logging
import gzip
import pysam
from .alignment_file_utils import get_bamfile
from parse_vcf import VcfReader

class BlameVariants(object):
    ''' Identify the number of contaminant reads versus
        non-contaminant carrying given variants.
    '''

    def __init__(self, contaminants, variants, bam, bed=None, output=None,
                 mapq=0, contam_ratio=0.0, quiet=False, debug=False):
        self.logger = self._get_logger(quiet, debug)
        self.mapq = mapq
        self.contam_ratio = contam_ratio
        self.read_contam_file(contaminants)
        self.vcf = VcfReader(variants)
        self.bamfile = get_bamfile(bam)
        self.bed = bed
        self.output = output
        self.outfile = None

    def read_contam_file(self, contamfile):
        self.contam_ids = set()
        self.contam_below_mapq = set()
        with open(contamfile, 'rt') as infile:
            header = next(infile).rstrip().split()
            indices = dict()
            for field in ('#ID', 'EXPECT', 'OLDEXPECT', 'OLD_MAPQ'):
                if field not in header:
                    raise RuntimeError("Could not find required header field" +
                                       " '{}' in {}".format(field, contamfile))
                i = header.index(field)
                field = field.replace('#', '')
                indices[field] = i
            for line in infile:
                cols = line.rstrip().split()
                if self.mapq and int(cols[indices['OLD_MAPQ']]) < self.mapq:
                    self.contam_below_mapq.add(cols[indices['ID']])
                else:
                    self.contam_ids.add(cols[indices['ID']])
        infile.close()

    def _get_logger(self, quiet=False, debug=False):
        logger = logging.getLogger("Blame Variants")
        if debug:
            logger.setLevel(logging.DEBUG)
        elif quiet:
            logger.setLevel(logging.WARNING)
        else:
            logger.setLevel(logging.INFO)
        formatter = logging.Formatter(
                        '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
        ch = logging.StreamHandler()
        ch.setLevel(logger.level)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        return logger

    def cleanup(self):
        if self.outfile:
            self.outfile.close()

    def assess_variants(self):
        if self.output is not None:
            self.outfile = open(self.output, 'wt')
        else:
            self.outfile = sys.stdout
        self.write_header()
        if self.bed is not None:
            self.parse_with_bed()
        else:
            for var in self.vcf.parser:
                self.check_variant(var)

    def write_header(self):
        self.vcf.header.add_header_field(
                name='BactaContamSupport',
                dictionary={'Number' : 'A',
                             'Type' : 'Integer',
                             'Description' : 'Number of contaminant reads ' +
                                             'supporting ALT allele according'+
                                             ' to BACTA'''},
                field_type='INFO')
        self.vcf.header.add_header_field(
                name='BactaNonContamSupport',
                dictionary={'Number' : 'A',
                             'Type' : 'Integer',
                             'Description' : 'Number of non-contaminant reads'+
                                             ' supporting ALT allele '+
                                             'according to BACTA'},
                field_type='INFO')
        self.outfile.write(str(self.vcf.header))

    def parse_with_bed(self):
        if self.bed.endswith((".gz", ".bgz")):
            bedfile = gzip.open(self.bed, mode='rt', errors='replace')
        else:
            bedfile = open(self.bed, 'rt')
        for line in bedfile:
            if line[0] == '#':
                continue
            cols = line.rstrip().split()
            self.vcf.set_region(cols[0], int(cols[1]), int(cols[2]))
            for var in self.vcf.parser:
                self.check_variant(var)
        bedfile.close()

    def check_variant(self, var):
        write_record = False
        i = 0
        contam = []
        non_contam = []
        for alt in var.DECOMPOSED_ALLELES:
            if len(alt.REF) != 1 or len(alt.ALT) != 1: #SNVs only
                continue
            pileup_iter = self.bamfile.pileup(alt.CHROM, alt.POS-1, alt.POS)
            supporting_contam = 0
            supporting_non_contam = 0
            for x in pileup_iter:
                #pysam positions are 0-based
                if x.reference_pos != var.POS - 1:
                    continue
                for p in x.pileups:
                    if p.is_del or p.is_refskip:
                        continue
                    if p.alignment.is_duplicate or p.alignment.is_secondary:
                        continue
                    seq = p.alignment.seq[p.query_position:p.query_position + 1]
                    if seq == alt.ALT:
                        query_name = p.alignment.query_name
                        if query_name in self.contam_ids:
                            supporting_contam += 1
                        elif query_name not in self.contam_below_mapq:
                            supporting_non_contam += 1
                break #can skip now we've hit our SNV position
            if supporting_contam:
                dp = supporting_non_contam + supporting_contam
                if supporting_contam/dp >= self.contam_ratio:
                    write_record=True
            contam.append(supporting_contam)
            non_contam.append(supporting_non_contam)
        if write_record:
            var.add_info_fields({"BactaContamSupport"   : contam,
                                 "BactaNonContamSupport": non_contam
                                })
            self.outfile.write(str(var) + "\n")
