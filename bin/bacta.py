#!/usr/bin/env python3
import sys
import os
import argparse
sys.path.insert(0, '')
from bacta.find_contaminants import BamAnalyzer

def get_parser():
    '''Get ArgumentParser'''
    parser = argparse.ArgumentParser(
                  description='Identify and remove bacterial reads from BAMS.')
    required_args = parser.add_argument_group('Required Arguments')
    optional_args = parser.add_argument_group('Optional Arguments')
    required_args.add_argument('-i', '--bam', '--input', required=True, 
                               metavar='BAM', help='''Input BAM filename''')
    required_args.add_argument('-r', '--ref', required=True, metavar='REF', 
                               help='''Reference fasta file containing 
                                        potential contaminant sequences. This 
                                        file  should be indexed with BWA.''')
    optional_args.add_argument('-o', '--output', metavar='BAM', 
                               help='''Write cleaned output BAM to this file. 
                                       If specified, a new BAM file without 
                                       contaminating reads will be created. 
                                    ''')
    optional_args.add_argument('-c', '--contaminants', metavar='PREFIX', 
                               help='''Prefix for contaminant output files. 
                                       Defaults to the basename of the input + 
                                       "_contaminants".
                                    ''')
    optional_args.add_argument('-v', '--vcf', metavar='VCF', 
                               help='''VCF file of candidate variants 
                                       introduced by contamination. If provided 
                                       the input BAM will be scanned for reads 
                                       that overlap these reads plus the value 
                                       of --flanks instead of processing the 
                                       whole BAM file.''')
    optional_args.add_argument('--flanks', metavar='FLANKS', type=int, 
                               default=500,
                               help='''Amount of bp up and downstream of 
                                       candidate variants (from --vcf argument) 
                                       to scan for reads. Default=500.''')
    optional_args.add_argument('--regions', metavar='REGIONS', nargs='+', 
                               default=[],
                               help='''List of regions (in format chr1:1-1000) 
                                       to scan rather than processing the whole
                                       BAM file.''')
    optional_args.add_argument('-b', '--bwa', metavar='BWA', 
                               help='''Location of bwa executable. Only 
                                       required if not in your PATH.''')
    optional_args.add_argument('-s', '--samtools', metavar='SAMTOOLS', 
                               help='''Location of samtools executable. Only 
                                       required if not in your PATH.''')
    optional_args.add_argument('-m', '--min_fraction_clipped', metavar='FLOAT', 
                               type=float, default=0.2,
                               help='''Minimum proportion of a read that is 
                                       hard or soft-clipped in the BAM file
                                       for a read to be analyzed as a potential 
                                       contaminant. Default = 0.2''')
    optional_args.add_argument('-n', '--min_bases_clipped', metavar='INT', 
                               type=int, default=None, 
                               help='''Minimum number of hard or soft-clipped
                                       bases for a read to be analyzed as a 
                                       potential contaminant. This overrides
                                       the --min_fraction_clipped argument.''')
    optional_args.add_argument('-f', '--fastqs', metavar='PREFIX', 
                               help='''Prefix for fastq files created from 
                                       reads that are tested for contamination.
                                       By default temporary files are created 
                                       and deleted after the program exits. If 
                                       this argument is given, these fastq 
                                       files will be named 
                                       "<prefix>_r1.fastq.gz" and 
                                       "<prefix>_r2.fastq.gz" and will persist
                                       after the program exits.''')
    optional_args.add_argument('-t', '--tmp', metavar='TMPDIR', 
                               help='''Directory to use for temporary files. 
                                       If not specified, the system default 
                                       temporary directory will be used.''')
    optional_args.add_argument('--quiet', action='store_true',  
                               help='''Only output error and warning logging 
                                       information.''')
    optional_args.add_argument('--debug', action='store_true', 
                               help='''Output debugging information.''')
    optional_args.add_argument('--no_caching', action='store_true', 
                               help='''Do not cache any reads to disk and hold 
                                       all as yet unpaired reads in memory. By 
                                       default, reads with a mate mapped to 
                                       another chomosome are cached to disk and 
                                       processed after reading the input file. 
                                       Set this option to True to disable this 
                                       caching at the expense of more RAM 
                                       usage.''')

    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    runner = BamAnalyzer(**vars(args))
    try:
        runner.read_bam()
        runner.align_candidates()
        if runner.output is not None:
            runner.clean_bam()
    finally:
        runner.cleanup()
