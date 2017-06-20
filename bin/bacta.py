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
                               help='''Output BAM filename. Defaults to the 
                                       basename of the input + "_bacta.bam".
                                    ''')
    optional_args.add_argument('-c', '--contaminants', metavar='PREFIX', 
                               help='''Prefix for contaminant output files. 
                                       Defaults to the basename of the input + 
                                       "_contaminants".
                                    ''')
    optional_args.add_argument('-b', '--bwa', metavar='BWA', 
                               help='''Location of bwa executable. Only 
                                       required if not in your PATH.''')
    optional_args.add_argument('-m', '--min_fraction_clipped', metavar='FLOAT', 
                               default=0.2, 
                               help='''Minimum proportion of a read that is 
                                       hard or soft-clipped in the BAM file
                                       for a read to be analyzed as a potential 
                                       contaminant. Default = 0.2''')
    optional_args.add_argument('-n', '--min_bases_clipped', metavar='INT', 
                               default=None, 
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

    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    runner = BamAnalyzer(**vars(args))
    runner.read_bam()
