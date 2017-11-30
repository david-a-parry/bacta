import sys
import os
import pandas as pd
from io import StringIO

def write_table(counts, fh):
    #create table as string
    #first the header...
    k, *_ = counts.keys()
    table = "Threshold\t" + str.join("\t", (x for x in sorted(counts[k])))+"\n"
    #...then results
    for t in sorted(counts.keys()):
        table += "{:g}\t".format(t)
        table += str.join("\t", (str.format("{:.3g}", counts[t][x]) 
                                 for x in sorted(counts[t]))) + "\n"
    #use table string as input to pandas to create dataframe
    tab = StringIO(table)
    df = pd.read_table(tab)
    #calculate stats
    total_ref = df['REF'].sum()
    total_contam = df['CONTAM'].sum()
    tc = 0
    tr = 0
    cum_contam = []
    cum_ref = []
    for ref in reversed(df.REF): 
        tr += ref
        cum_ref.append(tr)
    for contam in reversed(df.CONTAM): 
        tc += contam
        cum_contam.append(tc)
    cum_contam.reverse()
    cum_ref.reverse()
    df['TP'] = cum_contam
    df['FP'] = cum_ref
    df['Sensitivity'] = df.TP/total_contam
    df['1-Specificity'] = df.FP/total_ref
    df['Precision'] = df.TP/(df.TP + df.FP)
    df.to_csv(fh, sep='\t', index=False)
    fh.close()
    return df

def report_precision_recall(df, logger):
    for row in df.itertuples():
        logger.info("At threshold={} precision={:g}, recall={:g}"
                    .format(row.Threshold, row.Precision, row.Sensitivity))

def contigs_from_fasta(fasta):
    contigs = set()
    fai = fasta + '.fai'
    if not os.path.exists(fai):
        sys.exit("ERROR: could not find fasta index ('{}') for fasta "
                 .format(fai) + "reference. Please index ('samtools faidx"+
                 " {}') before running this program." .format(fasta))
    with open(fai, 'r') as fh:
        for line in fh:
            c = line.split()[0]
            if c in contigs:
                sys.exit("ERROR: Fasta index '{}' appears to be malformed."
                         .format(fai) + " Contig {} was found more than once" 
                         .format(c))
            contigs.add(c)
    return contigs



