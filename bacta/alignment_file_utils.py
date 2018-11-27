import sys
import os
import pysam
import tempfile

def get_bamfile(bam):
    bmode = 'rb'
    if bam.endswith(('.sam', '.SAM')):
        bmode = 'r'
    elif bam.endswith(('.cram', '.CRAM')):
        bmode = 'rc'
    return pysam.AlignmentFile(bam, bmode)

def get_align_output(output, template=None, header=None):
    if output is not None:
        wbmode = 'wb'
        if output.endswith(('.sam', '.SAM')):
            wbmode = 'w'
        elif output.endswith(('.cram', '.CRAM')):
            wbmode = 'wc'
        return pysam.AlignmentFile(output, wbmode, template=template,
                                   header=header)
    return  pysam.AlignmentFile('-', "w", template=template, header=header)

def seek_back_til_reads(bamfile, start_tid=None, ref=None):
    if start_tid is None:
        start_tid = bamfile.nreferences - 1
    #sys.stderr.write("Seeking backwards through contigs until we find reads.\n")
    got_one = 0
    for i in range(start_tid, -1, -1):
        #sys.stderr.write("Trying contig {}\n"
        #                 .format(bamfile.get_reference_name(i)))
        for read in bamfile.fetch(tid=i, reference=ref):
         #   if not got_one:
         #       sys.stderr.write("At contig {}\n"
         #                        .format(bamfile.get_reference_name(i)))
            got_one += 1
        if got_one:
            return

def get_tmp_bam(tmpdir=None, template=None, header=None):
    '''
        Returns a writeable AlignmentFile with a path created by
        tempfile.mkstemp, using bamfile as a template.
    '''
    (f, tmp_bam) = tempfile.mkstemp(suffix='.bam', dir=tmpdir)
    os.close(f)
    return (tmp_bam,
            pysam.AlignmentFile(tmp_bam, "w", template=template, header=header)
           )
