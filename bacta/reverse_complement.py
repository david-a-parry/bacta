import sys
if sys.version_info[0] < 3:
    from string import maketrans
else:
    maketrans = str.maketrans
comp = maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(dna):
    return complement(dna[::-1])

def complement(dna):
    #global trans
    return dna.translate(comp)
