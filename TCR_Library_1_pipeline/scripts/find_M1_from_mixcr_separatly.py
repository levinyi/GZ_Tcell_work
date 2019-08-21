import sys
import os
import gzip
from Bio import SeqIO


def usage():
    """
    python find_M1_from_mixcr_separatly.py G37E2L1_R1.mixcr.out.clonotypes.TRB.txt >G37E2L1_R1.TRB.reads
    
    """

    
mixcr = sys.argv[1]

sample_name = os.path.basename(mixcr).split(".")[0]
with open(mixcr, "r") as f:
    for line in f:
        line = line.rstrip("\t")
        if line.startswith("cloneId"):
            continue
        if len(line.split("\t")) > 15:
            cloneId = line.split("\t")[0]  # cloneId
            bestVGene = line.split("\t")[5].split("*")[0]  # allVHitsWithScore
            bestJGene = line.split("\t")[7].split("*")[0]  # allJHitsWithScore
            aaSeqCDR3 = line.split("\t")[32]  #
        else:
            cloneId, cloneCount, bestVGene, bestJGene, nSeqCDR3, aaSeqCDR3, qualCDR3 = line.split("\t")
        fastq_file = "{0}/{0}.{1}.fastq.gz".format(sample_name, cloneId)
        fq = gzip.open(fastq_file, "r")
        for record in SeqIO.parse(fq, "fastq"):
            print("{0}\t{1}".format(record.id, cloneId))
