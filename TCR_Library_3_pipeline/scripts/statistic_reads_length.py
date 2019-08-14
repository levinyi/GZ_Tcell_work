import sys
from Bio import SeqIO

adict = {}
index = 1
for record in SeqIO.parse(sys.argv[1],"fastq"):
    print("%s\t%s"%(index, len(record.seq)))
    index += 1

