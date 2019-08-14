import sys
import os
import gzip
from Bio import SeqIO

def usage():
    '''
    python pairs2umi.py G13E3L8.pairs >G13E3L8.pairs.umi.txt
    '''

def get_handle(file):
    if os.path.basename(file).endswith("gz"):
        handle = gzip.open(file, "rU")
    elif os.path.basename(file).endswith(("fq", "fastq")):
        handle = open(file, "rU")
    return handle

def deal_umi_file(afastq):
    handle = get_handle(afastq)
    umi_dict = {}
    for record in SeqIO.parse(handle,"fastq"):
        umi_dict[record.id] = str(record.seq)[:8]
    return umi_dict

def main():
    pairs_file = sys.argv[1]

    name = os.path.basename(pairs_file).split(".")[0]
    umi_file = '../data/'+ name+'_R2.fq.gz'
    if not os.path.exists(umi_file) :
        print("umi_file is not exists : %s " % umi_file)
        sys.exit()

    umi_dict = deal_umi_file(umi_file)

    # deal pairs file.
    with open(pairs_file,"r") as f:
        for line in f:
            line =  line.rstrip("\n")
            readid,cloneA,cloneB = line.split("\t")
            cloneA = 'a'+cloneA
            cloneB = 'b'+cloneB
            print('{pairID}\t{readID}\t{UMI}'.format(pairID=cloneA +'_'+ cloneB, readID=readid,UMI=umi_dict[readid]))

if __name__ == '__main__':
    main()