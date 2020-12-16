import sys
import os
import argparse
from Bio import SeqIO


def usage():
    '''
    see: cd /cygene2/work/P0000-Blackbird/2103-BL001/BL001004/BL001004_WES_RNAseq/BL001004
    example: python extract_sequence_TSS.py -i candidate.info.txt -l 200 -r /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta
    '''


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-i', '--input', action='store', dest='input_file',   help="count file, usually named G260E2L1.pair.count.txt")
    parser.add_argument('-l', '--length',  action='store', dest="extract_len", type=int,  help="frequency file, usually named G260E2L1.pair.acc_freq.txt")
    parser.add_argument('-r', '--reference',action='store', dest="reference_fa", help="reference file")
    return parser.parse_args()


def deal_fasta(fasta):
    fasta_dict = {}
    for record in SeqIO.parse(fasta, "fasta"):
        fasta_dict[str(record.id)] = str(record.seq)
    return fasta_dict


def main():
    parser = _argparse()
    input_file = parser.input_file
    extract_len = parser.extract_len
    reference = parser.reference_fa
    
    ref_dict = deal_fasta(reference)
    with open(input_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            name, chrome, pos = line.split()
            print("{}\t{}\t{}\t{}".format(name, chrome, pos, ref_dict[chrome][int(pos)-extract_len:int(pos)+extract_len]))


if __name__ == '__main__':
    main()
