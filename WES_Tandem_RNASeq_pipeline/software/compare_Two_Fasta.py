import multiprocessing
import sys
from Bio import SeqIO
import argparse


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-r', '--ref', action='store', dest='fasta1', help="this is fasta1 file")
    parser.add_argument('-t', '--treat', action='store', dest='fasta2', help="this is fasta2 file")
    parser.add_argument('-f', '--flag', action='store', dest='flag', default="short", help="[long] or [short]")
    parser.add_argument('-o', '--out', action='store', dest='output', help="output file")
    return parser.parse_args()


def deal_fasta(fasta_file, ff):
    fasta_dict = {}
    if ff == 'compare':
        for record in SeqIO.parse(fasta_file, "fasta"):
            fasta_dict[str(record.id)] = record.seq
    elif ff == 'ref':
        for record in SeqIO.parse(fasta_file,"fasta"):
            fasta_dict[str(record.id)] = str(record.seq).upper()
    else:
        sys.exit("wrong flag:[compare] or [ref] ")
    return fasta_dict

def comparing_long(recordid, dict1, dict2, f):
    if len(dict1[recordid]) != len(dict2[recordid]):
        sys.exit("different length: {}".format(recordid))

    seq1 = dict1[recordid]
    seq2 = dict2[recordid]
    for i in range(0, len(seq1)):
        if seq1[i] != seq2[i]:
            f.write("{}\t{}\t{}\t{}\t{}\n".format(recordid, i+1, seq1[i], seq2[i], seq1[i]+str(i+1)+seq2[i]))

def comparing(recordid, seq1, seq2):
    """docstring for comparing"""
    if len(seq1) != len(seq2):
        sys.exit("different length: {}".format(recordid))
    
    bedlist = []
    for i in range(0,len(seq1)):
        if seq1[i] != seq2[i]:
            bedlist.append([recordid,i+1,seq1[i],seq2[i],seq1[i]+str(i+1)+seq2[i]])
    return bedlist

def mycallback(alist):
    """docstring for mycallback"""
    output = _argparse().output
    with open(output,"a+") as f:
        for each in alist:
            f.write("{}\n".format("\t".join(map(str,each))))
    

def main():
    """docstring for main"""
    parser = _argparse()
    fasta1 = parser.fasta1
    fasta2 = parser.fasta2
    flag = parser.flag

    if flag == 'long':
        print("yes long start!")
        fasta1_dict = deal_fasta(fasta1, "ref")
        fasta2_dict = deal_fasta(fasta2, "compare")
        chrome_list = fasta1_dict.keys()
        
        output = _argparse().output
        with open(output, "w") as f:
            for chrome in chrome_list:
                comparing_long(chrome, fasta1_dict, fasta2_dict, f)
        print("yes finished long")

    else:
        pool = multiprocessing.Pool()
        p1 = pool.apply_async(deal_fasta, args=(fasta1,"ref",))
        p2 = pool.apply_async(deal_fasta, args=(fasta2,"compare",))
        fasta1_dict = p1.get()
        fasta2_dict = p2.get()
        pool.close()
        pool.join()
        chrome_list = fasta1_dict.keys()
        p = multiprocessing.Pool()
        for chrome in chrome_list:
            p.apply_async(comparing, args=(chrome, fasta1_dict[chrome], fasta2_dict[chrome],), callback=mycallback)
    
        p.close()
        p.join()


if __name__ == '__main__':
    main()
