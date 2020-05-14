import multiprocessing
import sys
from Bio import SeqIO
import argparse


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-r', '--ref', action='store', dest='comp_file', help="this is fasta1 file")
    parser.add_argument('-t', '--stf', action='store', dest='STF_file', help="this is fasta2 file")
    parser.add_argument('-f', '--flag',action='store', dest='flag', help='[genome] or [CDS]')
    parser.add_argument('-o', '--out', action='store', dest='output', help="output file")
    return parser.parse_args()

def deal_comp_file(afile):
    adict = {}
    # NM_015102.5     2482    G       A       G2482A
    with open(afile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            a,b,c,d,e = line.split("\t")
            a = a.split(".")[0]
            adict.setdefault(a,[]).append(e)
    return adict


def main():
    parser = _argparse()
    com_file = parser.comp_file
    STF_file = parser.STF_file
    output = parser.output
    flag = parser.flag

    adict = deal_comp_file(com_file)

    with open(STF_file, "r") as f, open(output, "w") as f2:
        if flag == "genome":
            for line in f:
                line = line.rstrip("\n")
                chrome = line.split("\t")[0]
                startpoint = line.split("\t")[1]
                ref = line.split("\t")[3]
                alt = line.split("\t")[4]
                nt_change = ref+startpoint+alt

                if chrome in adict:
                    if nt_change in adict[chrome]:
                        f2.write("{}\tYES\t{}\n".format(line,nt_change))
                    else:
                        f2.write("{}\tNO\t{}\n".format(line, adict[chrome]))
                else:
                    f2.write("line not in compare file\t{}\n".format(line))

        elif flag == "CDS":
        # chr1  22172744    22172744    C   T   missense SNV    HSPG2   NM_005529   exon64  c.G8321A    p.R2774H
            for line in f:
                line = line.rstrip("\n")
                NM = line.split("\t")[7]
                cm = line.split("\t")[9].lstrip("c.")
            
                if NM in adict:
                    if cm in adict[NM]:
                        f2.write("{}\tYES\t{}\n".format(line, cm))
                    else:
                        f2.write("{}\tNO\t{}\n".format(line,adict[NM]))
                else:
                    f2.write("line not in compare file\t{}\n".format(line))
        else:
            sys.exit("wrong flag for input.")


if __name__ == '__main__':
    main()