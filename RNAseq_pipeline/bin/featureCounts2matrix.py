''' these modules were used.'''
import os
import sys
import argparse

def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict


def main():
    input_files = sys.argv[1:]
    gene_dict = {}
    gene_list = []
    sample_list = []
    for each_file in input_files:
        sample_name = os.path.basename(each_file).split(".")[0]
        sample_list.append(sample_name)
        with open(each_file, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    continue
                if line.startswith("Geneid"):
                    continue
                Geneid,Chr,Start,End, Strand,Length,counts = line.split()
                addtwodimdict(gene_dict,Geneid,sample_name,counts)
                gene_list.append(Geneid)
    gene_list = set(gene_list)
    print("Geneid\t{}".format("\t".join(sample_list)))
    for gene in gene_list:
        gene_count = []
        for sample in sample_list:
            v = gene_dict[gene][sample]
            gene_count.append(v)
        print("{}\t{}".format(gene, "\t".join(gene_count)))

if __name__ == '__main__':
    main()
