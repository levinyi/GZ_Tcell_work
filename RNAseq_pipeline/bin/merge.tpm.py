import sys
import os

def usage():
    '''
    this script only used for XXX.exon.gene_name.counts.txt.TPM.txt file.

    usage:
    python merge.tpm.py /cygene2/work/Tumor_Expression_Data/*.TPM.txt >test.txt

    update:
    20210204. added usage description.
    2020xxxx. created.
    '''
    pass

def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict

gene_dict = {}
gene_list = []
sample_list = []
for each in sys.argv[1:]:
    sample_name = os.path.basename(each).split(".")[0]
    sample_list.append(sample_name)
    with open(each, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            gene,tpm = line.split()
            addtwodimdict(gene_dict,gene,sample_name,tpm)
            gene_list.append(gene)

gene_list = set(gene_list)
print("Geneid\t{}".format("\t".join(sample_list)))
for gene in gene_list:
    gene_tpm = []
    for sample in sample_list:
        v = gene_dict[gene][sample]
        gene_tpm.append(v)
    print("{}\t{}".format(gene, "\t".join(gene_tpm)))
