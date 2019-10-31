import sys


class clonotypes(object):
    """docstring for clonotypes"""
    def __init__(self, barcode,contig_ids,chain,v_gene,j_gene, cdr3,cdr3nt,reads,umis,):
        # super(clonotypes, self).__init__()
        self.barcode, = barcode,
        self.contig_ids
        
def two_dim_dict(thedict, key_a, key_b, value):
    if key_a in thedict:
        if key_b in thedict[key_a]:
            value = thedict[key_a][key_b] + value
            thedict[key_a].update({key_b: value})
        else:
            thedict[key_a].update({key_b: value})
    else:
        thedict.update({key_a: {key_b: value}})
    return thedict

file1 = 'filtered_contig_annotations.csv'

reads_number =0
clonotype_dict = {}
with open(file1, "r") as f:
    for line in f:
        if line.startswith("barcode"):
            continue
        # clonotype_id,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis = line.split(",")
        barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,raw_consensus_id = line.split(",")
        if v_gene.startswith(("TRA","TRB")) and cdr3 != 'None' and productive == 'True':
            if chain == 'TRA':
                reads_number += int(reads)
                two_dim_dict(clonotype_dict,barcode, 'TRA', umis)
            elif chain == 'TRB':
                reads_number += int(reads)
                two_dim_dict(clonotype_dict,barcode, 'TRB', umis)
# print reads_number
# # print clonotype_dict            
for k in clonotype_dict:
    if 'TRA' in clonotype_dict[k]:
        print("{}\tTRA\t{}".format(k,clonotype_dict[k]['TRA']))
    if 'TRB' in clonotype_dict[k]:
        print("{}\tTRB\t{}".format(k,clonotype_dict[k]['TRB']))
