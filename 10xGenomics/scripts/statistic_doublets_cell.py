import sys
import json


def usage():
    """docstring for usage
    python statistic_doublets_cell.py
    input file is : filtered_contig_annotations.csv
    """


def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].setdefault(key_b, []).append(val)
    else:
        thedict.update({key_a: {key_b: [val]}})
    return thedict


inputfile = "filtered_contig_annotations.csv"

cell_object = {}
droplet_object = {}
with open(inputfile, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("barcode"):
            continue
        barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,raw_consensus_id = line.split(",")
        if chain == "TRA":
            addtwodimdict(cell_object, barcode, "TRA", v_gene)
            addtwodimdict(droplet_object, raw_clonotype_id, barcode, v_gene)
        elif chain == "TRB":
            addtwodimdict(cell_object, barcode, "TRB", v_gene)
            addtwodimdict(droplet_object, raw_clonotype_id, barcode, v_gene)

#print(json.dumps(cell_object, indent=1))
print(json.dumps(droplet_object, indent=1))
