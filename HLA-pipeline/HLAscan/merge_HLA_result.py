import sys
import os

input_files = sys.argv[1:]

def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict


def deal_hla(afile):
    type1 = 'NULL'
    type2 = 'NULL'
    name = afile.split(".")[2]
    with open(afile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("HLA gene"):
                a,name = line.split(":")
                name = name.lstrip()
                continue
            if line.startswith("[Type 1]"):
                c = line.split("\t")
                type1 = c[1]
                continue
            if line.startswith("[Type 2]"):
                c = line.split("\t")
                type2 = c[1]
                continue
    return name, type1, type2


results = {}
hla_types = ["HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","MICA","MICB","HLA-DMA","HLA-DMB",
        "HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1",
        "HLA-DRB5","TAP1","TAP2"]
for each_file in input_files:
    name, type1, type2 = deal_hla(each_file)
    addtwodimdict(results,name,'type1',type1)
    addtwodimdict(results,name,'type2',type2)

for each in hla_types:
    if each in results:
        print("{}-1\t{}".format(each, results[each].get("type1",'NULL')))
        print("{}-2\t{}".format(each, results[each].get("type2",'NULL')))
    else:
        print("{}-1\t{}".format(each, "NULL"))
        print("{}-1\t{}".format(each, 'NULL'))
