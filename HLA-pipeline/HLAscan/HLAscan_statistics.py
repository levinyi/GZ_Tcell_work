import os
import sys


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

def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict

# read data and split path to folder,
# sample_name, sample_type,hla_type
import pandas as pd
from collections import namedtuple
Item = namedtuple('sample_name',['sample_name','sample_type','hla_type','type1','type2'])
items = []

results = {}
input_files = sys.argv[1:]
for each_file in input_files:
    file_basename = os.path.basename(each_file)
    sample_name,sample_type, hla_type, out, txt  = file_basename.split(".")
    # print(sample_name, sample_type, hla_type)

    # get hla type
    name, type1, type2 = deal_hla(each_file)
    items.append(Item(sample_name, sample_type, hla_type, type1,type2))
    '''
    # merge sample_name and sample_type
    sample_name_and_type = sample_name + '.' + sample_type
    # store type1 and type2 to a dict 
    addtwodimdict(results,sample_name_and_type, 'type1', type1)
    addtwodimdict(results,sample_name_and_type, 'type2', type2)
    '''

df = pd.DataFrame.from_records(items, columns=['sample_name','sample_type','hla_type','type1','type2'])
print(df)
df = df.
results = df.to_dict('inndex')

# df2 = df.pivot(index=['hla_type'], columns=['sample_name','sample_type'],values=['type1','type2'])
# df2.to_csv("test.txt", sep="\t", index=True)
# df2= df.pivot(index=['hla_type'],columns=['sample_name','sample_type'])
# print(df2)
# print(df2.melt())
# print hla types
hlascan_types = ["HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","MICA","MICB","HLA-DMA","HLA-DMB",
        "HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1",
        "HLA-DRB5","TAP1","TAP2"]
import json
print(json.dumps(results, indent=4))
# print(json.dumps(dict(list(results.items())[:5]),indent=2))
# output results:
'''
for each_hla in hlascan_types:
    # print("{},".format(each_hla))
    for k, v in results.items():
        print("{}\t{}\t{}".format(each_hla, k, v))
'''
