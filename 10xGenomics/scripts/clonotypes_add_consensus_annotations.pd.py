# %%
import os,sys
import pandas as pd
from sklearn.tree import plot_tree


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


def read_consensus(file2):
    len_dict = {}
    detail_dict = {}
    with open(file2, "r") as f2:
        for line in f2:
            line = line.rstrip("\n")
            if line.startswith("clonotype_id"):
                continue
            a = line.split(",")
            clonotype_id = a[0]
            chain = a[3]

            two_dim_dict(len_dict,clonotype_id,chain,1)
            if clonotype_id in detail_dict:
                if chain in detail_dict[clonotype_id]:
                    new_chain = chain + "_1"
                    detail_dict.setdefault(clonotype_id,{})[new_chain] = a[1:]
                else:
                    detail_dict.setdefault(clonotype_id,{})[chain] = a[1:]
            else:
                detail_dict.setdefault(clonotype_id,{})[chain] = a[1:]
    return len_dict, detail_dict

# %%
file1 = "clonotypes.csv"
file2 = "consensus_annotations.csv"

clonotype = pd.read_csv(file1)
print(clonotype.head())

annotate = pd.read_csv(file2)
print(annotate.head())

clonotype_dict = clonotype.set_index("clonotype_id").to_dict('index')
print(dict(list(clonotype_dict.items())[:3]))

#anno_dict = annotate.set_index("clonotype_id").to_dict('index')
anno_dict = annotate.set_index(["clonotype_id",'chain']).to_dict('index')
# anno_dict = {k: g.to_dict('records') for k, g in annotate.groupby(level=0)}
print(dict(list(anno_dict.items())[:3]))
# anno_dict = annotate.groupby(by='clonotype_id', sort=False).apply(lambda x:x.to_dict(orient='index'))
# print(dict(list(anno_dict.items())[:3]))


def recur_dictify(frame):
    if len(frame.columns) == 1:
        if frame.values.size == 1: return frame.values[0][0]
        return frame.values.squeeze()
    grouped = frame.groupby(frame.columns[0])
    d = {k: recur_dictify(g.iloc[:,1:]) for k,g in grouped}
    return d
# anno_dict = recur_dictify(annotate)
# print(dict(list(anno_dict.items())[:3]))

'''
flag = 0 
# default 1 TRA, 1 TRB if sum(TRA + TRB) >2
if len() >2:
    index += 1
    for each in TRA:
        for each in TRB:
            pd.dataframe()
'''
