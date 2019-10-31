import sys

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

file1 = 'consensus_annotations.csv'


clonotype_dict = {}
with open(file1, "r") as f:
    for line in f:
        if line.startswith("clonotype_id"):
            continue
        clonotype_id,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis = line.split(",")
        if chain == 'TRA':
            two_dim_dict(clonotype_dict,clonotype_id, 'TRA', 1)
        else:
            two_dim_dict(clonotype_dict,clonotype_id, 'TRB', 1)
# print clonotype_dict

TRA_1_TRB_1 = 0
TRA_0_TRB_1 = 0
TRA_1_TRB_0 = 0
TRA_2_TRB_1 = 0
TRA_1_TRB_2 = 0
TRA_2_TRB_2 = 0
others =0
for k in clonotype_dict:
    if 'TRA' in clonotype_dict[k] and 'TRB' in clonotype_dict[k]:
        if clonotype_dict[k]['TRA'] == 1 and clonotype_dict[k]['TRB'] ==1:
            TRA_1_TRB_1 += 1
        elif clonotype_dict[k]['TRA'] == 2 and clonotype_dict[k]['TRB'] ==1:
            TRA_2_TRB_1 += 1
        elif clonotype_dict[k]['TRA'] == 1 and clonotype_dict[k]['TRB'] ==2:
            TRA_1_TRB_2 += 1
        elif clonotype_dict[k]['TRA'] == 2 and clonotype_dict[k]['TRB'] ==2:
            TRA_2_TRB_2 += 1
        else:
            others += 1
    elif 'TRA' not in clonotype_dict[k] and 'TRB' in clonotype_dict[k]:
        TRA_0_TRB_1 += 1 
    elif 'TRA' in clonotype_dict[k] and 'TRB' not in clonotype_dict[k]:
        TRA_1_TRB_0 += 1
    else:
        others += 1

total  = TRA_1_TRB_1 + TRA_0_TRB_1 + TRA_1_TRB_0 + TRA_2_TRB_1 + TRA_1_TRB_2 + TRA_2_TRB_2 + others 
print("Total_Types:\t{}\t-\nTRA_1_TRB_1:\t{}\t{:.2%}\nTRA_0_TRB_1:\t{}\t{:.2%}\nTRA_1_TRB_0:\t{}\t{:.2%}\nTRA_2_TRB_1:\t{}\t{:.2%}\nTRA_1_TRB_2:\t{}\t{:.2%}\nTRA_2_TRB_2:\t{}\t{:.2%}\nOthers_Types:\t{}\t{:.2%}".format(\
    total, TRA_1_TRB_1, float(TRA_1_TRB_1)/total, \
    TRA_0_TRB_1, float(TRA_0_TRB_1)/total, \
    TRA_1_TRB_0, float(TRA_1_TRB_0)/total, \
    TRA_2_TRB_1, float(TRA_2_TRB_1)/total, \
    TRA_1_TRB_2, float(TRA_1_TRB_2)/total, \
    TRA_2_TRB_2, float(TRA_2_TRB_2)/total,\
    others, float(others)/total))