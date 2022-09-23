import sys
import pandas as pd
from itertools import product

def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].setdefault(key_b, []).append(val)
    else:
        thedict.update({key_a: {key_b: [val]}})
    return thedict


tcr_table = sys.argv[1]  # LC25TIL1_merged_tcr_table.csv
output_prefix = sys.argv[2] # LC25TIL1

# read data
data = pd.read_csv(tcr_table, header=0, index_col=0)

# deal with data, convert to a dict.
clonotype_dict = {}
for index, row in data.iterrows():
    clonotype = row['original_clone_id']
    chain = row['chain']
    v_gene = row['v_gene']
    j_gene = row['j_gene']
    cdr3 = row['cdr3']
    addtwodimdict(clonotype_dict, clonotype, chain, v_gene+"_"+cdr3+"_"+j_gene)

# deal with dict, calculate clonotype type, and combine AB pair.
cd4posfile = open(output_prefix+"_CD4pos.TCRpairs.xls", "w")
cd8posfile = open(output_prefix+"_CD8pos.TCRpairs.xls", "w")
type_list = []
for clonotype, chains in clonotype_dict.items():
    # how many TRA and TRB are there in the dict?
    A = len(set(chains.get('TRA',[])))
    B = len(set(chains.get('TRB',[])))
    type_list.append(str(A)+'A'+str(B)+'B') # eg.2A2B
    
    # permutation and combination
    loop_val = [list(set(chains.get('TRA',''))), list(set(chains.get('TRB','')))]
    
    # output separately.
    if 'CD4pos' in clonotype: # need check
        for i in product(*loop_val):
            cd4posfile.write("{}\t{}\n".format(clonotype, "_".join(x for x in i)))
    else:
        for i in product(*loop_val):
            cd8posfile.write("{}\t{}\n".format(clonotype, "_".join(x for x in i)))
cd8posfile.close()
cd4posfile.close()

####### output TCR types
import collections
frequency = collections.Counter(type_list)
# print(dict(frequency))
for key, value in dict(frequency).items():
    print("{}\t{}".format(key, value))

