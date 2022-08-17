import sys
import os
import pandas as pd

maf_file = sys.argv[1]

# read in maf file

data = pd.read_csv(maf_file, sep="\t",)
print(data.head())

for index, row in data.iterrows():
    # if chromosome is null pass
    
    if row['Chromosome'] == 'chrM':
        continue
    else