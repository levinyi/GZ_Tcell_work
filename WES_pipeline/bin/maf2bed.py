import sys
import pandas as pd

maf_file = sys.argv[1]
data = pd.read_csv(maf_file, sep="\t",dtype=str)

for index, row in data.iterrows():
    if pd.isnull(row["Chromosome"]):
        continue
    
    start = int(row['Start_Position'])-1
    if row['Chromosome'] == 'MT':
        chrome = "chrM"
    else:
        chrome = row['Chromosome']
    
    print("{}\t{}\t{}".format(chrome, start, row['End_Position']))