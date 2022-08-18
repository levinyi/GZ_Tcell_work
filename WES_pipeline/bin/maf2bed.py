import sys
import pandas as pd

maf_file = sys.argv[1]

# read in maf file
data = pd.read_csv(maf_file, sep="\t",)
print(data.head())

for index, row in data.iterrows():
    # if chromosome is null pass
    if pd.isnull(row["Chromosome"]):
        continue
    if row['Chromosome'] == 'chrM':
        start = row['Start_Position']-1
        end = row['End_Position'] -1

print(row['Chromosome'], start, end, row['Tumor_Sample_Barcode'], sep="\t")