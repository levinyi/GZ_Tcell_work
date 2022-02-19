import sys
import pandas as pd

add_depth_maf = sys.argv[1]
new_maf = sys.argv[2]

list1 = []

# read add_depth_maf and extract all the gene_position info as data 1
# and extract all the tumor data that not equal 0, add these data to data2.
# data = pd.read_csv(add_depth_maf, sep="\t",)
# data2 = data[data['muta_reads_in_Tumor'] != 0]
# print(data2)

with open(add_depth_maf, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("Hugo_Symbol"):
            continue
        a = line.split("\t")
        muta_reads_in_tumor = a[-2]
        if int(muta_reads_in_tumor) > 0:
            list1.append(a[0]+':'+a[4]+":"+a[5])
print(len(list1))
# read all the data 2 and plus data1 and unique.
with open(new_maf) as f2:
    for line in f2:
        line = line.rstrip("\n")
        if line.startswith("Hugo_Symbol"):
            continue
        a = line.split("\t")
        name = a[0]+':'+a[4]+":"+a[5] 
        if name not in list1:
            list1.append(name)

with open(sys.argv[3], "w") as output:
    for line in list1:
        output.write("{}\n".format(line))
# print data 1 and data 2 
# using data 2 venn data 2
print(len(list1))
