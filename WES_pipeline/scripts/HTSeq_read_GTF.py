#! coding:UTF-8
import pandas as pd


gtf_file = 'gencode.v33lift37.annotation.gtf'

gene_dict = {}
with open(gtf_file, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        a = line.split()
        gene_name = a[9].split(".")[0].lstrip("\"")
        if a[2] == 'exon':
            gene_dict.setdefault(gene_name, []).append([int(a[3]), int(a[4])])

gene_length_dict = {}
for gene_name, intervals in gene_dict.items():
    if len(intervals) < 2:
        res = intervals
    res = []
    intervals.sort(key=lambda x: x[0])
    res.append(intervals[0])
    for i in range(1, len(intervals)):
        a = res[-1]
        b = intervals[i]
        if a[1] >= b[0]:
            # a[1] = max(a[1], b[1])  
            # 找到最大的 end 用下面这种可以减少一半的运行
            if a[1] < b[1]:
                a[1] = b[1]
        else:
            res.append(b)
    
    total_exon_length = 0
    for each in res:
        start, end = each
        total_exon_length += end - start

    gene_length_dict[gene_name] = total_exon_length

# from itertools import islice
gene_length = pd.DataFrame.from_dict(gene_length_dict, orient='index', columns=["gene_length"])
#print(gene_length.head(10))

data = pd.read_table("count.txt", sep="\t", header=0, index_col=0)
#print(data.head(10))

big_data = pd.concat([gene_length, data], axis=1)
#print(big_data.head(10))

big_data2 = big_data.loc[:,"GSM4054837":].div(big_data['gene_length'], axis=0)
#print(big_data2.head(10))
#print(big_data2['GSM4054837'].sum())
big_data3 = big_data2*1000000

big_data4 = big_data3.loc[:,:].div(big_data2.sum(axis=0), axis=1)
#print(big_data3.head(10))
#print(big_data4.head(10))
# big_data4.dropna(axis=0, how='all', inplace=True)
big_data4.to_csv("count.TPM.txt",sep="\t")
