import sys
import pandas as pd

tcr_stats = sys.argv[1]   # tcr_stats_barcodes.csv
merged_tcr = sys.argv[2]   # merged_tcr_table.csv

noGex_clone = []
Gex_clone = []

for index, row in pd.read_csv(tcr_stats).iterrows():
    if pd.isnull(row['YostEx_AUCscore']) and pd.isnull(row['TCRx_sig_AUCscore']):
        noGex_clone.append(row['renamed_clone'])
    else:
        Gex_clone.append(row['renamed_clone'])


total_clonotype_dict = {}
for index, row in pd.read_csv(merged_tcr).iterrows():
    total_clonotype_dict[row['raw_clonotype_id']] = row["cdr3s_aa"]

noGex_cdr3_dict = {}
for clone in noGex_clone:
    if clone in total_clonotype_dict:
        noGex_cdr3_dict[clone] = total_clonotype_dict[clone]

Gex_cdr3_dict = {}
for clone in Gex_clone:
    if clone in total_clonotype_dict:
        Gex_cdr3_dict[clone] = total_clonotype_dict[clone]

noGex_cdr3_dict = {k:v for k,v in noGex_cdr3_dict.items() if k not in Gex_cdr3_dict}

# print(len(total_clonotype_dict))
# print(len(Gex_cdr3_dict))
# print(len(noGex_cdr3_dict))

# compare two dict
# print(Gex_cdr3_dict)
# 创建一个空列表来存储匹配的键和值
matches = []

# 检查 noGex_cdr3_dict 中的每个值是否在 Gex_cdr3_dict 中,
for key1, value1 in noGex_cdr3_dict.items():
    for cdr3 in [x.strip() for x in value1.split(';')]:
        if cdr3 in [x.strip() for v in Gex_cdr3_dict.values() for x in v.split(';')]:
            # 如果找到匹配，则将匹配的键和值添加到 matches 列表中
            for key2, value2 in Gex_cdr3_dict.items():
                if cdr3 in [x.strip() for x in value2.split(';')]:
                    matches.append((key1, cdr3, key2, value2))
        else:
            matches.append((key1,cdr3))

# 输出表格
# print("noGex_cdr3_key\tnoGex_cdr3_value\tGex_cdr3_key\tGex_cdr3_value")
for match in matches:
    print("\t".join(match))
