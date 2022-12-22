import sys
import pandas as pd

gex_df = pd.DataFrame()
vdj_df = pd.DataFrame()
with open("metrics.summary.txt","r") as f:
    for line in f:
        line = line.rstrip("\n")
        a = line.split("/")
        print(a)
        name = a[-3]
        # print(name)
        df = pd.read_csv(line)
        df['name'] = name[:-4]
        if name.endswith("GEX"):
            gex_df = pd.concat([gex_df, df], ignore_index=True)
        elif name.endswith("VDJ"):
            vdj_df = pd.concat([vdj_df, df], ignore_index=True)

# print(gex_df)
# print(vdj_df)

df2 = pd.merge(gex_df, vdj_df, how="outer",on="name",suffixes=("_GEX","_VDJ"))
df2.to_csv("summary.csv",header=True,index=None)