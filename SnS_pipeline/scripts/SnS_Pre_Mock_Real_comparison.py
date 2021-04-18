import sys
import pandas as pd

def deal_data(afile):
    df = pd.read_csv(afile)
    df_col2list = df.columns.tolist()
    df_pair = df_col2list[0]
    df_real_f = df_col2list[7]  # real_freq
    df_mp = df_col2list[10]  # mock/pre
    df_rp = df_col2list[11] # real/pre
    df_rm = df_col2list[12] # real/mock
    # print(df_pair,df_real_f,df_mp,df_rp,df_rm)
    df = df[(df[df_rp] > 3) & (df[df_real_f]>0.0001)]
    df = df[[df_pair,df_real_f,df_mp,df_rp,df_rm]]
    return df

df1 = deal_data(sys.argv[1])
df2 = deal_data(sys.argv[2])
df = pd.merge(df1,df2,on="TCR_id(Real_samples_union)", how="inner")
# print(df)
df.to_csv(sys.argv[3], sep=",",index=False)
print("done")
