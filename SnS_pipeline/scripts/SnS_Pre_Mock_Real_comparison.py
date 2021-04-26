import sys
import pandas as pd
import argparse

def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-a', '--fileA', action='store', dest='filea', help="this is config file")
    parser.add_argument('-b', '--fileB', action='store', dest='fileb', help='this is project')
    parser.add_argument('-r', '--realFreqColumn', action='store', type = int, dest='real_freq_col', help='the column number for real_freq_column')
    parser.add_argument('-c', '--realPreFoldChange,', action='store', type=int, dest='realPreFoldChange',default=3, help='Real/Pre fold change threshold: default: 3')
    parser.add_argument('-f', '--filterRealFreq', action='store', type=float, dest='filter_real_f',default=0.0001, help='filter real_freq greater than : default: 0.0001')
    parser.add_argument('-o', '--outputPrefix', action='store', type=str, dest='outputPrefix', help="out put file name prefix.")
    return parser.parse_args()

def deal_data(afile, real_freq_col,realPreFoldChange,filter_real_f):
    df = pd.read_csv(afile)
    df_col2list = df.columns.tolist()
    
    df_pair = df_col2list[0]
    df_real_f = df_col2list[real_freq_col]  # real_freq
    df_mp = df_col2list[-3]  # mock/pre
    df_rp = df_col2list[-2]  # real/pre
    df_rm = df_col2list[-1]  # real/mock
    # print(df_pair,df_real_f,df_mp,df_rp,df_rm)
    df = df[(df[df_rp] > realPreFoldChange) & (df[df_real_f]>filter_real_f)]
    df = df[[df_pair, df_real_f, df_mp, df_rp, df_rm]]
    return df, df_raw


def main():
    parser = _argparse()
    df1, df1_raw = deal_data(parser.filea,parser.real_freq_col,parser.realPreFoldChange,parser.filter_real_f)
    df2, df2_raw = deal_data(parser.fileb,parser.real_freq_col,parser.realPreFoldChange,parser.filter_real_f)
    df_merge_inner = pd.merge(df1,df2, on="TCR_id(Real_samples_union)", how="inner")
    df_merge_inner.to_csv(parser.outputPrefix + '.csv', sep=",",index=False)
    print("done")

if __name__ == '__main__':
    main()