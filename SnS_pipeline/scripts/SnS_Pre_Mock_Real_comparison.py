import sys,os
import pandas as pd
import argparse


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-a', '--fileA', action='store', dest='filea', help="this is config file")
    parser.add_argument('-b', '--fileB', action='store', dest='fileb', help='this is project')
    parser.add_argument('-a1', '--realFreqColumn_a', action='store', type=int, dest='real_freq_col_a', help='column number for real_freq_column, 0-based. 0 for absent')
    parser.add_argument('-a2', '--mockPreColumn_a', action='store', type=int, dest='df_mp_col_a', help='column number for mock/pre column, 0-based. 0 for absent.')
    parser.add_argument('-a3', '--realPreColumn_a', action='store', type=int, dest='df_rp_col_a', help='column number for real/pre column, 0-based. 0 for absent.')
    parser.add_argument('-a4', '--realMockColumn_a', action='store', type=int, dest='df_rm_col_a', help='column number for real/mock column, 0-based. 0 for absent.')
    parser.add_argument('-b1', '--realFreqColumn_b', action='store', type=int, dest='real_freq_col_b', help='column number for real_freq_column, 0-based. 0 for absent')
    parser.add_argument('-b2', '--mockPreColumn_b', action='store', type=int, dest='df_mp_col_b', help='column number for mock/pre column, 0-based. 0 for absent.')
    parser.add_argument('-b3', '--realPreColumn_b', action='store', type=int, dest='df_rp_col_b', help='column number for real/pre column, 0-based. 0 for absent.')
    parser.add_argument('-b4', '--realMockColumn_b', action='store', type=int, dest='df_rm_col_b', help='column number for real/mock column, 0-based. 0 for absent.')
    parser.add_argument('-c', '--realPreFoldChange', action='store', type=int, dest='realPreFoldChange',default=3, help='Real/Pre fold change threshold: default: 3')
    parser.add_argument('-f', '--filterRealFreq', action='store', type=float, dest='filter_real_f',default=0.0001, help='filter real_freq greater than : default: 0.0001')
    parser.add_argument('-o', '--outputPrefix', action='store', type=str, dest='outputPrefix', help="out put file name prefix.")
    return parser.parse_args()


def deal_data(afile, real_freq_col, realPreFoldChange, filter_real_f, df_mp_col, df_rp_col, df_rm_col):
    df = pd.read_csv(afile)
    df_col2list = df.columns.tolist()
    print("total_len header: {} ".format(len(df_col2list)))
    
    df_pair = df_col2list[0]
    df_real_f = df_col2list[real_freq_col]  # real_freq
    
    extract_item = [df_pair, df_real_f]
    if df_mp_col != 0:
        df_mp = df_col2list[df_mp_col]  # mock/pre
        extract_item.append(df_mp)
    else:
        df_mp = "NULL"
    df_rp = df_col2list[df_rp_col]  # real/pre
    extract_item.append(df_rp)
    if df_rm_col != 0:
        df_rm = df_col2list[df_rm_col]  # real/mock
        extract_item.append(df_rm)
    else:
        df_rm = "NULL"

    # print(df_pair,df_real_f,df_mp,df_rp,df_rm)
    print("checking column number...")
    print("{}: real_freq:{}:{}\tmock/pre:{}:{}\treal/pre:{}:{}\treal/mock:{}:{}".format(
        os.path.basename(afile), 
        real_freq_col, df_real_f, 
        df_mp_col, df_mp,
        df_rp_col, df_rp,
        df_rm_col, df_rm,))
    df = df[(df[df_rp] > realPreFoldChange) & (df[df_real_f]>filter_real_f)]
    #df = df[[df_pair, df_real_f, df_mp, df_rp, df_rm]]
    df = df[extract_item]
    return df


def main():
    p = _argparse()
    if not p.real_freq_col_b:
        p.real_freq_col_b = p.real_freq_col_a
    if not p.df_mp_col_b:
        p.df_mp_col_b = p.df_mp_col_a
    if not p.df_rp_col_b:
        p.df_rp_col_b = p.df_rp_col_a
    if not p.df_rm_col_b:
        p.df_rm_col_b = p.df_rm_col_a
    df1 = deal_data(p.filea, p.real_freq_col_a, p.realPreFoldChange, p.filter_real_f, p.df_mp_col_a, p.df_rp_col_a, p.df_rm_col_a)
    df2 = deal_data(p.fileb, p.real_freq_col_b, p.realPreFoldChange, p.filter_real_f, p.df_mp_col_b, p.df_rp_col_b, p.df_rm_col_b)
    df_merge_inner = pd.merge(df1,df2, on="TCR_id(Real_samples_union)", how="inner")
    df_merge_inner.to_csv(p.outputPrefix + '.csv', sep=",",index=False)
    print("done")

if __name__ == '__main__':
    main()
