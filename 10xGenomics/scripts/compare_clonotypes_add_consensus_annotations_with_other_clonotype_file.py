import sys,os
import pandas as pd
from venn import venn
from matplotlib import pyplot as plt


def deal_with_clonotypes(afile, sep, sample_name, assigned_field, trim):
    """read file and trim TRAV,"""
    if sep == "comma":
        df = pd.read_csv(afile, sep=",", header=0)
    elif sep == "tab":
        df = pd.read_csv(afile, sep="\t", header=0)
    
    if trim:
        for each_field in assigned_field:
            df[each_field] = df[each_field].str.replace(".0.$", "", regex=True)
    
    df['clonotype_pair_id'] = df[assigned_field].apply("_".join, axis=1)

    df = df[['clonotype_pair_id','clonotype_id']]
    df.rename(columns = {"clonotype_id":'clonotype_id_'+sample_name})
    
    return df


def main():
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    file1_assigned_field = ['v_gene','cdr3','v_gene.1','cdr3.1']
    file2_assigned_field = ['TRAV','CDR3A','TRBV','CDR3B']
    sample1_name = sys.argv[3]
    sample2_name = sys.argv[4]
    sep1 = "comma"
    sep2 = "tab"

    df1 = deal_with_clonotypes(afile=file1, sep=sep1, sample_name=sample1_name, assigned_field=file1_assigned_field, trim=False)
    df2 = deal_with_clonotypes(afile=file2, sep=sep2, sample_name=sample2_name, assigned_field=file2_assigned_field, trim=True )
    
    #big_DataFrame = pd.DataFrame()
    big_DataFrame = pd.merge(df1, df2, on='clonotype_pair_id', how='outer')
    
    # to big dataframe
    sample_name_list = [sample1_name,sample2_name] # for sample name
    big_DataFrame.to_csv("_".join(sample_name_list)+ ".xls", sep="\t", index=False)
    
    # to venn 
    venn_dict = {}
    venn_dict.setdefault(sample1_name, set(df1['clonotype_pair_id'].tolist()))
    venn_dict.setdefault(sample2_name, set(df2['clonotype_pair_id'].tolist()))

    # print(big_DataFrame.head())
    #print(venn_dict)
    print("write xls file: {}.xls".format("_".join(sample_name_list)))

    venn(venn_dict)
    plt.savefig("_".join(sample_name_list)+".venn.png", dpi=300) # transparent background.
    print("draw plot: {}.venn.png".format("_".join(sample_name_list)))
    

if __name__ == '__main__':
    main()
