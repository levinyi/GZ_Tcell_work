import sys,os
import pandas as pd
from venn import venn
from matplotlib import pyplot as plt



def main():
    venn_dict = {} # for venn diagram
    big_DataFrame = pd.DataFrame()
    sample_name_list = []
    for each in sys.argv[1:]:
        sample_name = os.path.basename(each).split("_")[0]
        print("read file: {}".format(sample_name))
        sample_name_list.append(sample_name)
        # df1 = pd.read_table(each, sep="\t", header=None, names=["Gene",sample_name])
        # read table 
        df1 = pd.read_csv(each, sep=",", header=0)
        # print(df1.head())
        # add new clonotype_tra_id,clonotype_trb_id,clonotype_pair_id
        df1['clonotype_tra_id'] = df1[['v_gene', 'cdr3', 'j_gene']].apply('_'.join, axis=1)
        df1['clonotype_trb_id'] = df1[['v_gene.1', 'cdr3.1', 'j_gene.1']].apply('_'.join, axis=1)
        df1['clonotype_pair_id'] = df1[['clonotype_tra_id','clonotype_trb_id']].apply("_".join, axis=1)
        # select data clonotype_tra_id,clonotype_trb_id,clonotype_pair_id, freq
        df1 = df1[['clonotype_pair_id','clonotype_id','proportion']]
        # print(df1.head())
        # add suffixes name for each sample:
        df1 = df1.rename(columns={'clonotype_id':'clonotype_id_'+sample_name,'proportion':'proportion_'+sample_name},)
        if big_DataFrame.empty :
            big_DataFrame = df1
        else:
            big_DataFrame = pd.merge(big_DataFrame, df1, on='clonotype_pair_id', how='outer')
        # to venn list:
        venn_dict.setdefault(sample_name, set(df1['clonotype_pair_id'].tolist()))
    big_DataFrame.to_csv("_".join(sample_name_list)+ ".xls", sep="\t", index=False)
    # print(big_DataFrame.head())
    # print(venn_dict)
    print("write xls file: {}.xls".format("_".join(sample_name_list)))

    venn(venn_dict)
    plt.savefig("_".join(sample_name_list)+".venn.png",dpi=300)
    print("draw plot: {}.venn.png".format("_".join(sample_name_list)))


if __name__ == '__main__':
    main()
