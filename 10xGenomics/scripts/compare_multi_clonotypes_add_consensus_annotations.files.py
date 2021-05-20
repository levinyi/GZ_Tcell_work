import sys
import os
import pandas as pd


def usage():
    print('''
    usage:
    
    20210504    restructured.
    20200914    created.
    '''.format(os.path.basename(sys.argv[0])))


def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict


def deal_clonotypes_add_consensus_file(afile):
    pair_dict = {}
    with open(afile, "r") as f:
        for line in f:
            if line.startswith("clonotype_id"):
                continue
            line = line.rstrip("\n")
            c = line.split(",")
            clonotype_pair = c[4].split("*")[0].replace("/","") + c[10] + c[6].split("*")[0] + "_" + c[17].split("*")[0].replace("/","") + c[23] + c[19].split("*")[0]
            pair_dict[clonotype_pair] = c[-1]
    return pair_dict
def deal_clonotypes_add_consensus_files(afile):
    df = pd.read_csv(afile, sep=",")
    df['clonotype_pair'] = df['v_gene'].str.split("*").str[0].str.replace("/","") + df['cdr3'] + df['j_gene'].str.split("*").str[0].str.replace("/","")+ df['v_gene.1'].str.split("*").str[0].str.replace("/","")+ df['cdr3.1'] + df['j_gene.1'].str.split("*").str[0].str.replace("/","")
    print(df.head(10))
    return df


def main():
    data_dict = {}
    for each_file in sys.argv[1:]:
        file_name = os.path.basename(each_file)
        sample_name = file_name.split(".")[0]
        df = deal_clonotypes_add_consensus_files(each_file)
    data_dict[sample_name] = df
    '''
    big_dict = {}
    big_list = []
    sample_list = []
    for each_file in sys.argv[1:]:
        file_name = os.path.basename(each_file)
        sample_name = file_name.split(".")[0]
        sample_list.append(sample_name)
        pair_dict = deal_clonotypes_add_consensus_file(each_file)
        for each_pair in pair_dict:
            addtwodimdict(big_dict, each_pair, sample_name, pair_dict[each_pair])
            if each_pair not in big_list:
                big_list.append(each_pair)
    
    output  = open("{}.frequency.xls".format("_".join(sample_list))  ,"w")
    output.write("clonotype_pairs\t{}\n".format("\t".join(sample_list)))
    for each_pair in big_list:
        output.write("{}".format(each_pair))
        for each_sample in sample_list:
            if each_pair in big_dict:
                output.write("\t{}".format(big_dict[each_pair].get(each_sample, "NULL")))
            else:
                output.write("\tNULL")
        output.write("\n")

    output.close()
    '''


if __name__ == '__main__':
    main()
