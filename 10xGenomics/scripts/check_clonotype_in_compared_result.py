import sys
import pandas as pd
import argparse


def _argparse():
    parser = argparse.ArgumentParser(description="This is metagenomics pipeline")
    parser.add_argument('-f1', '--clonotype_file', action='store', dest='clonotype_file', default='config.txt', help = "a list file.")
    parser.add_argument('-f2', '--total_file',     action='store', dest='total_file', default='total.clonotypes.frequency.csv', help = "total file")
    parser.add_argument('-l', '--left_column_name', action='store', dest='left_column_name', help = "left side column name")
    parser.add_argument('-r', '--right_column_name', action='store', dest='right_column_name', help = "right side column name")
    parser.add_argument('-o', '--output', action='store', dest='output_file_name',default='clonotype.checked.out.list.csv', help="output file name")
    return parser.parse_args()

def main():
    parser = _argparse()

    left_column_name = parser.left_column_name
    right_column_name = parser.right_column_name
    total_file = parser.total_file
    clonotype_file = parser.clonotype_file

    # deal with clonotype file, return a list.
    clonotype_data = pd.read_table(clonotype_file, header=None, names=['clonotype'])
    clonotype_list = clonotype_data['clonotype'].tolist()
    
    # deal with total file .
    total_data = pd.read_table(total_file, sep=",")
    rmdup_total_data = total_data.drop_duplicates([left_column_name])
    select = rmdup_total_data.loc[rmdup_total_data[left_column_name].isin(clonotype_list),[left_column_name,right_column_name]]
    if parser.output_file_name:
        select.to_csv(parser.output_file_name, sep="\t", index=False)
    else:
        select.to_csv("clonotype.checked.out.list.csv", sep="\t", index=False)
        

    left_list = select.loc[select[right_column_name].isna(),[left_column_name]][left_column_name].tolist()
    middle_list = select.loc[select[right_column_name].notnull(),[left_column_name]][left_column_name].tolist()

    print("left: [{}]\nTotal number in left: {}".format(",".join(left_list),len(left_list)))
    print("middle: [{}]\nTotal number in middle: {}".format(",".join(middle_list), len(middle_list)))

if __name__ == '__main__':
    main()
