import sys
import os
import pandas as pd

def usage():
    '''
    this script only used for XXX.exon.gene_name.counts.txt.TPM.txt file.

    usage:
    python merge.tpm.py /cygene2/work/Tumor_Expression_Data/*.TPM.txt >test.txt

    update:
    20210713. fix bug: GB001003-Tumor-pc-P6 key error.
    20210204. added usage description.
    2020xxxx. created.
    '''
    pass


def main():
    big_DataFrame = pd.DataFrame()
    for each in sys.argv[1:]:
        sample_name = os.path.basename(each).split(".")[0]
        df1 = pd.read_table(each, sep="\t", header=None, names=["Gene",sample_name])
        if big_DataFrame.empty :
            big_DataFrame = df1
        else:
            big_DataFrame = pd.merge(big_DataFrame, df1, on='Gene', how='outer')
    print(big_DataFrame.to_csv(sep="\t", index=False))


if __name__ == '__main__':
    main()
