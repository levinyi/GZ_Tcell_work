import sys
import os
import pandas as pd
import numpy as np


def usage():
    print("""
        python {} <xls> <tpm file>
Update:
        20220818    fix bugs.
        20220725    re-write the structure of the print function.
        20220719    re-write the function using pandas.
        20210410    fix bugs.
        20200703    Created.
    """.format(os.path.basename(sys.argv[0])))


def deal_tpm_file(tpm_file):
    adict = {}
    with open(tpm_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            transcript, tpm = line.split()
            transcript = transcript.split(".")[0]  # bug?
            adict[transcript] = tpm
    return adict


def main():
    if len(sys.argv) == 1:
        usage()
        sys.exit("Error: Two input files and one argment are needed!")

    if sys.argv[3] in ["gene_id", "transcript_id"]:
        arg = {"gene_id": "Hugo_Symbol", "transcript_id": "Annotation_Transcript"}
    else:
        sys.exit("Error: The third argment must be gene_id or transcript_id!")
    
    afile = sys.argv[1]
    tpm_file = sys.argv[2]
    output_file = open(os.path.basename(afile).rstrip("xls") + 'add.TPM.xls', "w")

    TPM_dict = deal_tpm_file(tpm_file)
 
    data = pd.read_table(afile, sep="\t", dtype=str)
    header = list(data.columns.values)
    header.append("TPM")
    output_file.write("{}\n".format("\t".join(header)))
    data = data.replace(np.nan, "")

    for index, row in data.iterrows():
        # remain the original data
        # if row["Chromosome"] =="":
        #     output_file.write("{}\n".format("\t".join([str(i) for i in row])))
        #     continue

        transcript = row[arg[sys.argv[3]]].split(".")[0]

        if transcript in TPM_dict:
            row["TPM"] = TPM_dict[transcript]
        else:
            row["TPM"] = "NA"
        output_file.write("{}\n".format("\t".join([str(x) for x in row])))
    output_file.close()


if __name__ == '__main__':
    main()
