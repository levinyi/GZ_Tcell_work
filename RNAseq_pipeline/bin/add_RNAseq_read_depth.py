#!/usr/bin/python3
import sys
import os
import re
import pandas as pd
import numpy as np


def usage():
    print("""
        python {} <xls> <tpm file>
    Update:
        20230116    fix mpileup count error bug.
        20220818    fix format bug.
        20220725    re-write the function to print NAN to empty.
        20210410    add parmaters.
        20200703    Created.
    """.format(os.path.basename(sys.argv[0])))

def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict


def deal_mpileup(RNAseq_mpileup_file):
    a_dict = {}
    with open(RNAseq_mpileup_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            chrom, pos, mapped_base, reads_count, mapping_info, mapping_qul = line.split("\t")
            # remove ^and qualitybase. and remove $. or <<><<<>>><><>
            mapping_info = re.sub(r'\^\S', '', mapping_info)
            mapping_info = re.sub(r'\$',   '', mapping_info)
            mapping_info = re.sub(r'>+|<+','', mapping_info)

            if '+' in mapping_info:
                ins_type = mapping_info.upper().count("+")
                mutatype = ins_type
            else:
                ins_type = 0

            if '*' in mapping_info:
                del_type = mapping_info.upper().count("*")
                mutatype = del_type
            else:
                del_type = 0
            
            wildtype = mapping_info.upper().count(".") + mapping_info.upper().count(",")
            mutatype = len(mapping_info) - ins_type - del_type - wildtype
            
            addtwodimdict(a_dict, pos, 'W', wildtype)
            addtwodimdict(a_dict, pos, 'M', mutatype)
            addtwodimdict(a_dict, pos, 'I', ins_type)
            addtwodimdict(a_dict, pos, 'D', del_type)
    return a_dict


def main():
    xls_file = sys.argv[1]
    RNAseq_mpileup_file = sys.argv[2]

    mpileup_dict = deal_mpileup(RNAseq_mpileup_file)

    output_file = open(os.path.basename(xls_file).rstrip("xls")+'RNA_Depth.xls', "w")
    output_file2 = open(os.path.basename(xls_file).rstrip("xls")+'RNA_Depth.Error.Position.txt', "w")
    output_file2.write("## Position not in mpileup file. please check!\n")
    output_file2.write("Hugo_Symbol\tChromosome\tStart_Position\t\End_Position\tflag\tTPM\n")
    data = pd.read_table(xls_file, sep="\t", dtype=str)
    header = list(data.columns.values)
    header.extend(["MutatedReads(RNA)", "Wild-typeReads(RNA)"])
    output_file.write("{}\n".format("\t".join(header)))
    
    data = data.replace(np.nan, "")
    # print(data.head())
    
    for index, row in data.iterrows():
        # if chromosome is empty, print the row
        if row['Chromosome'] == "":
            output_file.write("{}\n".format("\t".join([str(i) for i in row])))
            continue

        Position_start = str(row["Start_Position"])
        Variant_Type = str(row['Variant_Type'])
        if Position_start in mpileup_dict:
            RNA_wild_reads = int(mpileup_dict[Position_start].get("W", 0))
            if Variant_Type == 'SNP' or Variant_Type == 'DNP':
                RNA_muta_reads = int(mpileup_dict[Position_start].get("M", 0))
            elif Variant_Type == 'DEL':
                RNA_muta_reads = int(mpileup_dict[Position_start].get("D", 0))
            elif Variant_Type == 'INS':
                RNA_muta_reads = int(mpileup_dict[Position_start].get("I", 0))
        else:
            if str(row['End_Position']) in mpileup_dict:
                flag="End_Position_in_mpileup"
            elif str(int(Position_start) -1) in mpileup_dict:
                flag = "start_position-1_in_mpileup"
            else:
                flag = "both_Not_in_mpileup"
            output_file2.write("{}\t{}\t{}\t{}\t{}\t{}\t\n".format(row["Hugo_Symbol"], row["Chromosome"], Position_start, row['End_Position'], flag, row["TPM"])) # this is a bug, should be fixed in the future.
            RNA_wild_reads = "NA"
            RNA_muta_reads = "NA"
        
        row["MutatedReads(RNA)"] = RNA_muta_reads
        row["Wild-typeReads(RNA)"] = RNA_wild_reads
        output_file.write("{}\n".format("\t".join([str(row[i]) for i in header])))
    output_file.close()


if __name__ == "__main__":
    main()
