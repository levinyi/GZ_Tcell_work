#!/usr/bin/python3
import sys
import os
import re
import pandas as pd


def usage():
    print("""
        python {} <xls> <tpm file>
    Update:
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
            match_mut = re.search(r'[ATCG]', mapping_info.upper())
            wildtype = mapping_info.upper().count(".") + mapping_info.upper().count(",") + mapping_info.upper().count("$")
            if '+' in mapping_info:
                ins_type = mapping_info.upper().count("+")
            else:
                ins_type = 0
            if '*' in mapping_info:
                del_type = mapping_info.upper().count("*")
            else:
                del_type = 0
            if match_mut:
                mutatype = mapping_info.upper().count(match_mut.group())
            else:
                mutatype = 0
            #wildtype = int(reads_count) - mutatype
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

    data = pd.read_table(xls_file, sep="\t")
    header = list(data.columns.values)
    header.extend(["MutatedReads(RNA)", "Wild-typeReads(RNA)"])
    output_file.write("{}\n".format("\t".join(header)))
    
    for index, row in data.iterrows():
        Position_start = str(row["Start_Position"])
        Variant_Type = str(row['Variant_Type'])
        if Position_start in mpileup_dict:
            RNA_wild_reads = int(mpileup_dict[Position_start].get("W", 0))
            if Variant_Type == 'SNP':
                RNA_muta_reads = int(mpileup_dict[Position_start].get("M", 0))
            elif Variant_Type == 'DEL':
                RNA_muta_reads = int(mpileup_dict[Position_start].get("D", 0))
            elif Variant_Type == 'INS':
                RNA_muta_reads = int(mpileup_dict[Position_start].get("I", 0))
        else:
            print("position start not in mpileup dict. please check!")
            RNA_wild_reads = 0
            RNA_muta_reads = 0
        
        row["MutatedReads(RNA)"] = RNA_muta_reads
        row["Wild-typeReads(RNA)"] = RNA_wild_reads
        output_file.write("{}\n".format("\t".join([str(row[i]) for i in header])))
    
    output_file.close()


if __name__ == "__main__":
    main()