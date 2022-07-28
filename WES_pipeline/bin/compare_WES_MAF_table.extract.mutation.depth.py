#!/usr/bin/python3
import sys
import os
import re
import pandas as pd


def usage():
    print("""
        python {} <xls> <mpileup1> <mpileup2>
    Update:
        20220217    Created.
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
            addtwodimdict(a_dict, pos, 'W', wildtype)
            addtwodimdict(a_dict, pos, 'M', mutatype)
            addtwodimdict(a_dict, pos, 'I', ins_type)
            addtwodimdict(a_dict, pos, 'D', del_type)
    return a_dict

xls_file = sys.argv[1]
RNAseq_mpileup_file1 = sys.argv[2]
RNAseq_mpileup_file2 = sys.argv[3]

mpileup_dict1 = deal_mpileup(RNAseq_mpileup_file1)
mpileup_dict2 = deal_mpileup(RNAseq_mpileup_file2)
# print(mpileup_dict1)

contain_fields = [
        "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
        "Start_Position", "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele", 
        "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "Genome_Change", "Annotation_Transcript", 
        "Transcript_Strand", "Transcript_Exon", "cDNA_Change", "Codon_Change", "Protein_Change",
        "Refseq_mRNA_Id", "tumor_f", "t_alt_count", "t_ref_count", "n_alt_count",
        "n_ref_count", "DP", "Mutated_Minigene", "Wild-Type_Minigene",
        "MutatedReads(RNA)", "Wild-typeReads(RNA)", "TPM",
        "muta_reads_in_Normal","wild_reads_in_Normal","muta_reads_in_Tumor","wild_reads_in_Tumor"
    ]

output_file = open(os.path.basename(xls_file).rstrip("xls")+'compare_Depth.xls', "w")
output_file.write("{}\n".format("\t".join(contain_fields)))

data = pd.read_table(xls_file, sep="\t")
for index, row in data.iterrows():
    Position_start = str(row["Start_Position"])
    Variant_Type = str(row['Variant_Type'])
    
    if Position_start in mpileup_dict1:
        wild_reads_in_Normal = int(mpileup_dict1[Position_start].get("W", 0))
        if Variant_Type == 'SNP':
            muta_reads_in_Normal = int(mpileup_dict1[Position_start].get("M", 0))
        elif Variant_Type == 'DEL':
            muta_reads_in_Normal = int(mpileup_dict1[Position_start].get("D", 0))
        elif Variant_Type == 'INS':
            muta_reads_in_Normal = int(mpileup_dict1[Position_start].get("I", 0))
    else:
        print("position start not in mpileup dict. please check! :{}".format(row))
        wild_reads_in_Normal = 0
        muta_reads_in_Normal = 0

    if Position_start in mpileup_dict2:
        wild_reads_in_Tumor = int(mpileup_dict2[Position_start].get("W", 0))
        if Variant_Type == 'SNP':
            muta_reads_in_Tumor = int(mpileup_dict2[Position_start].get("M", 0))
        elif Variant_Type == 'DEL':
            muta_reads_in_Tumor = int(mpileup_dict2[Position_start].get("D", 0))
        elif Variant_Type == 'INS':
            muta_reads_in_Tumor = int(mpileup_dict2[Position_start].get("I", 0))
    else:
        print("position start not in mpileup dict. please check!\n{} {} {}".format(row['Hugo_Symbol'],row["Chromosome"],row["Start_Position"]))
        wild_reads_in_Tumor = 0
        muta_reads_in_Tumor = 0
    
    output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
        row["Hugo_Symbol"], row["Entrez_Gene_Id"], row["Center"], row["NCBI_Build"], row["Chromosome"],
        row["Start_Position"], row["Strand"], row["Variant_Classification"], row["Variant_Type"], row["Reference_Allele"],
        row["Tumor_Seq_Allele1"], row["Tumor_Seq_Allele2"], row["dbSNP_RS"], row["Genome_Change"], row["Annotation_Transcript"],
        row["Transcript_Strand"], row["Transcript_Exon"], row["cDNA_Change"], row["Codon_Change"], row["Protein_Change"],
        row["Refseq_mRNA_Id"], row["tumor_f"], row["t_alt_count"], row["t_ref_count"], row["n_alt_count"], 
        row["n_ref_count"], row["DP"],row["Mutated_Minigene"],row["Wild-Type_Minigene"],row["MutatedReads(RNA)"],row["MutatedReads(RNA)"],row["TPM"],
        muta_reads_in_Normal, wild_reads_in_Normal,
        muta_reads_in_Tumor, wild_reads_in_Tumor
        ))
    
