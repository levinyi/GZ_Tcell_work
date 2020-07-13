#!/usr/bin/python3
import sys
import os
import pandas as pd


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
            match = re.search(r'[ATCG]',mapping_info.upper())
            if match:
                match.group()
                mutatype = mapping_info.upper().count(match.group())
            else:
                mutatype = 0
            wildtype = int(reads_count) - mutatype 
            addtwodimdict(a_dict, pos, 'W', wildtype)
            addtwodimdict(a_dict, pos, 'M', mutatype)
    return a_dict

xls_file = sys.argv[1]
RNAseq_mpileup_file = sys.argv[2]

mpileup_dict = deal_mpileup(RNAseq_mpileup_file)

data = pd.read_table(xls_file, sep="\t")

contain_fields = [
        "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
        "Start_Position", "End_Position","Strand", "Variant_Classification", "Variant_Type", "Reference_Allele",
        "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "Genome_Change", "Annotation_Transcript", 
        "Transcript_Strand", "Transcript_Exon", "cDNA_Change", "Codon_Change", "Protein_Change",
        "Refseq_mRNA_Id","tumor_f", "t_alt_count", "t_ref_count", "n_alt_count", 
        "n_ref_count", "DP", "Wild-Type_Minigene","Mutated_Minigene",
    ]