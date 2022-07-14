from curses import raw
from dis import dis
from functools import partial
from hashlib import new
from operator import ne
import sys
import os
import re
from tokenize import endpats
from tracemalloc import stop
from turtle import st
from webbrowser import Chrome
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import warnings
from Bio import BiopythonWarning
from regex import R
from tables import Unknown
warnings.simplefilter('ignore', BiopythonWarning)


def usage():
    print("""
Usage:
    python {0} <maf file > <transcription file> <output file>
Example:
    python {0} 
Notes:
    This script is used for extract minigene from a cds file by a position.
    maf file without annotation line (drop lines which contains #).
    Function:
        1. filter useless sites.
            Mandatory fields: Missense_Mutation, Frame_Shift_Ins, Frame_Shift_Del, Nonesense_Mutation, Splice_site, In_Frame_Del, In_Frame_Ins.
            Silent_Mutation may used for check mutation position in transcription sequence.
        2. filter useless annotation items. 
            Mandatory fields: Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Variant_Type and Tumor_Sample_Barcode.
        3. standardization output format to a xls format.

Updates:
    20220711    Add 5 columns to output file: doNotsyn,MutMG_ID_for_order, Mut_AA29_for_order, wtMG_ID_for_order, wt_AA29_for_order.
    20220217    Fix bugs: DtypeWarning: Columns (56,88,101,103,107,108,111,147) have mixed types.Specify dtype option on import or set low_memory=False. [at pandas read file]
    20210423    Fix bugs: if mutation occurs in Mitochondrial, use coding_dna.translate(table="Vertebrate Mitochondrial")
    20200709    Add "End_Position" field to result.
    20200702    Updated. fix cDNA_Change bugs: c.11119_11120CC>AT.
    20200628    Created.
    """.format(os.path.basename(sys.argv[0])))


def extract_minigene(raw_seq, new_seq, Protein_Change, Chromosome, Protein_Change2=None):  
    if Protein_Change2:
        start1 = int(re.findall(r'\d+', Protein_Change)[0])
        start2 = int(re.findall(r'\d+', Protein_Change2)[0])

        # sort Protein_Change and Protein_Change2 by position.
        if start1 > start2:
            Protein_Change, Protein_Change2 = Protein_Change2, Protein_Change
    # if Chromosome is Mitochondrial. Use Seq(seq).translate(table="Vertebrate Mitochondrial") otherwise it will raise error.
    if Chromosome == "MT":
        raw_seq_aa = str(Seq(raw_seq).translate(table="Vertebrate Mitochondrial"))
        new_seq_aa = str(Seq(new_seq).translate(table="Vertebrate Mitochondrial"))
    else:
        raw_seq_aa = str(Seq(raw_seq).translate())
        new_seq_aa = str(Seq(new_seq).translate())

    aa_position = int(re.findall(r'\d+', Protein_Change)[0]) # if aa_position occurs before 14aa, it means the mutation is in the first 14aa.
    if aa_position < 14:
        start_position = 0
        end_position = 29
    else:
        start_position = aa_position - 15
        end_position = aa_position + 14
    
    # store the minigene in a list.
    raw_minigene_list = []
    new_minigene_list = []
    
    # print("1st,new_minigene_list: {}".format(new_minigene_list))
    # if new_seq_aa contain "*", it means the stop codon occured. and find the position of "*" in new_seq_aa.
    if "*" in new_seq_aa and not new_seq_aa.endswith("*"):
        stop_codon_position = new_seq_aa.find("*") # find first stop codon position in new_seq_aa.
        while stop_codon_position - end_position > 0:
            print("stop_codon_position: {}, end_position: {} {}".format(stop_codon_position, end_position, stop_codon_position - end_position))
            a = new_seq_aa[start_position:end_position]
            b = raw_seq_aa[start_position:end_position]
            if a not in new_minigene_list:
                new_minigene_list.append(a)
            if b not in raw_minigene_list:
                raw_minigene_list.append(b)
            start_position += 15
            end_position += 15
            # print("While : new_minigene_list: {}".format(new_minigene_list))
        else:
            print("only KAZN gene can print this line {}".format(Protein_Change))
            a = new_seq_aa[stop_codon_position-29: stop_codon_position]
            b = raw_seq_aa[stop_codon_position-29: stop_codon_position]
            if a not in new_minigene_list:
                new_minigene_list.append(a)
            if b not in raw_minigene_list:
                raw_minigene_list.append(b)
        print("stop_codon_position: {}".format(stop_codon_position))
        print("Out new_minigene_list: {}".format(new_minigene_list))
    else:
        new_minigene_list.append(new_seq_aa[start_position:end_position])
        raw_minigene_list.append(raw_seq_aa[start_position:end_position])
        # print("new_minigene_list: {}".format(new_minigene_list))
        # print("raw_minigene_list: {}".format(raw_minigene_list))
    
    if Protein_Change2:
        aa_position2 = int(re.findall(r'\d+', Protein_Change2)[0])
        if aa_position2 < 14:
            start_position = 0
            end_position = 29
        else:
            start_position = aa_position2 - 15
            end_position = aa_position2 + 14
        print("pro1:{},pro2:{}, aa_position2: {}, start_position: {}, end_position: {}".format(Protein_Change, Protein_Change2, aa_position2, start_position, end_position))
        print("raw_seq_aa: {}".format(raw_seq_aa))
        print("new_seq_aa: {}".format(new_seq_aa))
        if "*" in new_seq_aa and not new_seq_aa.endswith("*"):
            stop_codon_position = new_seq_aa.find("*")
            # print("stop_codon_position: {}, start: {}, end :{} ,{}, aa_position: {}".format(stop_codon_position,start_position, end_position,stop_codon_position-end_position, aa_position2))
            while stop_codon_position - end_position > 0:
                a = new_seq_aa[start_position:end_position]
                b = raw_seq_aa[start_position:end_position]
                if a not in new_minigene_list:
                    new_minigene_list.append(a)
                if b not in raw_minigene_list:
                    raw_minigene_list.append(b)
                start_position += 15
                end_position += 15
                # print("Protein_Change2 While: new_minigene_list: {}".format(new_minigene_list))
            else:
                a = new_seq_aa[stop_codon_position-29: stop_codon_position]
                b = raw_seq_aa[stop_codon_position-29: stop_codon_position]
                if a not in new_minigene_list:
                    new_minigene_list.append(a)
                if b not in raw_minigene_list:
                    raw_minigene_list.append(b)
                # print("Protein_Change2 else: new_minigene_list: {}".format(new_minigene_list))
        else:
            new_minigene_list.append(new_seq_aa[start_position: end_position])
            raw_minigene_list.append(raw_seq_aa[start_position: end_position])
            # print("Protein_Change2 IFelse: new_minigene_list: {}".format(new_minigene_list))
    return raw_minigene_list, new_minigene_list


def get_cDNA_Change_Info(cDNA_Change):
    # Examples for cDNA_Change:
    # c.1518_1519insAAACAGACCA 
    # c.2953_2972delCTAAATCACACTCCTGTATC , c.4113delG 
    # c.3335A>G , c.11119_11120CC>AT
    cDNA_Change = cDNA_Change.lstrip("c.")
    if "ins" in cDNA_Change:
        position_components, middle_seq = cDNA_Change.split("ins") 
        start, end = position_components.split("_")
        return "ins", int(start), int(end), middle_seq
    elif "del" in cDNA_Change:
        position_components, middle_seq = cDNA_Change.split("del")
        # deletion_length = len(middle_seq) # interger.
        if "_" in position_components:
            start, end = position_components.split("_")
            start = int(start)
            end = int(end)
        else:
            start = position_components
            start = int(start)
            end = int(start)
        return "del", start, end, middle_seq
    elif ">" in cDNA_Change:
        # c.11119_11120CC>AT
        position_components, after_mutation = cDNA_Change.split(">")
        starts = re.findall(r'\d+', position_components)
        if len(starts) == 2:
            start = int(starts[0])
            end = int(starts[1])
        elif len(starts) == 1:
            start = int(starts[0])
            end = int(starts[0])
        else:
            print("Ops! Unexpected format in cDNA_Change: {}".format(cDNA_Change))
        return ">", start, end, after_mutation
    else:
        print("Ops! Unexpected Condition in cDNA_Change: {}".format(cDNA_Change))
        return "unknown"

def single_position_fix(raw_seq, cDNA_Change):
    c1_type, c1_start, c1_end, c1_seq = get_cDNA_Change_Info(cDNA_Change)
    if c1_type == "ins":
        full_seq = raw_seq[:c1_start] + c1_seq + raw_seq[c1_start:]
    elif c1_type == "del":
        full_seq = raw_seq[:c1_start-1] + raw_seq[c1_end:]
    elif c1_type == ">":
        full_seq = raw_seq[:c1_start-1] + c1_seq + raw_seq[c1_start:]
    return full_seq

def multi_position_fix(raw_seq, cDNA_Change1, cDNA_Change2):
    start1 = int(re.findall(r'\d+', cDNA_Change1)[0])
    start2 = int(re.findall(r'\d+', cDNA_Change2)[0])

    # sort cDNA_Change1 and cDNA_Change2 by position.
    if start1 > start2:
        cDNA_Change1, cDNA_Change2 = cDNA_Change2, cDNA_Change1
    
    # get the Type of cDNA_Change1 and cDNA_Change2.
    c1_type, c1_start, c1_end, c1_seq = get_cDNA_Change_Info(cDNA_Change1)
    c2_type, c2_start, c2_end, c2_seq = get_cDNA_Change_Info(cDNA_Change2)
    # print("cDNA_Change1: {}, c1_type: {}, c1_start: {}, c1_end: {}, c1_seq: {}".format(cDNA_Change1, c1_type, c1_start, c1_end, c1_seq))
    # print("cDNA_Change2: {}, c2_type: {}, c2_start: {}, c2_end: {}, c2_seq: {}".format(cDNA_Change2, c2_type, c2_start, c2_end, c2_seq))
    
    if c1_type == "ins":
        partial_seq = raw_seq[:c1_start] + c1_seq + raw_seq[c1_start:c2_start]
    elif c1_type == "del":
        partial_seq = raw_seq[:c1_start-1] + raw_seq[c1_end: c2_start]
    elif c1_type == ">":
        partial_seq = raw_seq[:c1_start-1] + c1_seq + raw_seq[c1_start:c2_start]
    # print("partial_seq: {}".format(partial_seq))
    
    if c2_type == "ins":
        full_seq = partial_seq + c2_seq + raw_seq[c2_start:]
    elif c2_type == "del":
        full_seq = partial_seq[:-1] + raw_seq[c2_end:]
    elif c2_type == ">":
        full_seq = partial_seq[:-1] + c2_seq + raw_seq[c2_start:]
    # print("full_seq: {}".format(full_seq))
    return full_seq
            
def save_minigene(output_file, raw_minigene_list, new_minigene_list, row, gene_id, cDNA_Change, Protein_Change):
    for index, raw_minigene in enumerate(raw_minigene_list):
        if new_minigene_list[index] == raw_minigene:
            doNotSyn = "doNotSyn"
        else:
            doNotSyn = ""
        
        Mut_AA29_for_order = new_minigene_list[index]
        WT_AA29_for_order = raw_minigene
        
        if index+1 > 1:
            MutMG_ID_for_order = gene_id+"-"+ Protein_Change.lstrip("p.") + str(index+1)
            wtMG_ID_for_order = gene_id + "-" + str(re.findall(r'\d+', Protein_Change)[0]) + "WT"+str(index+1)
        else:
            MutMG_ID_for_order = gene_id+"-"+ Protein_Change.lstrip("p.")
            wtMG_ID_for_order = gene_id + "-" + str(re.findall(r'\d+', Protein_Change)[0]) + "WT"
        
        if "*" in new_minigene_list[index]:
            stop_position = new_minigene_list[index].index("*")
            add_GSG_number = 29 - stop_position   # 29 is the length of the amino acid.
            add_GSG_aa = "GSG" * add_GSG_number
            Mut_AA29_for_order = new_minigene_list[index][:stop_position] + add_GSG_aa[:add_GSG_number]  # add GSG to the end of the amino acid.    
        if "*" in raw_minigene:
            stop_position = raw_minigene.index("*")
            add_GSG_number = 29 - stop_position   # 29 is the length of the amino acid.
            add_GSG_aa = "GSG" * add_GSG_number
            WT_AA29_for_order = raw_minigene[:stop_position] + add_GSG_aa[:add_GSG_number]
            
        output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            row["Hugo_Symbol"], row["Entrez_Gene_Id"], row["Center"], row["NCBI_Build"], row["Chromosome"],  # 1-5
            row["Start_Position"], row["End_Position"], row["Strand"], row["Variant_Classification"], row["Variant_Type"],  # 6-10
            row["Reference_Allele"], row["Tumor_Seq_Allele1"], row["Tumor_Seq_Allele2"], row["Tumor_Sample_Barcode"], row["Matched_Norm_Sample_Barcode"],   # 11-15
            row["dbSNP_RS"], row["Genome_Change"], row["Annotation_Transcript"], row["Transcript_Strand"], row["Transcript_Exon"],  # 16-20
            cDNA_Change, row["Codon_Change"], row["Protein_Change"], row["Refseq_mRNA_Id"], row["tumor_f"],     # 21-25
            row["t_alt_count"], row["t_ref_count"], row["n_alt_count"], row["n_ref_count"], row["DP"],  # 26-30
            new_minigene_list[index], raw_minigene, doNotSyn, MutMG_ID_for_order, Mut_AA29_for_order,  # 31-35
            wtMG_ID_for_order, WT_AA29_for_order # 36-37
            ) 
        )


def main():
    cds_file = open(sys.argv[1], "r")
    maf_file = sys.argv[2]
    output_file = open(sys.argv[3], "w") 

    adict = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        record_id = str(record.id).split("|")[0].split(".")[0]
        if record_id.endswith("PAR_Y"): # skip xxx.*PAR_Y
            continue
        adict[record_id] = str(record.seq)
    cds_file.close()

    # read MAF file as a dataframe and first line stored as header.
    maf_data = pd.read_table(maf_file, sep="\t", dtype='str').fillna(value="NA")
    contain_fields = [
        "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", # 1-5
        "Start_Position", "End_Position","Strand", "Variant_Classification", "Variant_Type",  # 6-10
        "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",  # 11-15
        "dbSNP_RS", "Genome_Change", "Annotation_Transcript", "Transcript_Strand", "Transcript_Exon",  # 16-20
        "cDNA_Change", "Codon_Change", "Protein_Change", "Refseq_mRNA_Id","tumor_f", # 21-25
        "t_alt_count", "t_ref_count", "n_alt_count", "n_ref_count", "DP",  # 26-30
        "Mutated_Minigene", "Wild-Type_Minigene", "doNotSyn","MutMG_ID_for_order","Mut_AA29_for_order", # 31-35
        "wtMG_ID_for_order","WT_AA29_for_order", # 36-37
    ]
    output_file.write("{}\n".format("\t".join(contain_fields)))

    saved_Variant_Classification = {
        "Nonsense_Mutation": 1,
        "Missense_Mutation": 1,
        "Frame_Shift_Del": 1,
        "Frame_Shift_Ins": 1,
        "In_Frame_Del": 1, 
        "In_Frame_Ins": 1,
        "Splice_Site": 1,
        "Translation_Start_Site": 1,
        "Nonstop_Mutation":1,
        "Nonstop_Mutation": 1,
        # "Silent": 1,
    }

    # some basic information about the MAF file.
    useful_mutation = 0
    Unknown_Gene = 0
    unknown_type_count = 0

    # store the gene_name and mutation info into a dictionary.
    gene_dict = {} # key: gene_name, value: list of mutation info(whole line)
    for index, row in maf_data.iterrows():
        gene_id = row["Hugo_Symbol"]
        var_type = row["Variant_Classification"]
        if gene_id != "Unknown":
            if var_type in saved_Variant_Classification and row['Protein_Change'] != "NA":
                useful_mutation += 1
                gene_dict.setdefault(gene_id,[]).append(row)
            else:
                unknown_type_count += 1
        else:
            Unknown_Gene += 1
    print("Total number of useful mutations: {}".format(useful_mutation))
    print("Total number of Unknown_type_count: {}".format(unknown_type_count))
    print("Total number of Unknown_Gene: {}".format(Unknown_Gene))
    print("Total number of uniq genes: {}".format(len(gene_dict)))


    # start to extract minigene for each gene(line in the MAF file).
    for gene_id, rows in gene_dict.items():
        cDNA_Change = rows[0]["cDNA_Change"]
        Protein_Change = rows[0]["Protein_Change"]
        raw_seq = adict[rows[0]["Annotation_Transcript"].split(".")[0]]
        Chromosome = rows[0]["Chromosome"]

        # print(rows[0]['Hugo_Symbol'], cDNA_Change, rows[0]['Annotation_Transcript'], Protein_Change, Chromosome, rows[0]['Variant_Classification'])
        if len(rows) == 1: # only one mutation in this gene
            new_seq = single_position_fix(raw_seq, cDNA_Change)
            raw_minigene_list, new_minigene_list = extract_minigene(raw_seq, new_seq, Protein_Change, Chromosome)
            save_minigene(output_file, raw_minigene_list, new_minigene_list, rows[0], gene_id, cDNA_Change, Protein_Change)
        else: # multiple mutations in this gene
            print("Multiple mutations in this gene: {}".format(gene_id))
            positions = [row["Start_Position"] for row in rows] # get all the positions of this gene
            if len(positions) == 2: # only two mutations in this gene
                if abs(int(positions[0]) - int(positions[1])) > 45: # if the distance between two mutations is larger than 45, they are independent event.
                    # deal with the first mutation
                    new_seq = single_position_fix(raw_seq, cDNA_Change)
                    raw_minigene_list, new_minigene_list = extract_minigene(raw_seq, new_seq, Protein_Change, Chromosome)
                    # print("raw_minigene list: {}".format(raw_minigene_list))
                    # print("new_minigene list: {}".format(new_minigene_list))
                    save_minigene(output_file, raw_minigene_list, new_minigene_list, rows[0], gene_id, cDNA_Change, Protein_Change)
                    # deal with the second mutation
                    cDNA_Change2 = rows[1]["cDNA_Change"]
                    print("they are not disturbed by position {} {} {} {}".format(positions[0], positions[1],cDNA_Change, cDNA_Change2))
                    Protein_Change2 = rows[1]["Protein_Change"]
                    new_seq = single_position_fix(raw_seq, cDNA_Change2)
                    raw_minigene_list, new_minigene_list = extract_minigene(raw_seq, new_seq, Protein_Change2, Chromosome)
                    # print("raw_minigene list: {}".format(raw_minigene_list))
                    # print("new_minigene list: {}".format(new_minigene_list))
                    save_minigene(output_file, raw_minigene_list, new_minigene_list, rows[0], gene_id, cDNA_Change2, Protein_Change2)
                else: # if the distance between two mutations is less than 45bp, they are disturbed by position. they are simultaneous event.
                    cDNA_Change2 = rows[1]["cDNA_Change"] # c.2666_2667insT
                    print("they are disturbed by position {} {} {} {}".format(positions[0], positions[1], cDNA_Change, cDNA_Change2))
                    Protein_Change2 = rows[1]["Protein_Change"]
                    new_seq = multi_position_fix(raw_seq, cDNA_Change, cDNA_Change2)
                    
                    raw_minigene_list, new_minigene_list = extract_minigene(raw_seq, new_seq, Protein_Change, Chromosome, Protein_Change2)
                    print("raw_minigene list: {}".format(raw_minigene_list))
                    print("new_minigene list: {}".format(new_minigene_list))
                    # bug: if there are multiple mutations in the same gene, the first one will be saved
                    save_minigene(output_file, raw_minigene_list, new_minigene_list, rows[0], gene_id, cDNA_Change, Protein_Change)
                    
                    ######## fix bug: if there are multiple mutations in the same gene, the second one will be saved
                    row = rows[1]
                    output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        row["Hugo_Symbol"], row["Entrez_Gene_Id"], row["Center"], row["NCBI_Build"], row["Chromosome"],  # 1-5
                        row["Start_Position"], row["End_Position"], row["Strand"], row["Variant_Classification"], row["Variant_Type"],  # 6-10
                        row["Reference_Allele"], row["Tumor_Seq_Allele1"], row["Tumor_Seq_Allele2"], row["Tumor_Sample_Barcode"], row["Matched_Norm_Sample_Barcode"],   # 11-15
                        row["dbSNP_RS"], row["Genome_Change"], row["Annotation_Transcript"], row["Transcript_Strand"], row["Transcript_Exon"],  # 16-20
                        cDNA_Change, row["Codon_Change"], row["Protein_Change"], row["Refseq_mRNA_Id"], row["tumor_f"],     # 21-25
                        row["t_alt_count"], row["t_ref_count"], row["n_alt_count"], row["n_ref_count"], row["DP"],  # 26-30
                        ) 
                )
            else:
                print("more than 3 positions in this gene: {}. please manually check.".format(gene_id))
         
    output_file.close()
    print("Finished extract minigene!")

if __name__ == '__main__':
    main()
