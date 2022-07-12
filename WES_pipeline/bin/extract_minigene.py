from dis import dis
from functools import partial
import sys
import os
import re
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


def get_codon_minigene(Chromosome, seq, cDNA_Change, Protein_Change):
    # Four examples of cDNA_Change item:
    # c.1518_1519insAAACAGACCA
    # c.550_551insGGGGCCGCC
    # c.2953_2972delCTAAATCACACTCCTGTATC 
    # c.4113delG 
    # c.3335A>G , c.11119_11120CC>AT
    cDNA_Change = cDNA_Change.lstrip("c.")

    if "ins" in cDNA_Change:
        position_components, insertion_seq = cDNA_Change.split("ins") 
        # c.204_205insT   c.(205-207)tttfs        p.E70fs 
        start, end = position_components.split("_")
        start = int(start)
        end = int(end)
        new_seq = seq[:start] + insertion_seq + seq[end-1:]
    elif "del" in cDNA_Change:
        position_components, deletion_seq = cDNA_Change.split("del")
        deletion_length = len(deletion_seq) # interger.
        if "_" in position_components:
            start, end = position_components.split("_")
            start = int(start)
            end = int(end)
            new_seq = seq[:start-1] + seq[end:]
        else:
            start = position_components
            start = int(start)
            new_seq = seq[:start-1] +seq[start:]
    elif ">" in cDNA_Change:
        # c.26G   >   C
        # c.11119_11120CC>AT
        position_components, after_mutation = cDNA_Change.split(">")
        starts = re.findall(r'\d+', position_components)
        if len(starts) == 2:
            start = int(starts[0])
            end = int(starts[1])
            new_seq = seq[:start-1] + after_mutation + seq[end:]
        elif len(starts) == 1:
            start = int(starts[0])
            new_seq = seq[:start-1] + after_mutation + seq[start:]
        else:
            print("Ops! Unexpected format in cDNA_Change: {}".format(cDNA_Change))
    else:
        print("Ops! Unexpected Condition in cDNA_Change: {}".format(cDNA_Change))

    # Four examples of Protein_Change item:
    # p.C725*
    # p.I373fs
    # p.299_300insGLD
    # p.LQREKLQREK1479del
    result = re.findall(r'\d+', Protein_Change)
    aa_position = int(result[0])

    # if Chromosome is Mitochondrial. Use Seq(seq).translate(table="Vertebrate Mitochondrial") otherwise it will raise error.
    if Chromosome == "MT":
        seq_aa = Seq(seq).translate(table="Vertebrate Mitochondrial")
        new_seq_aa = Seq(new_seq).translate(table="Vertebrate Mitochondrial")
    else:
        seq_aa = Seq(seq).translate()
        new_seq_aa = Seq(new_seq).translate()
    
    # if new_seq_aa contain "*", it means the stop codon occured. and find the position of "*" in new_seq_aa.
    if "*" in new_seq_aa:
        stop_codon_position = new_seq_aa.find("*")

    # if aa_position occurs before 14aa, it means the mutation is in the first 14aa.
    if aa_position < 14:
        old_minigene = seq_aa[: aa_position] + seq_aa[aa_position: aa_position + 14]
        new_minigene = new_seq_aa[: aa_position] + new_seq_aa[aa_position: aa_position + 14]
    else:
        old_minigene = seq_aa[aa_position -1 - 14: aa_position] + seq_aa[aa_position: aa_position + 14]
        new_minigene = new_seq_aa[aa_position -1 - 14: aa_position] + new_seq_aa[aa_position: aa_position + 14]
    return old_minigene, new_minigene


def get_cDNA_Change_Info(cDNA_Change):
    # Four examples of cDNA_Change item:
    # c.1518_1519insAAACAGACCA
    # c.550_551insGGGGCCGCC
    # c.2953_2972delCTAAATCACACTCCTGTATC 
    # c.4113delG 
    # c.3335A>G , c.11119_11120CC>AT
    cDNA_Change = cDNA_Change.lstrip("c.")
    if "ins" in cDNA_Change:
        position_components, middle_seq = cDNA_Change.split("ins") 
        # c.204_205insT   c.(205-207)tttfs        p.E70fs 
        start, end = position_components.split("_")
        start = int(start)
        end = int(end)
        return "ins", start, end, middle_seq
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

def multi_position_fix(raw_seq, cDNA_Change1, cDNA_Change2):
    result = re.findall(r'\d+', cDNA_Change1)
    start1 = int(result[0])
    result = re.findall(r'\d+', cDNA_Change2)
    start2 = int(result[0])
    # sort cDNA_Change1 and cDNA_Change2 by position.
    if start1 > start2:
        cDNA_Change1, cDNA_Change2 = cDNA_Change2, cDNA_Change1
    
    # get the Type of cDNA_Change1 and cDNA_Change2.
    c1_type, c1_start, c1_end, c1_seq = get_cDNA_Change_Info(cDNA_Change1)
    c2_type, c2_start, c2_end, c2_seq = get_cDNA_Change_Info(cDNA_Change2)
    print("cDNA_Change1: {}, c1_type: {}, c1_start: {}, c1_end: {}, c1_seq: {}".format(cDNA_Change1, c1_type, c1_start, c1_end, c1_seq))
    print("cDNA_Change2: {}, c2_type: {}, c2_start: {}, c2_end: {}, c2_seq: {}".format(cDNA_Change2, c2_type, c2_start, c2_end, c2_seq))
    
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
            

def main():
    if len(sys.argv) == 1:
        usage()
        sys.exit("Error: Wrong input file")

    cds_file = open(sys.argv[1], "r")
    maf_file = sys.argv[2]
    output_file = open(sys.argv[3], "w") 

    adict = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        record_id = str(record.id).split("|")[0]
        adict[record_id] = str(record.seq)
    cds_file.close()

    # read maf file as a dataframe and first line stored as header.
    maf_data = pd.read_table(maf_file, sep="\t", dtype='str').fillna(value="NA")
    contain_fields = [
        "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
        "Start_Position", "End_Position","Strand", "Variant_Classification", "Variant_Type", "Reference_Allele",
        "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode","dbSNP_RS", 
        "Genome_Change", "Annotation_Transcript", 
        "Transcript_Strand", "Transcript_Exon", "cDNA_Change", "Codon_Change", "Protein_Change",
        "Refseq_mRNA_Id","tumor_f", "t_alt_count", "t_ref_count", "n_alt_count", 
        "n_ref_count", "DP", "Mutated_Minigene", "Wild-Type_Minigene",
        "doNotSyn","MutMG_ID_for_order","Mut_AA29_for_order","wtMG_ID_for_order","WT_AA29_for_order",
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

    gene_dict = {}
    total_number_of_mutations = 0
    Unknown_Gene = 0
    unknown_type_count = 0
    for index, row in maf_data.iterrows():
        gene_id = row["Hugo_Symbol"]
        var_type = row["Variant_Classification"] = row["Variant_Classification"]
        if gene_id != "Unknown":
            if var_type in saved_Variant_Classification:
                total_number_of_mutations += 1
                gene_dict.setdefault(gene_id,[]).append(row)
            else:
                unknown_type_count += 1
        else:
            Unknown_Gene += 1
    print("Total number of mutations: {}".format(total_number_of_mutations))
    print("Total number of Unknown_type_count: {}".format(unknown_type_count))
    print("Total number of Unknown_Gene: {}".format(Unknown_Gene))
    print("Total number of uniq genes: {}".format(len(gene_dict)))


    for gene_id, rows in gene_dict.items():
        if len(rows) == 1:
            pass
        else:
            positions = [row["Start_Position"] for row in rows]
            if len(positions) == 2:
                if abs(int(positions[0]) - int(positions[1])) > 45:
                    print("they are not disturbed by position")
                else:
                    print("they are disturbed by position")
                    cDNA_Change1 = rows[0]["cDNA_Change"] # c.2666G>T
                    cDNA_Change2 = rows[1]["cDNA_Change"] # c.2666_2667insT
                    # print("Test: {} must be same as {}".format(rows[0]["Annotation_Transcript"], rows[1]["Annotation_Transcript"]))
                    # raw_seq = adict[rows[0]["Annotation_Transcript"]]

                    raw_seq = "GAATGCTGGGAGAGTCCGACGAGCGCTGCACTAACGCAGGATCCGGCTGCCGAAGGTCCTCGCCAGCAGGATGAAGTTAAAGGAAGTAGATCGTACAGCC"
                    print(raw_seq)
                    # test1:
                    cDNA_Change1 = "c.30_31insC"
                    cDNA_Change2 = "c.80_81insG"
                    new_seq = multi_position_fix(raw_seq, cDNA_Change1, cDNA_Change2)
                    # test2:
                    cDNA_Change1 = "c.30_31insC"
                    cDNA_Change2 = "c.50_53delCCGA" 
                    new_seq = multi_position_fix(raw_seq, cDNA_Change1, cDNA_Change2)
                    # test3:
                    cDNA_Change1 = "c.30_31insC"
                    cDNA_Change2 = "c.50C>G"
                    new_seq = multi_position_fix(raw_seq, cDNA_Change1, cDNA_Change2)
                    # test4:
                    cDNA_Change1 = "c.50delC"
                    cDNA_Change2 = "c.80_82insGC"
                    new_seq = multi_position_fix(raw_seq, cDNA_Change1, cDNA_Change2)
                    # test5:
                    cDNA_Change1 = "c.50delC"
                    cDNA_Change2 = "c.80_82delAAG"
                    new_seq = multi_position_fix(raw_seq, cDNA_Change1, cDNA_Change2)
                    # test6:
                    cDNA_Change1 = "c.50delC"
                    cDNA_Change2 = "c.80A>G"
                    new_seq = multi_position_fix(raw_seq, cDNA_Change1, cDNA_Change2)
                    # test7:
                    cDNA_Change1 = "c.50C>G"
                    cDNA_Change2 = "c.80_81insG"
                    new_seq = multi_position_fix(raw_seq, cDNA_Change1, cDNA_Change2)
                    # test8:
                    cDNA_Change1 = "c.50C>G"
                    cDNA_Change2 = "c.80_82delAAG"
                    new_seq = multi_position_fix(raw_seq, cDNA_Change1, cDNA_Change2)
                    # test9:
                    cDNA_Change1 = "c.50C>G"
                    cDNA_Change2 = "c.80A>G"
                    new_seq = multi_position_fix(raw_seq, cDNA_Change1, cDNA_Change2)
                print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
            else:
                print("Error: {}".format(gene_id))

            print("gene:{}, positions: {}".format(gene_id, positions))
         

    '''
    for index, row in maf_data.iterrows():
        if row["Variant_Classification"] in saved_Variant_Classification and row["Protein_Change"] != "NA" :
            if row["Annotation_Transcript"] in adict:
                gene_name = row["Hugo_Symbol"]
                Chromosome = row["Chromosome"]
                seq = adict[row["Annotation_Transcript"]]
                cDNA_Change = row["cDNA_Change"]
                Protein_Change = row["Protein_Change"]
                old_minigene, new_minigene = get_codon_minigene(Chromosome, seq, cDNA_Change, Protein_Change)
                output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    row["Hugo_Symbol"], row["Entrez_Gene_Id"], row["Center"], row["NCBI_Build"], row["Chromosome"],
                    row["Start_Position"], row["End_Position"], row["Strand"], row["Variant_Classification"], row["Variant_Type"], row["Reference_Allele"],
                    row["Tumor_Seq_Allele1"], row["Tumor_Seq_Allele2"], row["Tumor_Sample_Barcode"], row["Matched_Norm_Sample_Barcode"], row["dbSNP_RS"], 
                    row["Genome_Change"], 
                    row["Annotation_Transcript"], row["Transcript_Strand"], row["Transcript_Exon"], cDNA_Change,
                    row["Codon_Change"], row["Protein_Change"], row["Refseq_mRNA_Id"], row["tumor_f"], row["t_alt_count"],
                    row["t_ref_count"], row["n_alt_count"], row["n_ref_count"], row["DP"], new_minigene, old_minigene,
                    )
                )
            else:
                print("Annotation_Transcript not found {} {}".format(row["Annotation_Transcript"],row["Hugo_Symbol"]))
                print("please check the annotation file version, better use the same version as the cds file!")
    '''
    output_file.close()
    print("Finished extract minigene!")

if __name__ == '__main__':
    main()
