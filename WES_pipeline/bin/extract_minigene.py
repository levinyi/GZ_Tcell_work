#! coding: utf-8
import sys
import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import warnings
from Bio import BiopythonWarning
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
    20220811    fixed a bug in two continuous mutations in one codon. LC24:ARID3A:c.(1159-1161)aAT>aGC:p.N387S
    20220801    fixed nonstop mutation bug.
    20220728    fixed ID_for_order too long in save_minigene function.
    20220711    Add 5 columns to output file: doNotsyn,MutMG_ID_for_order, Mut_AA29_for_order, wtMG_ID_for_order, wt_AA29_for_order.
    20220217    Fix bugs: DtypeWarning: Columns (56,88,101,103,107,108,111,147) have mixed types.Specify dtype option on import or set low_memory=False. [at pandas read file]
    20210423    Fix bugs: if mutation occurs in Mitochondrial, use coding_dna.translate(table="Vertebrate Mitochondrial")
    20200709    Add "End_Position" field to result.
    20200702    Updated. fix cDNA_Change bugs: c.11119_11120CC>AT.
    20200628    Created.
    """.format(os.path.basename(sys.argv[0])))

def single_position_fix(raw_seq, cDNA_Change):
    c1_type, c1_start, c1_end, c1_seq = get_cDNA_Change_Info(cDNA_Change)
    print("c1_type: {}, c1_start: {}, c1_end: {}, c1_seq: {}".format(c1_type, c1_start, c1_end, c1_seq))
    if c1_type == "ins":
        full_seq = raw_seq[:c1_start] + c1_seq + raw_seq[c1_start:]
    elif c1_type == "del":
        full_seq = raw_seq[:c1_start-1] + raw_seq[c1_end:]
    elif c1_type == ">":
        full_seq = raw_seq[:c1_start-1] + c1_seq + raw_seq[c1_end:]
    return full_seq


def multi_position_fix(raw_seq, cDNA_Change1, cDNA_Change2):
    # Examples for cDNA_Change:
    # ins: c.1518_1519insAAACAGACCA 
    # del: c.2953_2972delCTAAATCACACTCCTGTATC , c.4113delG 
    # >  : c.3335A>G , c.11119_11120CC>AT
    start1 = int(re.findall(r'\d+', cDNA_Change1)[0])
    start2 = int(re.findall(r'\d+', cDNA_Change2)[0])

    # sort cDNA_Change1 and cDNA_Change2 by position.
    if start1 > start2:
        cDNA_Change1, cDNA_Change2 = cDNA_Change2, cDNA_Change1
    
    # get the Type of cDNA_Change1 and cDNA_Change2.
    c1_type, c1_start, c1_end, c1_seq = get_cDNA_Change_Info(cDNA_Change1)
    c2_type, c2_start, c2_end, c2_seq = get_cDNA_Change_Info(cDNA_Change2)
    
    # combine partial sequence.
    if c1_type == "ins":
        partial_seq = raw_seq[:c1_start] + c1_seq + raw_seq[c1_start:c2_start]
    elif c1_type == "del":
        partial_seq = raw_seq[:c1_start-1] + raw_seq[c1_end: c2_start]
    elif c1_type == ">":
        partial_seq = raw_seq[:c1_start-1] + c1_seq + raw_seq[c1_start:c2_start]

    # combine full sequence.
    if c2_type == "ins":
        full_seq = partial_seq + c2_seq + raw_seq[c2_start:]
    elif c2_type == "del":
        full_seq = partial_seq[:-1] + raw_seq[c2_end:]
    elif c2_type == ">":
        full_seq = partial_seq[:c2_start-1] + c2_seq + partial_seq[c2_end:]+ raw_seq[c2_start:]
    return full_seq


def get_cDNA_Change_Info(cDNA_Change):
    # Examples for cDNA_Change:
    # ins: c.1518_1519insAAACAGACCA 
    # del: c.2953_2972delCTAAATCACACTCCTGTATC , c.4113delG 
    # >  : c.3335A>G , c.11119_11120CC>AT
    cDNA_Change = cDNA_Change.lstrip("c.")
    if "ins" in cDNA_Change:
        position_components, middle_seq = cDNA_Change.split("ins") 
        start, end = position_components.split("_")
        return "ins", int(start), int(end), middle_seq
    elif "del" in cDNA_Change:
        position_components, middle_seq = cDNA_Change.split("del")
        if "_" in position_components:
            start, end = position_components.split("_")
            start = int(start)
            end = int(end)
        else:
            start = position_components
            start = int(start)
            end = int(start) # end is the same as start.
        return "del", start, end, middle_seq
    elif ">" in cDNA_Change:
        position_components, after_mutation = cDNA_Change.split(">")
        starts = re.findall(r'\d+', position_components)
        if len(starts) == 2:
            start = int(starts[0])
            end = int(starts[1])
        elif len(starts) == 1:
            start = int(starts[0])
            end = int(starts[0]) # end is the same as start.
        else:
            print("Ops! Unexpected format in cDNA_Change: {}".format(cDNA_Change))
        return ">", start, end, after_mutation
    else:
        print("Ops! Unexpected Condition in cDNA_Change: {}. should be { ins, del, > }".format(cDNA_Change))
        sys.exit(1)

def loop_minigene(aa_position, start_position, end_position, stop_codon_position, raw_seq_aa, new_seq_aa):
    # all the parameters that from here are 0-based.
    new_minigene_list = []
    raw_minigene_list = []

    # loop to get the minigene.
    if stop_codon_position < end_position:
        a = new_seq_aa[start_position: stop_codon_position+1] 
        b = raw_seq_aa[start_position: end_position]
        new_minigene_list.append(a)
        raw_minigene_list.append(b)
    else:
        while stop_codon_position > end_position: # if stop codon position is after end_position, it means the stop codon is in the next minigene.
            a = new_seq_aa[start_position:end_position]
            if a not in new_minigene_list:
                new_minigene_list.append(a)
            if start_position + 1 < aa_position: # wt 只输出第一个minigene，因为后面循环的minigene已经无法对应wt序列了，所以只输出第一个WT minigene.
                b = raw_seq_aa[start_position:end_position]
                if b not in raw_minigene_list:
                    raw_minigene_list.append(b)
            start_position += 14  # adjust the start position to the next minigene's start position. 14表示下一个minigene的起始位置包含了前一个的mutation位置，15表示下一个minigene的起始位置不包含前一个的mutation位置。
            end_position = start_position + 29
        else:
            a = new_seq_aa[stop_codon_position-29:stop_codon_position+1] # 末尾把stop codon位置的碱基也顺带输出：如果是*就输出*，如果是碱基就输出碱基。
            if len(a.rstrip("*"))>29:
                a = a[:29] # 防止超过29个碱基，只输出前29个碱基。
            if a not in new_minigene_list:
                new_minigene_list.append(a)
    return new_minigene_list, raw_minigene_list


def extract_minigene(raw_seq, new_seq, Protein_Change, Chromosome, Protein_Change2=None):  
    if Protein_Change2: # optional
        start1 = int(re.findall(r'\d+', Protein_Change)[0])
        start2 = int(re.findall(r'\d+', Protein_Change2)[0])
        # sort Protein_Change and Protein_Change2 by position.
        if start1 > start2:
            Protein_Change, Protein_Change2 = Protein_Change2, Protein_Change
    # if Chromosome is Mitochondrial. Use Seq(seq).translate(table="Vertebrate Mitochondrial") otherwise it will raise error.
    if Chromosome == "MT":
        raw_seq_aa = str(Seq(raw_seq).translate(table="Vertebrate Mitochondrial")) +"*" # add * to the end of the raw_seq_aa. fake the stop codon.
        new_seq_aa = str(Seq(new_seq).translate(table="Vertebrate Mitochondrial")) # if no * in the end of the new_seq_aa. turn to manualy check. otherwise turn to ...
    else:
        raw_seq_aa = str(Seq(raw_seq).translate())
        new_seq_aa = str(Seq(new_seq).translate())

    aa_position = int(re.findall(r'\d+', Protein_Change)[0]) -1
    print("aa_position: {}".format(aa_position))
    print("raw_seq_aa: {}, new_seq_aa: {}".format(raw_seq_aa, new_seq_aa))
    if aa_position <= 14: # if aa_position occurs before 14aa, it means the mutation is in the first 14aa.
        start_position = 0
    else:
        start_position = aa_position - 14
    end_position = aa_position + 15 # string cut rule
    print("start_position: {}, end_position: {}".format(start_position, end_position))
    # if new_seq_aa contain "*", it means the stop codon occured. and find the position of "*" in new_seq_aa.
    # 下面的条件中所有的position都必须调整为0-based.
    if "*" in new_seq_aa[:-1] :  # contain * but not the last one.
        stop_codon_position = new_seq_aa.find("*") # find first stop codon position in new_seq_aa. 0-based.
        # print("stop_codon_position: {}".format(stop_codon_position)) # 86
        new_minigene_list, raw_minigene_list = loop_minigene(aa_position, start_position, end_position, stop_codon_position, raw_seq_aa, new_seq_aa)
    elif "*" == new_seq_aa[-1]: # last one is *.
        # print("new_seq_aa: {}, raw_seq_aa: {}".format(new_seq_aa, raw_seq_aa))
        # protein_change 中只有ins是有两个数字的，p.234_235insGESTSFFFFF. 如果插入的很长，就需要loop输出minigene.
        # 并且要注意，aaposition的位置是从后面的数字开始的,也就是要加1位。
        if "ins" in Protein_Change: 
            aa_position = aa_position +1
            if aa_position <= 14:
                start_position = 0
            else:
                start_position = aa_position - 14
            end_position = aa_position + 15
            insertion = Protein_Change.lstrip("p.").split("ins")[-1]
            stop_codon_position = aa_position + len(insertion) + 14
            # print("aa_position: {}, start_position: {}, end_position: {}, stop_codon_position: {}".format(aa_position, start_position, end_position, stop_codon_position))
            new_minigene_list, raw_minigene_list = loop_minigene(aa_position,start_position,end_position,stop_codon_position,raw_seq_aa,new_seq_aa)
        else:
            new_minigene_list = [new_seq_aa[start_position:end_position]]
            raw_minigene_list = [raw_seq_aa[start_position:end_position]]
        
        if Protein_Change2:
            aa_position = int(re.findall(r'\d+', Protein_Change2)[0]) -1 # 0-based.
            if aa_position < 14:
                start_position = 0
            else:
                start_position = aa_position - 14
            new_minigene_list.append(new_seq_aa[start_position:aa_position + 15])
            raw_minigene_list.append(raw_seq_aa[start_position:aa_position + 15])
    else: # not contain *
        stop_codon_position = raw_seq_aa.find("*")
        # 这里不再需要判断是否有Protein_Change2，因为如果有的话，也是循环输出，跟protein_change2没关系。
        new_minigene_list, raw_minigene_list = loop_minigene(aa_position, start_position, end_position, stop_codon_position, raw_seq_aa, new_seq_aa)
        new_minigene_list[-1] = new_minigene_list[-1] + "*?"
    while len(raw_minigene_list) != len(new_minigene_list):
        raw_minigene_list.append(raw_minigene_list[0])
    return raw_minigene_list, new_minigene_list


def save_minigene(output_file, raw_minigene_list, new_minigene_list, row, gene_id, cDNA_Change, Protein_Change):
    for index, new_minigene in enumerate(new_minigene_list):
        # if nonsense_mutation:
        if new_minigene.rstrip("*") in raw_minigene_list[index]:
            doNotSyn = "doNotSyn"
        else:
            doNotSyn = ""
        
        Mut_AA29_for_order = new_minigene
        WT_AA29_for_order = raw_minigene_list[index]

        if len(gene_id) >10:
            gene_id = gene_id[:10]+"*"
        
        wtMG_ID_for_order = gene_id + "-" + str(re.findall(r'\d+', Protein_Change)[0]) + "WT"
        ### for MutMG_ID for order ###
        # if the length of the string before number or after number is great than 4, neet to replace the string as "X1,X2,X3,X4"
        # p.I266L, p.N179KGLKNISNLAY*P
        # p.20_21insFRGTDYENVQIHMDGTHMDLYI, 
        # p.ANPTGPAANPPATTANPPAPAANPSAPA236del, p.KRE866del
        # p.3480_3490LLL>FMC, p.301_302LE>F*
        if len(Protein_Change.lstrip("p.")) > 12:
            a = re.findall(r'\d+', Protein_Change.lstrip("p."))[0]
            if "ins" in Protein_Change:
                MutMG_ID_for_order = gene_id + "-" + Protein_Change.lstrip("p.").split("ins")[0] + "insX1"
            elif "del" in Protein_Change:
                MutMG_ID_for_order = gene_id + "-X1" +a + "del"
            elif ">" in Protein_Change:
                MutMG_ID_for_order =gene_id + "-" + a+ "X1>X1"
            else:
                a = re.findall(r'.\d+', Protein_Change.lstrip("p."))[0]
                MutMG_ID_for_order = gene_id + "-" +a + "X1"
        else:
            MutMG_ID_for_order = gene_id + "-" + Protein_Change.lstrip("p.")
        # print("MutMG_ID_for_order: {}".format(MutMG_ID_for_order))

        if "?" in new_minigene:
            # manually check the mutation.
            doNotSyn = "ManualCheck"

        if "*" in new_minigene:
            stop_position = new_minigene.index("*")
            add_GSG_number = 29 - stop_position   # 29 is the length of the amino acid.
            add_GSG_aa = "GSG" * add_GSG_number
            Mut_AA29_for_order = new_minigene[:stop_position] + add_GSG_aa[:add_GSG_number]  # add GSG to the end of the amino acid.

        if "*" in raw_minigene_list[index]:
            stop_position = raw_minigene_list[index].index("*")
            add_GSG_number = 29 - stop_position   # 29 is the length of the amino acid.
            add_GSG_aa = "GSG" * add_GSG_number
            WT_AA29_for_order = raw_minigene_list[index][:stop_position] + add_GSG_aa[:add_GSG_number]

        # chech minigene is shorter than 29 AA. and "*" is not in the end of the amino acid.
        if len(Mut_AA29_for_order) < 29 and "*" not in Mut_AA29_for_order:
            a = 29 - len(Mut_AA29_for_order)
            Mut_AA29_for_order = ("GSG" *a)[:a][::-1] + Mut_AA29_for_order
            WT_AA29_for_order = ("GSG" * a)[:a][::-1] + WT_AA29_for_order
        
        if index+1 > 1:
            MutMG_ID_for_order += str(index+1)
            wtMG_ID_for_order += str(index+1)
            output_file.write("{}\t\t\t\t\t\t\t\tLongChangedPeptide\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                row["Hugo_Symbol"], new_minigene, raw_minigene_list[index], doNotSyn, MutMG_ID_for_order, Mut_AA29_for_order,
                wtMG_ID_for_order, WT_AA29_for_order))
        else:
            print("MutMG_ID_for_order: {}".format(MutMG_ID_for_order))
            output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                row["Hugo_Symbol"], row["Entrez_Gene_Id"], row["Center"], row["NCBI_Build"], row["Chromosome"],  # 1-5
                row["Start_Position"], row["End_Position"], row["Strand"], row["Variant_Classification"], row["Variant_Type"],  # 6-10
                row["Reference_Allele"], row["Tumor_Seq_Allele1"], row["Tumor_Seq_Allele2"], row["Tumor_Sample_Barcode"], row["Matched_Norm_Sample_Barcode"],   # 11-15
                row["dbSNP_RS"], row["Genome_Change"], row["Annotation_Transcript"], row["Transcript_Strand"], row["Transcript_Exon"],  # 16-20
                cDNA_Change, row["Codon_Change"], row["Protein_Change"], row["Refseq_mRNA_Id"], row["tumor_f"],     # 21-25
                row["t_alt_count"], row["t_ref_count"], row["n_alt_count"], row["n_ref_count"], row["DP"],  # 26-30
                new_minigene, raw_minigene_list[index], doNotSyn, MutMG_ID_for_order, Mut_AA29_for_order,  # 31-35
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
    contain_fields = [ # from MAF table.
        "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", # 1-5
        "Start_Position", "End_Position","Strand", "Variant_Classification", "Variant_Type",  # 6-10
        "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",  # 11-15
        "dbSNP_RS", "Genome_Change", "Annotation_Transcript", "Transcript_Strand", "Transcript_Exon",  # 16-20
        "cDNA_Change", "Codon_Change", "Protein_Change", "Refseq_mRNA_Id","tumor_f", # 21-25
        "t_alt_count", "t_ref_count", "n_alt_count", "n_ref_count", "DP",  # 26-30
        "Mutated_Minigene", "Wild-Type_Minigene", "doNotSyn","MutMG_ID_for_order","Mut_AA29_for_order", # 31-35
        "wtMG_ID_for_order","WT_AA29_for_order", # 36-37
    ]
    output_file.write("{}\n".format("\t".join(contain_fields))) # write header.

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
        if gene_id == "Unknown":
            Unknown_Gene += 1
        elif var_type in saved_Variant_Classification and row['Protein_Change'] != "NA":
            useful_mutation += 1
            gene_dict.setdefault(gene_id,[]).append(row)
        else:
            unknown_type_count += 1
    print("Total number of useful mutations: {}".format(useful_mutation))
    print("Total number of Unknown_type_count: {}".format(unknown_type_count))
    print("Total number of Unknown_Gene: {}".format(Unknown_Gene))
    print("Total number of uniq genes: {}".format(len(gene_dict)))

    # start to extract minigene for each gene(line in the MAF file).
    for gene_id, rows in gene_dict.items():
        cDNA_Change = rows[0]["cDNA_Change"]
        Protein_Change = rows[0]["Protein_Change"]
        raw_seq = adict[rows[0]["Annotation_Transcript"].split(".")[0]] # warning: may be bug here. need to check.
        Chromosome = rows[0]["Chromosome"]
        # print(rows[0]['Hugo_Symbol'], cDNA_Change, rows[0]['Annotation_Transcript'], Protein_Change, Chromosome, rows[0]['Variant_Classification'])
        if len(rows) == 1: # only one mutation in this gene
            new_seq = single_position_fix(raw_seq, cDNA_Change)
            print("gene_id:{}, raw_seq:{}, new_seq:{}".format(gene_id, raw_seq, new_seq))
            raw_minigene_list, new_minigene_list = extract_minigene(raw_seq, new_seq, Protein_Change, Chromosome)
            print("Protein_Change: {}, raw_minigene_list: {}, new_minigene_list: {}".format(Protein_Change, raw_minigene_list, new_minigene_list))
            save_minigene(output_file, raw_minigene_list, new_minigene_list, rows[0], gene_id, cDNA_Change, Protein_Change)
        elif len(rows) == 2: # two mutations in this gene
            positions = [row["Start_Position"] for row in rows] # get all the positions of this gene
            # print("Detect Multiple({}) mutations in this gene: {}, they are in {}".format(len(positions),gene_id, positions))
            if abs(int(positions[0]) - int(positions[1])) > 45: # if the distance between two mutations is larger than 45, they are independent event.
                # deal with the first mutation
                new_seq = single_position_fix(raw_seq, cDNA_Change)
                raw_minigene_list, new_minigene_list = extract_minigene(raw_seq, new_seq, Protein_Change, Chromosome)
                save_minigene(output_file, raw_minigene_list, new_minigene_list, rows[0], gene_id, cDNA_Change, Protein_Change)
                # deal with the second mutation
                cDNA_Change2 = rows[1]["cDNA_Change"]
                Protein_Change2 = rows[1]["Protein_Change"]
                new_seq = single_position_fix(raw_seq, cDNA_Change2)
                raw_minigene_list, new_minigene_list = extract_minigene(raw_seq, new_seq, Protein_Change2, Chromosome)
                save_minigene(output_file, raw_minigene_list, new_minigene_list, rows[0], gene_id, cDNA_Change2, Protein_Change2)
            else: # if the distance between two mutations is less than 45bp, they are disturbed by position. they are simultaneous event.
                cDNA_Change2 = rows[1]["cDNA_Change"] # c.2666_2667insT
                # print("they are disturbed by position {} {} {} {}".format(positions[0], positions[1], cDNA_Change, cDNA_Change2))
                Protein_Change2 = rows[1]["Protein_Change"]
                new_seq = multi_position_fix(raw_seq, cDNA_Change, cDNA_Change2)
                print("gene_id:{}, raw_seq:{}, new_seq:{}".format(gene_id, raw_seq, new_seq))
                raw_minigene_list, new_minigene_list = extract_minigene(raw_seq, new_seq, Protein_Change, Chromosome, Protein_Change2)
                save_minigene(output_file, raw_minigene_list, new_minigene_list, rows[0], gene_id, cDNA_Change, Protein_Change)
                ######## save the second mutation record ##############
                row = rows[1]
                output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tMerged\n".format(
                    row["Hugo_Symbol"], row["Entrez_Gene_Id"], row["Center"], row["NCBI_Build"], row["Chromosome"],  # 1-5
                    row["Start_Position"], row["End_Position"], row["Strand"], row["Variant_Classification"], row["Variant_Type"],  # 6-10
                    row["Reference_Allele"], row["Tumor_Seq_Allele1"], row["Tumor_Seq_Allele2"], row["Tumor_Sample_Barcode"], row["Matched_Norm_Sample_Barcode"],   # 11-15
                    row["dbSNP_RS"], row["Genome_Change"], row["Annotation_Transcript"], row["Transcript_Strand"], row["Transcript_Exon"],  # 16-20
                    cDNA_Change2, row["Codon_Change"], row["Protein_Change"], row["Refseq_mRNA_Id"], row["tumor_f"],     # 21-25
                    row["t_alt_count"], row["t_ref_count"], row["n_alt_count"], row["n_ref_count"], row["DP"],  # 26-30
                    ) 
                )
        else:
            print("more than 3 positions in this gene: {}. please check manually.".format(gene_id))
            # print rows information
            for row in rows:
                output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tMulti-mutations\t\tManualCheck\n".format(
                    row["Hugo_Symbol"], row["Entrez_Gene_Id"], row["Center"], row["NCBI_Build"], row["Chromosome"],  # 1-5
                    row["Start_Position"], row["End_Position"], row["Strand"], row["Variant_Classification"], row["Variant_Type"],  # 6-10
                    row["Reference_Allele"], row["Tumor_Seq_Allele1"], row["Tumor_Seq_Allele2"], row["Tumor_Sample_Barcode"], row["Matched_Norm_Sample_Barcode"],   # 11-15
                    row["dbSNP_RS"], row["Genome_Change"], row["Annotation_Transcript"], row["Transcript_Strand"], row["Transcript_Exon"],  # 16-20
                    row['cDNA_Change'], row["Codon_Change"], row["Protein_Change"], row["Refseq_mRNA_Id"], row["tumor_f"],     # 21-25
                    row["t_alt_count"], row["t_ref_count"], row["n_alt_count"], row["n_ref_count"], row["DP"],  # 26-30
                    )
                )
    output_file.close()
    print("Finished extract minigene!")

if __name__ == '__main__':
    main()
