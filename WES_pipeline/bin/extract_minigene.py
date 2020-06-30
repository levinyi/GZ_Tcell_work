import sys
import os
import re
import pandas as pd
from Bio import SeqIO

def usage():
    print("""
Usage:
    python {0} <maf file > <transcription file> <output file>
Example:
    python {0} 
Notes:
    This script is used for extract minigene from a pep file by a position.
    maf file without annotation line (drop lines which contains #).
    Function:
        1. filter useless sites.
            Mandatory fields: Missense_Mutation, Frame_Shift_Ins, Frame_Shift_Del, Nonesense_Mutation, Splice_site, In_Frame_Del, In_Frame_Ins.
            Silent_Mutation may used for check mutation position in transcription sequence.
        2. filter useless annotation items. 
            Mandatory fields: Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Variant_Type and Tumor_Sample_Barcode.
        3. standardization output format to a xls or csv format.

Updates:
    20200628    Created.
    """.format(os.path.basename(sys.argv[0])))


def get_codon_minigene(seq, cDNA_Change, prog):
    # c.1518_1519insAAACAGACCA
    # c.2953_2972delCTAAATCACACTCCTGTATC
    # c.4113delG
    # c.3335A>G
    if "ins" in cDNA_Change:
        pass
    elif "del" in cDNA_Change:
        pass
    elif ">" in cDNA_Change:
        pass
    
        start, end = result.group(1).split("-")
        start = int(start)
        end = int(end)

        codon = seq[start-1 : end]
        minigene = seq[start-1-42: end+42]

        # for 
        change = result.group(2)
        if change.endswith("fs"):
            mutated_minigene = seq[start-1-42: start] + "fs"
        elif change.endswith("del"):
            mutated_minigene = seq[start-1-42: start] + seq[end:end+42]
        else:
            a, b = change.split(">")
            if len(b) >= 42:
                mutated_minigene = seq[start-1-42: start-1] + b[:42]
            else:
                mutated_minigene = seq[start-1-42: start-1] +b + seq[end: end + 42 + len(a) -len(b)]

    return codon, minigene, mutated_minigene


def translate(dna):
    gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    amino_acid_sequence = ""
    for start in range(0,len(dna) - 2, 3):
        stop = start + 3
        codon = dna[start:stop]
        aa = gencode.get(codon.upper(),'X') #当指定键的值不存在时，返回X
        amino_acid_sequence = amino_acid_sequence + aa
    return(amino_acid_sequence)

def main():
    if len(sys.argv) == 1:
        usage()
        sys.exit("Error: Wrong input file")

    cds_file = open(sys.argv[1], "r")
    maf_file = sys.argv[2]

    adict = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        adict[str(record.id)] = str(record.seq)

    # read maf file as a dataframe and first line stored as header.
    # print("yes")

    maf_data = pd.read_table(maf_file, sep="\t",).fillna(value="NA")
    # print(maf_data.head())
    contain_fields = [
        "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
        "Start_Position", "End_Position", "Strand", "Variant_Classification",
        "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", 
        "dbSNP_RS", "Genome_Change","Annotation_Transcript", "Transcript_Strand", 
        "Transcript_Exon", "Transcript_Position", "cDNA_Change", "Codon_Change", 
        "Protein_Change","Refseq_mRNA_Id","tumor_f", "t_alt_count", 
        "t_ref_count", "n_alt_count", "n_ref_count", "DP",
    ]
    # maf_data_filter_column = maf_data.loc[:,contain_fields].fillna(value="NA")
    # print(maf_data_filter_column.head())

    saved_Variant_Classification = {
        "Splice_Site": 1,
        "Nonsense_Mutation": 1,
        "Missense_Mutation": 1,
        "Frame_Shift_Del": 1,
        "Frame_Shift_Ins": 1,
        "In_Frame_Del": 1, 
        "In_Frame_Ins": 1,
    }


    prog = re.compile(r'c\.\((.*)\)(.*)') # regix
    for index, row in maf_data.iterrows():
        if row["Variant_Classification"] in saved_Variant_Classification and row["Protein_Change"] != "NA" :
            if row["Annotation_Transcript"] in adict:
                seq = adict[row["Annotation_Transcript"]]
                cDNA_Change = row["cDNA_Change"]
                check, minigene, mutated_minigene = get_codon_minigene(seq, cDNA_Change, prog)
                minigene_aa = translate(minigene)
                mutated_minigene_aa = translate(mutated_minigene)
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    row["Hugo_Symbol"], 
                    row["Transcript_Strand"],
                    cDNA_Change,
                    row["Codon_Change"], 
                    row["Protein_Change"],
                    check,
                    minigene,
                    mutated_minigene,
                    minigene_aa,
                    mutated_minigene_aa)
                )


if __name__ == '__main__':
    main()
