import sys
import os
import re
import pandas as pd
from Bio import SeqIO
from Bio import Seq


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
    20200628    Created.
    """.format(os.path.basename(sys.argv[0])))


def get_codon_minigene(seq, cDNA_Change, Protein_Change):
    # c.1518_1519insAAACAGACCA 
    # c.2953_2972delCTAAATCACACTCCTGTATC 
    # c.4113delG 
    # c.3335A>G 
    cDNA_Change = cDNA_Change.lstrip("c.")
    protein_pattern1 = re.compile(r'p\.(\d+)\_(\d+)\w+')
    protein_pattern2 = re.compile(r'p\.\w+(\d+)\w+')

    if protein_pattern1.match(Protein_Change):
        result = protein_pattern1.match(Protein_Change)
    else:
        result = protein_pattern2.match(Protein_Change)
    aa_position = int(result.group(1))

    seq_aa = Seq.Seq(seq, IUPAC.unambiguous_dna).translate()
    old_minigene = seq_aa[aa_position-14:aa_position] + seq_aa[aa_position+14+1:]

    if "ins" in cDNA_Change:
        position_components, insertion_seq = cDNA_Change.split("ins") 
        start, end = position_components.split("_")
        start = int(start)
        end = int(end)
        new_seq = seq[:start] + insertion_seq + seq[end:]
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
        # 3335A>G 
        position_components, after_mutation = cDNA_Change.split(">")
        start = position_components[:-1]
        start = int(start)
        before_mutation = position_components[-1]
        new_seq = seq[:start-1] + after_mutation + seq[start:]
    else:
        print("Ops! Unexpected Condition in cDNA_Change: {}".format(cDNA_Change))
    new_seq_aa = Seq.Seq(new_seq, IUPAC.unambiguous_dna).translate()
    new_minigene = new_seq_aa[aa_position-14:aa_position] + new_seq_aa[aa_position +14 +1:]
    
    return old_minigene, new_minigene


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

    for index, row in maf_data.iterrows():
        if row["Variant_Classification"] in saved_Variant_Classification and row["Protein_Change"] != "NA" :
            if row["Annotation_Transcript"] in adict:
                seq = adict[row["Annotation_Transcript"]]
                cDNA_Change = row["cDNA_Change"]
                old_minigene, new_minigene = get_codon_minigene(seq, cDNA_Change, Protein_Change)
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    row["Hugo_Symbol"], 
                    row["Transcript_Strand"],
                    cDNA_Change,
                    row["Codon_Change"], 
                    row["Protein_Change"],
                    old_minigene,
                    # new_minigene,
                    # minigene_aa,
                    new_minigene,
                    )
                )


if __name__ == '__main__':
    main()

