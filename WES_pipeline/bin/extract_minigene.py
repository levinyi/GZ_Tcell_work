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
    20200709    Add "End_Position" field to result.
    20200702    Updated. fix cDNA_Change bugs: c.11119_11120CC>AT.
    20200628    Created.
    """.format(os.path.basename(sys.argv[0])))


def get_codon_minigene(seq, cDNA_Change, Protein_Change):
    # Four examples of cDNA_Change item:
    # c.1518_1519insAAACAGACCA 
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
    seq_aa = Seq(seq).translate()
    new_seq_aa = Seq(new_seq).translate()
    
    if aa_position < 14:
        old_minigene = seq_aa[: aa_position] + seq_aa[aa_position: aa_position + 14]
        new_minigene = new_seq_aa[: aa_position] + new_seq_aa[aa_position: aa_position + 14]
    else:
        old_minigene = seq_aa[aa_position -1 - 14: aa_position] + seq_aa[aa_position: aa_position + 14]
        new_minigene = new_seq_aa[aa_position -1 - 14: aa_position] + new_seq_aa[aa_position: aa_position + 14]

    return old_minigene, new_minigene
    # return old_minigene, new_minigene, aa_position


def main():
    if len(sys.argv) == 1:
        usage()
        sys.exit("Error: Wrong input file")

    cds_file = open(sys.argv[1], "r")
    maf_file = sys.argv[2]
    output_file = open(sys.argv[3], "w") 

    adict = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        adict[str(record.id)] = str(record.seq)
    cds_file.close()

    # read maf file as a dataframe and first line stored as header.
    maf_data = pd.read_table(maf_file, sep="\t",).fillna(value="NA")
    contain_fields = [
        "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
        "Start_Position", "End_Position","Strand", "Variant_Classification", "Variant_Type", "Reference_Allele",
        "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "Genome_Change", "Annotation_Transcript", 
        "Transcript_Strand", "Transcript_Exon", "cDNA_Change", "Codon_Change", "Protein_Change",
        "Refseq_mRNA_Id","tumor_f", "t_alt_count", "t_ref_count", "n_alt_count", 
        "n_ref_count", "DP", "Mutated_Minigene", "Wild-Type_Minigene",
    ]
    output_file.write("{}\n".format("\t".join(contain_fields)))

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
                Protein_Change = row["Protein_Change"]
                old_minigene, new_minigene = get_codon_minigene(seq, cDNA_Change, Protein_Change)
                output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    row["Hugo_Symbol"], row["Entrez_Gene_Id"], row["Center"], row["NCBI_Build"], row["Chromosome"],
                    row["Start_Position"], row["End_Position"], row["Strand"], row["Variant_Classification"], row["Variant_Type"], row["Reference_Allele"],
                    row["Tumor_Seq_Allele1"], row["Tumor_Seq_Allele2"], row["dbSNP_RS"], row["Genome_Change"], 
                    row["Annotation_Transcript"], row["Transcript_Strand"], row["Transcript_Exon"], cDNA_Change,
                    row["Codon_Change"], row["Protein_Change"], row["Refseq_mRNA_Id"], row["tumor_f"], row["t_alt_count"],
                    row["t_ref_count"], row["n_alt_count"], row["n_ref_count"], row["DP"], new_minigene, old_minigene,
                    )
                )
    output_file.close()


if __name__ == '__main__':
    main()
