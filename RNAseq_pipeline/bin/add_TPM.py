import sys
import os
import pandas as pd


def usage():
    print("""
        python {} <xls> <tpm file>
Update:
    20200703    Created.
    """.format(os.path.basename(sys.argv[0])))


def deal_tpm_file(tpm_file):
    adict = {}
    with open(tpm_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            transcript, tpm = line.split()
            transcript = transcript.split(".")[0]
            adict[transcript] = tpm
    return adict


def main():
    if len(sys.argv) == 1:
        usage()
        sys.exit("Error: Two input files are needed!")

    afile = sys.argv[1]
    tpm_file = sys.argv[2]
    output_file = open(os.path.basename(afile).rstrip("xls") + 'add.TPM.xls', "w")

    TPM_dict = deal_tpm_file(tpm_file)
    contain_fields = [
        "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
        "Start_Position", "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele",
        "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "Genome_Change", "Annotation_Transcript",
        "Transcript_Strand", "Transcript_Exon", "cDNA_Change", "Codon_Change", "Protein_Change",
        "Refseq_mRNA_Id","tumor_f", "t_alt_count", "t_ref_count", "n_alt_count", 
        "n_ref_count", "DP", "Mutated_Minigene", "Wild-Type_Minigene", "TPM",
    ]
    output_file.write("{}\n".format("\t".join(contain_fields)))
    data = pd.read_table(afile, sep="\t")
    if sys.argv[3] == "gene_id":
        for index, row in data.iterrows():
            output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                row["Hugo_Symbol"], row["Entrez_Gene_Id"], row["Center"], row["NCBI_Build"], row["Chromosome"],
                row["Start_Position"], row["Strand"], row["Variant_Classification"], row["Variant_Type"], row["Reference_Allele"],
                row["Tumor_Seq_Allele1"], row["Tumor_Seq_Allele2"], row["dbSNP_RS"], row["Genome_Change"], row["Annotation_Transcript"],
                row["Transcript_Strand"], row["Transcript_Exon"], row["cDNA_Change"], row["Codon_Change"], row["Protein_Change"],
                row["Refseq_mRNA_Id"],row["tumor_f"], row["t_alt_count"], row["t_ref_count"], row["n_alt_count"], 
                row["n_ref_count"], row["DP"], row["Mutated_Minigene"], row["Wild-Type_Minigene"],
                TPM_dict.get(row['HGNC_Ensembl_Gene_ID'].split(".")[0], 0),
                ))
    elif sys.argv[3] == "transcript_id":
        for index, row in data.iterrows():
            # print(type(row))
            output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                row["Hugo_Symbol"], row["Entrez_Gene_Id"], row["Center"], row["NCBI_Build"], row["Chromosome"],
                row["Start_Position"], row["Strand"], row["Variant_Classification"], row["Variant_Type"], row["Reference_Allele"],
                row["Tumor_Seq_Allele1"], row["Tumor_Seq_Allele2"], row["dbSNP_RS"], row["Genome_Change"], row["Annotation_Transcript"],
                row["Transcript_Strand"], row["Transcript_Exon"], row["cDNA_Change"], row["Codon_Change"], row["Protein_Change"],
                row["Refseq_mRNA_Id"], row["tumor_f"], row["t_alt_count"], row["t_ref_count"], row["n_alt_count"], 
                row["n_ref_count"], row["DP"], row["Mutated_Minigene"], row["Wild-Type_Minigene"],
                TPM_dict.get(row['Annotation_Transcript'].split(".")[0], 0)
                ))


if __name__ == '__main__':
    main()
