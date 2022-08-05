import sys
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict

def deal_RNA_tpm(RNAseq_tpm_file):
    a_dict = {}
    with open(RNAseq_tpm_file,"r") as f:
        for line in f:
            line = line.rstrip("\n")
            a,b = line.split("\t")
            a_dict[a.split(".")[0]] = b
    return a_dict


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


def deal_cds_file(cds_file):
    adict = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        adict[str(record.id)] = str(record.seq)
    return adict


def extract_minigene(cds_dict, Transcript_Ref, cDNA_Change, Protein_Change):
    if Transcript_Ref in cds_dict:
        seq = cds_dict[Transcript_Ref]
    else:
        return "transcript_ref_not_found", "transcript_ref_not_found"
    #c.10296_10297insTTTTTTAATGATA
    #c.1003_1004delinsGT
    #c.1003_1005del
    #c.1003delG
    #c.1005dupT
    #c.T521C

    cDNA_Change = cDNA_Change.lstrip("c.")
    if "delins" in cDNA_Change:
        position_components, insertion_seq = cDNA_Change.split("delins")
        # c.1030_1031delinsAA
        start, end = position_components.split("_")
        start = int(start)
        end = int(end)
        new_seq = seq[:start] + insertion_seq + seq[end-1:]
    elif "ins" in cDNA_Change:
        position_components, insertion_seq = cDNA_Change.split("ins")
        # c.204_205insT
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
    elif "dup" in cDNA_Change:
        # c.1065dupA
        starts, dup_seq = cDNA_Change.split("dup")
        start = int(starts)
        new_seq = seq[:start-1] + dup_seq + seq[start:]
    elif cDNA_Change.startswith(("A","T","C","G")):
        pos = re.findall(r"\d+",cDNA_Change)[0]
        subtitution = re.findall(r"[ATCG]", cDNA_Change)[-1]
        start = int(pos)
        new_seq = seq[:start-1] + subtitution + seq[start+1:]
    else:
        print("Ops! Unexpected Condition in cDNA_Change: {}".format(cDNA_Change))

    # Four examples of Protein_Change item:
    # p.R41Sfs*2
    # p.X442delinsX
    # p.S231L
    Protein_Change = Protein_Change.split(",")[0]  # p.Q74R,Tmem70       p.X442delinsX,Tfap2b

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

RNAseq_mpileup_file = sys.argv[2]
RNAseq_tpm_file = sys.argv[3]
wes_mpileup_file = sys.argv[4]
cds_file = sys.argv[5]

mpileup_dict = deal_mpileup(RNAseq_mpileup_file)
wes_mpileup_dict = deal_mpileup(wes_mpileup_file)
tpm_dict = deal_RNA_tpm(RNAseq_tpm_file)
cds_dict = deal_cds_file(cds_file)

with open(sys.argv[1], "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("Chr"):
            print("{}\tWild-TypeMinigene\tMutatedMinigene\tWild-typeReads(WES)\tMutatedReads(WES)\tWild-typeReads(RNA)\tMutatedReads(RNA)\tTPM\tVaf".format(line))
            continue
        items = line.split("\t")
        Position_start = items[1]

        if items[5] == "exonic" and not items[8].startswith("synonymous") and not items[8].startswith("unknown"):
            Transcript_Ref = items[9].split(":")[1]
            cDNA_Change = items[9].split(":")[3]
            Protein_Change = items[9].split(":")[4]
            # print(Transcript_Ref, cDNA_Change, Protein_Change)
            wildtype_minigene, mutated_minigene = extract_minigene(cds_dict, Transcript_Ref,cDNA_Change, Protein_Change)
            # print(wildtype_minigene,mutated_minigene)
            if Position_start in wes_mpileup_dict:
                wes_wild_reads = int(wes_mpileup_dict[Position_start].get("W", 0))
                wes_muta_reads = int(wes_mpileup_dict[Position_start].get("M", 0))
            else:
                wes_wild_reads = 0
                wes_muta_reads = 0
            try:
                Vaf = wes_muta_reads/float(wes_wild_reads+wes_muta_reads)
            except ZeroDivisionError:
                Vaf = 0

            if Position_start in mpileup_dict:
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    line,
                    wildtype_minigene,
                    mutated_minigene,
                    wes_wild_reads,
                    wes_muta_reads,
                    mpileup_dict[Position_start].get("W", 0),
                    mpileup_dict[Position_start].get('M', 0),
                    tpm_dict.get(Transcript_Ref,"NULL"),
                    Vaf,
                    )
                )
            else:
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    line,
                    wildtype_minigene,
                    mutated_minigene,
                    wes_wild_reads,
                    wes_muta_reads,
                    "-",
                    "-",
                    tpm_dict.get(Transcript_Ref,"NULL"),
                    Vaf,
                    )
                )

