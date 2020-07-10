import sys
import os
import re


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

RNAseq_mpileup_file = sys.argv[2]
RNAseq_tpm_file = sys.argv[3]
wes_mpileup_file = sys.argv[4]
mpileup_dict = deal_mpileup(RNAseq_mpileup_file)
wes_mpileup_dict = deal_mpileup(wes_mpileup_file)
tpm_dict = deal_RNA_tpm(RNAseq_tpm_file)

with open(sys.argv[1], "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("Chr"):
            print("{}\tWild-typeReads(WES)\tMutatedReads(WES)\tWild-typeReads(RNA)\tMutatedReads(RNA)\tTPM\tVaf".format(line))
            continue
        items = line.split("\t")
        Position_start = items[1]

        if items[5] == "exonic" and not items[8].startswith("synonymous") and not items[8].startswith("unknown"):
            Transcript_Ref = items[9].split(":")[1]
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
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    line,
                    wes_wild_reads,
                    wes_muta_reads,
                    mpileup_dict[Position_start].get("W", 0),
                    mpileup_dict[Position_start].get('M', 0),
                    tpm_dict.get(Transcript_Ref,"NULL"),
                    Vaf,
                    )
                )
            else:
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    line,
                    wes_wild_reads,
                    wes_muta_reads,
                    "-",
                    "-",
                    tpm_dict.get(Transcript_Ref,"NULL"),
                    Vaf,
                    )
                )
