import sys
import re
from itertools import islice


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

#file1 = 'GB001.TableS4.Mutated.Tandem.Minigenes.V2.xls'
file1 = sys.argv[1]
RNAseq_tpm_file = sys.argv[2]
RNAseq_mpileup_file = sys.argv[3]
wes_mpileup_file = sys.argv[4]

#RNAseq_tpm_file = '1MR000001VV2.gene.counts.TPM.O.txt'
#RNAseq_mpileup_file = '1MR000001VV2.mpileup_f.check.txt'
#wes_mpileup_file = 'GB001001.WES.mpileup.txt'



# if RNAseq_tpm:
tpm_dict = deal_RNA_tpm(RNAseq_tpm_file)
# if mpileup_file:
# print(tpm_dict)
mpileup_dict = deal_mpileup(RNAseq_mpileup_file)
wes_mpileup_dict = deal_mpileup(wes_mpileup_file)


zhibentitle="""OM_ID
Gene
Transcript_Ref
Classification
Alternation_type1
Alternation_type2
Chromosome
Position_start
Position_end
Genomic_DNA_change
Coding_DNA_change
AA_change
Vaf
Mean_coverage
Checked_AAchange
CheckFlag
Wild-TypeMinigene(AminoAcid)
MutatedMinigene(AminoAcid)
Wild-typeReads(WES)
MutatedReads(WES)
Wild-typeReads(RNA)
MutatedReads(RNA)
TPM(RNA)""".split("\n")

novotitle='''Chromosome
Position_start
Position_end
Genomic_DNA_Ref
Genemic_DNA_Alt
Annotation
Gene
Transcript_Ref
Exon
Coding_DNA_change
AA_change
Checked_AAchange
CheckFlag
Wild-TypeMinigene(AminoAcid)
MutatedMinigene(AminoAcid)
Wild-typeReads(WES)
MutatedReads(WES)
Wild-typeReads(RNA)
MutatedReads(RNA)
TPM(RNA)
Vaf'''.split("\n")

if len(open(file1,"r").readline()[0].split("\t")) == 15:
    print("\t".join(zhibentitle))
else:
    print("\t".join(novotitle))


with open(file1, "r") as f:
    for line in islice(f, 0, None):
        line = line.rstrip("\n")
        if len(line.split("\t")) == 15:
            Chromosome, Position_start, Position_end,Ref,Alt, annotation,\
                Gene,Transcript_Ref,exon,Coding_DNA_change,AA_change,\
                a,b,c,d = line.split("\t")
        else:
            OM_ID, Gene, Transcript_Ref, Classification, Alternation_type1,\
                Alternation_type2, Chromosome, Position_start, Position_end,\
                Genomic_DNA_change, Coding_DNA_change, AA_change, \
                Vaf, Mean_coverage,a,b,c,d = line.split("\t")
        # print(Transcript_Ref,Position_start)
        if AA_change != '-' and AA_change != "NA":
            wes_wild_reads = int(wes_mpileup_dict[Position_start].get("W", 0))
            wes_muta_reads = int(wes_mpileup_dict[Position_start].get("M", 0))
            if Position_start in mpileup_dict:
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    line, 
                    wes_wild_reads,
                    wes_muta_reads,
                    mpileup_dict[Position_start].get("W", 0),
                    mpileup_dict[Position_start].get('M', 0),
                    tpm_dict.get(Transcript_Ref,"NULL"),
                    wes_muta_reads/float(wes_wild_reads+wes_muta_reads),
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
                    wes_muta_reads/float(wes_wild_reads+wes_muta_reads),
                    )
                )
