import sys
import os
import re
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO


def usage():
    """
    usage: python $0 xxx.reads > xxx.reads.cloneid.xls
    example: python G13E3L8_R1.TRA.reads > G13E3L8_R1.TRA.reads.cloneId.xls
    filtered reads number >= 3
    """


def deal_fastq_file(afastq):
    fastq_dict = {}
    header = gzip.open(afastq, "r")
    try:
        for title, seq, qual in FastqGeneralIterator(header):
            fastq_dict[title.split()[0]] = seq
    finally:
        header.close()

    return fastq_dict


def deal_reads_file(readsfile):
    adict = {}
    with open(readsfile, "r") as f:
        for line in f:
            line = line.rstrip()
            readid, cloneid = line.split("\t")
            adict.setdefault(cloneid, []).append(readid)
    return adict


def two_dim_dict(thedict, key_a, key_b, value):
    if key_a in thedict:
        if key_b in thedict[key_a]:
            value = thedict[key_a][key_b] + value
            thedict[key_a].update({key_b: value})
        else:
            thedict[key_a].update({key_b: value})
    else:
        thedict.update({key_a: {key_b: value}})
    return thedict


def deal_mixcr_file(mixcr_file):
    mixcr_dict = {}
    with open(mixcr_file, "r") as f:
        for line in f:
            if line.startswith("cloneId"):
                continue
            cloneId, cloneCount, cloneFraction, targetSequences, targetQualities, allVHitsWithScore, allDHitsWithScore, allJHitsWithScore,\
            allCHitsWithScore, allVAlignments, allDAlignments, allJAlignments, allCAlignments, nSeqFR1, minQualFR1, nSeqCDR1, minQualCDR1,\
            nSeqFR2, minQualFR2, nSeqCDR2, minQualCDR2, nSeqFR3, minQualFR3, nSeqCDR3, minQualCDR3, nSeqFR4, minQualFR4, aaSeqFR1, aaSeqCDR1,\
            aaSeqFR2, aaSeqCDR2, aaSeqFR3, aaSeqCDR3, aaSeqFR4, refPoints = line.split("\t")
            allVHitsWithScore = allVHitsWithScore.split(",")[0]
            allJHitsWithScore = allJHitsWithScore.split(",")[0]
            V = re.match(r'^(TR.*)\*00\(.*', allVHitsWithScore)
            J = re.match(r'^(TR.*)\*00\(.*', allJHitsWithScore)
            if V and J:
                clonotype = V.group(1) + aaSeqCDR3 + J.group(1)
                mixcr_dict[cloneId] = clonotype
    return mixcr_dict


def main():
    reads_file = sys.argv[1]
    mixcr_file = sys.argv[2]
    raw_fastq = sys.argv[3]

    raw_fastq_dict = deal_fastq_file(raw_fastq)
    mixcr_dict = deal_mixcr_file(mixcr_file)

    raw_umi_dict = {}
    with open(reads_file, "r") as f:
        for line in f:
            line = line.rstrip()
            readid, cloneid = line.split("\t")
            umi = raw_fastq_dict[readid][:8]
            two_dim_dict(raw_umi_dict, mixcr_dict[cloneid], umi, 1)
    
    for clonotype in raw_umi_dict:
        for umi in raw_umi_dict[clonotype]:
            print('{0}\t{1}\t{2}'.format(clonotype,umi,raw_umi_dict[clonotype][umi]))

if __name__ == '__main__':
    main()
