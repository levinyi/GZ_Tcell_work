import sys
import re


def usage():
    '''
    usage: python reads2relaxed_clonotype.py <mixcr result> <.reads file> > <output file>
    example: python reads2relaxed_clonotype.py G13AB_R1.mixcr.out.clonotypes.TRA.txt G13E3L8_R1.TRA.reads  > G13E3L8_R1.TRA.relaxed.clonotype.reads.xls
    '''


def deal_reads(afile):
    readCount_dict = {}
    with open(afile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            readId, cloneId = line.split("\t")
            readCount_dict[cloneId] = readCount_dict.get(cloneId, 0) + 1
    return readCount_dict


def main():
    mixcr_file = sys.argv[1]
    reads_file = sys.argv[2]

    readCount_dict = deal_reads(reads_file)
    a_dict = {}
    with open(sys.argv[1], "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("cloneId"):
                continue
            cloneId, cloneCount, cloneFraction, targetSequences, targetQualities, allVHitsWithScore, allDHitsWithScore, allJHitsWithScore, allCHitsWithScore, allVAlignments, allDAlignments, allJAlignments, allCAlignments, nSeqFR1, minQualFR1, nSeqCDR1, minQualCDR1, nSeqFR2, minQualFR2, nSeqCDR2, minQualCDR2, nSeqFR3, minQualFR3, nSeqCDR3, minQualCDR3, nSeqFR4, minQualFR4, aaSeqFR1, aaSeqCDR1, aaSeqFR2, aaSeqCDR2, aaSeqFR3, aaSeqCDR3, aaSeqFR4, refPoints = line.split("\t")
            allVHitsWithScore = allVHitsWithScore.split(",")[0]
            allJHitsWithScore = allJHitsWithScore.split(",")[0]
            V = re.match(r'^(TR.*)\*00\(.*', allVHitsWithScore)
            J = re.match(r'^(TR.*)\*00\(.*', allJHitsWithScore)
            C = re.match(r'^(TR.*)\*00\(.*', allCHitsWithScore)
            if V and J:
                clonoty = V.group(1) + aaSeqCDR3 + J.group(1)
                if clonoty not in a_dict:
                    if cloneId in readCount_dict:
                        a_dict[clonoty] = readCount_dict[cloneId]
                else:
                    if cloneId in readCount_dict:
                        a_dict[clonoty] += readCount_dict[cloneId]

    # output
    for k in a_dict:
        if a_dict[k] >= 3:
            print("%s\t%s" % (k, a_dict[k]))


if __name__ == '__main__':
    main()
