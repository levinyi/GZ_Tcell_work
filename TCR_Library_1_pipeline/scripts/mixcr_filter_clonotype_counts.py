import sys
import re


def usage():
    """docstring for usage
    python mixcr2frequency.noUMI.py GXXX.mixcr.out.clonotypes.TRA.txt > xxx.raw.freq.xls

    20200107: created
    """

def deal_mixcr_file(mixcr_file):
    """docstring for deal_mixcr_file"""
    print("Clonotype\tTRV\tCDR3\tTRJ\tReadsNumber")
    with open(mixcr_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
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
                if float(cloneCount) >= 3:
                    print("{}\t{}\t{}\t{}\t{}".format(clonotype, V.group(1), aaSeqCDR3, J.group(1), cloneCount))


def main():
    """docstring for main"""
    mixcr_file = sys.argv[1]
    deal_mixcr_file(mixcr_file)
    
if __name__ == '__main__':
    main()
