import sys
import re

def usage():
    print("""
    python deal_mixcr2SnS_acc_freq.py G359E3L2.mixcr.out.clonotypes.TRA.txt > G359E3L2.pair.acc_freq.txt

    see: /cygene2/work/P0035-OTS-EBV_pepmix/LMP2
    """)

clonoFile = sys.argv[1]

print("pair\tCDR3_acc\tCDR3_freq\tperfect_count\tperfect_freq")
p = re.compile(r'^(TR.*)\*00\(.*')
with open(clonoFile) as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("cloneId"):
            continue
        cloneId, cloneCount, cloneFraction, targetSequences, targetQualities, allVHitsWithScore, allDHitsWithScore, allJHitsWithScore,\
            allCHitsWithScore, allVAlignments, allDAlignments, allJAlignments, allCAlignments, nSeqFR1, minQualFR1, nSeqCDR1, minQualCDR1,\
            nSeqFR2, minQualFR2, nSeqCDR2, minQualCDR2, nSeqFR3, minQualFR3, nSeqCDR3, minQualCDR3, nSeqFR4, minQualFR4, aaSeqFR1, aaSeqCDR1,\
            aaSeqFR2, aaSeqCDR2, aaSeqFR3, aaSeqCDR3, aaSeqFR4, refPoints = line.split("\t")
        if float(cloneCount) >= 3:
            allVHitsWithScore = allVHitsWithScore.split(",")[0]
            allJHitsWithScore = allJHitsWithScore.split(",")[0]
            V = p.match(allVHitsWithScore)
            J = p.match(allJHitsWithScore)
            clonotype = V.group(1) + aaSeqCDR3 + J.group(1)
            print("{}\t{}\t{}\t{}\t{}".format(clonotype, cloneCount, cloneFraction, cloneCount, cloneFraction))


