import sys
import re


a_dict = {}
with open(sys.argv[1],"r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("cloneId"):
            continue
        cloneId,cloneCount,cloneFraction,targetSequences,targetQualities,allVHitsWithScore,allDHitsWithScore,allJHitsWithScore,allCHitsWithScore,allVAlignments,allDAlignments,allJAlignments,allCAlignments,nSeqFR1,minQualFR1,nSeqCDR1,minQualCDR1,nSeqFR2,minQualFR2,nSeqCDR2,minQualCDR2,nSeqFR3,minQualFR3,nSeqCDR3,minQualCDR3,nSeqFR4,minQualFR4,aaSeqFR1,aaSeqCDR1,aaSeqFR2,aaSeqCDR2,aaSeqFR3,aaSeqCDR3,aaSeqFR4,refPoints = line.split("\t")
 #       print c[5],c[7],c[8]
        allVHitsWithScore = allVHitsWithScore.split(",")[0]
        allJHitsWithScore = allJHitsWithScore.split(",")[0]
        V = re.match(r'^(TR.*)\*00\(.*', allVHitsWithScore)
        J = re.match(r'^(TR.*)\*00\(.*', allJHitsWithScore)
        C = re.match(r'^(TR.*)\*00\(.*', allCHitsWithScore)
        if float(cloneCount)>=3 and V and J :
            #print("%s%s%s%s"% (targetSequences,V.group(1),J.group(1),C.group(1)))
            clonoty = V.group(1) + aaSeqCDR3 + J.group(1)
            if clonoty not in a_dict:
                a_dict[clonoty] = float(cloneCount)
            else:
                a_dict[clonoty] = a_dict[clonoty] + float(cloneCount)

for k in a_dict :
    print("%s\t%s"% (k,a_dict[k]))
