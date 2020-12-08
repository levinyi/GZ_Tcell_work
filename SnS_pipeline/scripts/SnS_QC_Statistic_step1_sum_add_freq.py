import sys
import re
sum_file = sys.argv[1]
freq_file = sys.argv[2]

adict = {}
with open(freq_file, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("pair"):
            continue
        # pair,CDR3_acc, CDR3_freq, perfect_count, perfect_freq = line.split()
        # adict[pair] = perfect_freq
        a = line.split()
        adict[a[0]] = a[4]


pattern = re.compile(r'.*\.(clonotype.*)\.[ATCG]{8}.*')
with open(sum_file, "r") as f:
    print("Molecule\tCDR3aJ.Identity\tCDR3aJ.Coverage%\tCDR3aJ.Mutation\tCDR3aJ.InDel\tconA.Identity\tconA.Coverage%\tconA.Mutation\tconA.InDel\tTRAV.Identity\tTRAV.Coverage%\tTRAV.Mutation\tTRAV.InDel\tCDR3bJ.Identity\tCDR3bJ.Coverage%\tCDR3bJ.Mutation\tCDR3bJ.InDel\tconB.Identity\tconB.Coverage%\tconB.Mutation\tconB.InDel\tTRBV.Identity\tTRBV.Coverage%\tTRBV.Mutation\tTRBV.InDel\tSum.of.5\tTCR_freq")
    for line in f:
        if line.startswith(("Molecule","molecule")):
            continue
        line = line.rstrip("\n")
        c = line.split()
        sum5 = int(c[3])+int(c[4])+int(c[7])+int(c[8])+int(c[11])+int(c[12])+int(c[15])+int(c[16])+int(c[19])+int(c[20])
        match = pattern.match(c[0])
        # print(match.group(1))
        tcr_freq = adict[match.group(1)]
        # print(tcr_freq)

        print("{}\t{}\t{}".format("\t".join(c),sum5,tcr_freq))
