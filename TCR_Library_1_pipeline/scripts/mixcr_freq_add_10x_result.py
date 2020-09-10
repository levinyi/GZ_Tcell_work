import sys
import os

def deal_freq_file(freq_file):
    adict = {}
    alength = 0
    aname = []
    with open(freq_file, "r") as f2:
        for line in f2:
            if line.startswith("Clonetype"):
                for e in line.split()[1:]:
                    aname.append(e)
                continue
            line = line.rstrip("\n")
            c = line.split()
            alength = len(c) -1
            adict[c[0]] = c[1:]
    return adict, alength, aname

adict, alength, aname = deal_freq_file(sys.argv[2])
bdict, blength, bname = deal_freq_file(sys.argv[3])
#print(sys.argv[2])
#print(sys.argv[3])
#print(aname, bname)

with open(sys.argv[1], "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("clonotype_id"):
            print("{},{},{}".format(line, ",".join(aname), ",".join(bname)))
            continue
        c = line.split(",")
        tra_clonotype = c[4].split("*")[0].replace("/", "") + c[10] + c[6].split("*")[0]
        trb_clonotype = c[17].split("*")[0].replace("/", "") + c[23] + c[19].split("*")[0]
        print("{},{},{}".format(line, ",".join(adict.get(tra_clonotype, ["NULL"]*alength)), ",".join(bdict.get(trb_clonotype, ["NULL"]*blength))))

