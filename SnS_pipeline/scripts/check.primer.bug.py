import sys

def usage():
    """docstring for usage
    python $1 pairs.awk6.eq2.distance.0.xls ../database/No_used.A.341.txt ../database/No_used.B.341.txt
    """

def deal_no_used(afile):
    """docstring for deal_no_used"""
    noused_list = []
    with open(afile,"r") as f:
        for line in f:
            line =  line.rstrip("\n")
            name,seq = line.split("\t")
            name = name.rstrip(("a|b"))
            noused_list.append(name)
    return noused_list

No_used_A = deal_no_used(sys.argv[2])
No_used_B = deal_no_used(sys.argv[3])

with open(sys.argv[1],"r") as f:
    for line in f:
        line = line.rstrip("\n")
        c = line.split("\t")
        if c[1] in No_used_A:
            flag = "YES"
        else:
            flag = "NO"
        if c[2] in No_used_B:
            flag2 = "YES"
        else:
            flag2 = "NO"
        print("{}\t{}\t{}".format("\t".join(c),flag,flag2))

