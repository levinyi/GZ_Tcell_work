import sys


def usage():
    """
    usage: python $0 xxx.reads > xxx.reads.cloneid.xls
    example: python G13E3L8_R1.TRA.reads > G13E3L8_R1.TRA.reads.cloneId.xls
    filtered reads number >= 3
    """


def main():
    """docstring for main"""
    reads_file = sys.argv[1]  # G22E1L1_R1.TRA.reads

    adict = {}
    with open(reads_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            read_id, clone_id = line.split("\t")
            adict[clone_id] = adict.get(clone_id, 0) + 1

    # output
    print("cloneId\treads")
    for k in adict:
        if adict[k] >= 3:
            print("%s\t%s" % (k, adict[k]))

if __name__ == '__main__':
    main()
