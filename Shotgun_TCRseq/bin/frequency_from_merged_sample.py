import sys,os


def usage():
    print('''
    python {0} <file>
    Example:
    python {0} G113E1L2.TRB.clone.reads.barcode.txt > xxx.result.xls

    Updates:
    20200610    Created.
    '''.format(os.path.basename(sys.argv[0])))


def main():
    if len(sys.argv) != 2:
        usage()
        sys.exit("Error: wrong input file")

    total_count = 0
    clone_freq = {}
    with open(sys.argv[1], "r") as f:
        for line in f:
            line = line.rstrip("\n")
            total_count += 1
            a = line.split()
            if len(a) == 3:
                clone, read_id, Well = line.split()
            else:
                clone, read_id = line.split()

            if clone not in clone_freq:
                clone_freq[clone] = 1
            else:
                clone_freq[clone] += 1
    
    for each in clone_freq:
        counts = clone_freq[each]
        print("{}\t{}\t{}".format(each, counts, counts/float(total_count)))


if __name__ == '__main__':
    main()
