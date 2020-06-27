import sys
import os

def usage():
    print('''
Example:
    python {0} 
Updates:
    20200610    Created.
    '''.format(os.path.basename(sys.argv[0])))


def deal_freq_file(afile):
    adict = {}
    with open(afile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            clone, count, freq = line.split()
            adict[clone] = freq
    return adict

def main():
    print(len(sys.argv))
    if len(sys.argv) == 1:
        usage()
        sys.exit("Error: Wrong input file.")
    freq1_file = sys.argv[2]
    freq2_file = sys.argv[3]
    freq3_file = sys.argv[4]
    freq4_file = sys.argv[5]
    freq1_dict = deal_freq_file(freq1_file)
    freq2_dict = deal_freq_file(freq2_file)
    freq3_dict = deal_freq_file(freq3_file)
    freq4_dict = deal_freq_file(freq4_file)

    with open(sys.argv[1], "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("Note"):
                print(line)
                continue
            tra_clone = line.split("\t")[1]
            trb_clone = line.split("\t")[4]
            print("{}\t{}\t{}\t{}\t{}".format(
                line,
                freq1_dict.get(tra_clone,0),freq2_dict.get(trb_clone,0),
                freq3_dict.get(tra_clone,0),freq4_dict.get(trb_clone,0)))

if __name__ == '__main__':
    main()
