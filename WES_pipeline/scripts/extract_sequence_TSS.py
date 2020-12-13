import sys
import os
import argparse


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-i', '--input', action='store', dest='input_file',   help="count file, usually named G260E2L1.pair.count.txt")
    parser.add_argument('-l', '--length',  action='store', dest="extract_len",   help="frequency file, usually named G260E2L1.pair.acc_freq.txt")
    parser.add_argument('-r', '--reference',action='store', dest="reference_fa", help="reference file")
    return parser.parse_args()

def main():
    parser = _argparse()
    input_file = parser.input_file
    extract_len = parse.extract_len
    reference = parser.reference_fa

    with open(input_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            name,chrome,pos = line.split()
            print()


if __name__ == '__main__':
    main()
