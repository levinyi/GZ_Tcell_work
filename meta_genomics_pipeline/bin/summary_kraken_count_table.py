#!/usr/bin/python3
import os
import re
import sys
import argparse
import configparser


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('--input', action='store', dest='input_path', type=str, help = "the input folder with kraken2 log files.")
    parser.add_argument('--output', action='store', dest='output_path',type=str,  help = "the output file to write.")
    return parser.parse_args()

def main():
    parser = _argparse()
    abs_path = os.path.abspath(parser.input_path)
    output = open(parser.output_path, "w")
    output.write("sample\ttotal_reads\tconsumed_time\tclassified_reads\tpct\tunclassified_reads\tpct\n")
    pattern1 = re.compile('^(\d+).*in (\d+\.\d+\w).*')
    pattern2 = re.compile('\s+?(\d+)\ssequences\sclassified\s[(](\d+\.\d+%)[)]')
    pattern3 = re.compile('\s+?(\d+)\ssequences\sunclassified\s[(](\d+\.\d+%)[)]')
    files = [abs_path + '/' + each for each in os.listdir(abs_path) if each.endswith("log")]
    for each_file in files:
        sample_name = os.path.basename(each_file).split(".")[0]
        output_content = []
        with open(each_file, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                match1 = pattern1.match(line)
                if match1:
                    output_content.append(match1.group(1))
                    output_content.append(match1.group(2))
                    continue
                match2 = pattern2.match(line)
                if match2:
                    output_content.append(match2.group(1))
                    output_content.append(match2.group(2))
                    continue
                match3 = pattern3.match(line)
                if match3:
                    output_content.append(match3.group(1))
                    output_content.append(match3.group(2))
                    continue
        output.write("{}\t{}\n".format(sample_name, "\t".join(output_content)))
    output.close()


if __name__ == '__main__':
    main()
