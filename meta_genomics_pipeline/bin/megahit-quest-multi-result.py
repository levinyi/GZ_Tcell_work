import sys
import os
import pandas as pd
import argparse


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('--input', action='store', dest='input_path', type=str, help = "the input folder with kraken2 log files.")
    parser.add_argument('--output', action='store', dest='output_file',type=str,  help = "the output file to write.")
    return parser.parse_args()


def main():
    parser = _argparse()
    file_dir = os.path.abspath(parser.input_path)

    all_file_dir = [each for each in os.listdir(os.path.abspath(file_dir)) if os.path.isdir(os.path.join(os.path.abspath(file_dir),each))]
    all_file_list =[each+'/'+each +'.megahit-quast-report/report.tsv' for each in all_file_dir]

    for single_file in all_file_list:
        print("Reading report file: {}".format(os.path.join(file_dir, single_file)))
        single_data_frame = pd.read_table(os.path.join(file_dir, single_file), sep="\t")
        name = single_file.split("/")[0]
        # change the column names
        single_data_frame.columns = ['Sample', name]
        if single_file == all_file_list[0]:
            all_data_frame = single_data_frame
        else:
            all_data_frame = pd.merge(all_data_frame, single_data_frame, on="Sample")

    all_data_frame = all_data_frame.T.sort_index()
    all_data_frame.to_csv(parser.output_file, index=True, header=False, sep="\t")
    print("Summary table written: {}".format(parser.output_file))


if __name__ == '__main__':
    main()
