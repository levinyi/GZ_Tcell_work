#!/usr/bin/python3
import sys
from Bio import SeqIO
import gzip
sys.path.append("/cygene/script/python_packages")
import edit_distance


def main():


    rd1_fastq_files = sys.argv[1:5]
    
    # barcode_dict = {}
    barcode_list_dict = {}
    rd1_dict = {}
    for each_file in rd1_fastq_files:
        fastq_file = gzip.open(each_file, "rt")  # in python3, open mode.
        for record in SeqIO.parse(fastq_file, "fastq"):
            barcode = str(record.seq)[:16]
            rd1_dict[str(record.id)] = record

            barcode_list_dict.setdefault(barcode, []).append(str(record.id))

            # if barcode not in barcode_dict:
            #   barcode_dict[barcode] = 1
            # else:
            #   barcode_dict[barcode] += 1
        fastq_file.close()

    rd2_fastq_files = sys.argv[5:]
    rd2_dict = {}
    for each_file in rd2_fastq_files:
        fastq_file = gzip.open(each_file, "rt")
        for record in SeqIO.parse(fastq_file, "fastq"):
            rd2_dict[str(record.id)] = record
        fastq_file.close()

    for barcode in barcode_list_dict:
        if len(barcode_list_dict[barcode]) >5000:
            with open(barcode + "_S1_L001_R1_001.fastq", "w") as f1, open(barcode + "_S1_L001_R2_001.fastq", "w") as f2:
                for read_id in barcode_list_dict[barcode]:
                    SeqIO.write(rd1_dict[read_id], f1, "fastq")
                    SeqIO.write(rd2_dict[read_id], f2, "fastq")


if __name__ == '__main__':
    main()
