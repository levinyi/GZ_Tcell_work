import sys, os
import gzip
from Bio import SeqIO


def usage():
    print("""
    usage:
        python {0} <barcode_file> <reads_file> <result_file_prefix>  
        
    example:
        python {0} /cygene/work/00.test/pipeline/Shotgun_TCRseq/database/Shotgun_TCRseq.TRB.Barcodes.txt\
                reads_file 
q
    update:
    20200723: change output format, add a new column for Barcode seq.
    20200604: updated.
    20200302: created.
    """.format(os.path.basename(sys.argv[0])))


def main():
    if len(sys.argv) != 4:
        usage()
        sys.exit()

    barcode_file = sys.argv[1]
    reads_file = sys.argv[2]
    prefix_name = sys.argv[3]

    barcode_dict = {}
    barcode_dict_reverse = {}
    with open(barcode_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            bc_name, bc_seq = line.split()
            barcode_dict[bc_seq] = bc_name
            barcode_dict_reverse[bc_name] = bc_seq


    barcode2reads_dict = {}
    total_reads = 0
    if reads_file.endswith(".gz"):
        fastq_file = gzip.open(reads_file, "rt")
    else:
        fastq_file = open(reads_file, "r")

    output = open(prefix_name + ".barcode.read.txt", "w")
    for record in SeqIO.parse(fastq_file, "fastq"):
        barcode = str(record.seq[:8])
        total_reads += 1
        if barcode in barcode_dict:
            bc_name = barcode_dict[barcode]
            output.write("{}\t{}\n".format(bc_name, str(record.id)))
            if bc_name not in barcode2reads_dict:
                barcode2reads_dict[bc_name] = 1
            else:
                barcode2reads_dict[bc_name] += 1

    with open(prefix_name + ".barcode.total.rate.xls", "w") as f:
        for bc_name in barcode2reads_dict:
            count = barcode2reads_dict[bc_name]
            f.write("{}\t{}\t{}\t{}\t{}\n".format(bc_name,
                barcode_dict_reverse[bc_name],
                count,
                total_reads,
                count/float(total_reads)))
    output.close()

if __name__ == '__main__':
    main()
