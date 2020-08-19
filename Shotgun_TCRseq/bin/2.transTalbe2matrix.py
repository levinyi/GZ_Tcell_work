#!/usr/bin/pythn3 
import sys,os
from collections import Counter


def usage():
    print("""Usage:
    python3 {0} <reads_barcode_file>

example:
    python3 {0} Total.G87E2.TRAB.clone.reads.barcode.txt
    
Updates:
    20200607: optimized code. replace output file name "Total.TRAB.96wells.boole.matrix.csv" to "Total.TRAB.wells.boole.matrix.csv".
    20200315: output read counts to a csv file.
    20200304: fix format from .txt to .csv to save more memory.
    20200301: created script.
    """.format(os.path.basename(sys.argv[0])))


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-i', '--input', action='store', dest='input_file',help="input file")
    parser.add_argument('-ff', '--filter_full_wells', action='store',dest='',help='')
    parser.add_argument('-mc', '--mixcr_clone', action='store',dest='',help='')
    parser.add_argument('-or', '--overall', action='store',dest='',help='')
    parser.add_argument('-cm', '--colum_mean', action='store',dest='',help='')
    parser.add_argument('-rm', '--row_mean', action='store',dest='',help='')


def two_dim_dict(thedict, key_a, key_b, value):
    if key_a in thedict:
        if key_b in thedict[key_a]:
            value = thedict[key_a][key_b] + value
            thedict[key_a].update({key_b: value})
        else:
            thedict[key_a].update({key_b: value})
    else:
        thedict.update({key_a: {key_b: value}})
    return thedict


def deal_table(input_file):
    """docstring for deal_table"""
    clonotype_dict = {}
    barcode_list = []  # for check barcode number!
    cloneid_list = []  # should equal to the input file column 1 unique
    read_count_dict = {}
    with open(input_file, "r") as f:
        for line in f:
            cloneid, readid, barcode = line.split()
            well = barcode.rstrip("\n")[4:]
            two_dim_dict(read_count_dict,cloneid,well,1)
            if well not in barcode_list:
                barcode_list.append(well)
            if cloneid not in cloneid_list:
                cloneid_list.append(cloneid)
            if cloneid in clonotype_dict:
                if well not in clonotype_dict[cloneid]:
                    clonotype_dict.setdefault(cloneid, []).append(well)
            else:
                clonotype_dict.setdefault(cloneid, []).append(well)
    return clonotype_dict, barcode_list, cloneid_list, read_count_dict


def main():
    """docstring for main"""
    if len(sys.argv) != 2:
        usage()
        sys.exit("Error: Put your input file here")

    input_file = sys.argv[1] # Total.G112E1.TRAB.clone.reads.barcode.txt
    clonotype_dict, barcode_list, cloneid_list, read_count_dict = deal_table(input_file)
    print("Barcode numbers: {}".format(len(barcode_list)))

    prefix = input_file.split(".")[1]
    outfile1 = "Total."+prefix+".TRAB.wells.boole.matrix.csv"
    outfile2 = "Total."+prefix+".TRAB.wells.count.matrix.csv"

    with open(outfile1, "w") as out1, open(outfile2, "w") as out2: 
        out1.write("clone_id,clone_type,{}\n".format(",".join(barcode_list)))
        out2.write("clone_id,clone_type,{}\n".format(",".join(barcode_list)))
        for clone in cloneid_list:
            if clone.startswith('a'):
                clone_type = "TRA"
            elif clone.startswith('b'):
                clone_type = "TRB"
            out1.write("{},{}".format(clone, clone_type))
            out2.write("{},{}".format(clone, clone_type))
            for barcode in barcode_list:
                if barcode in read_count_dict[clone]:
                    read_count = read_count_dict[clone][barcode]
                    out2.write(",{}".format(read_count))
                else:
                    out2.write(",0")
                if barcode in clonotype_dict[clone]:
                    out1.write(",1")
                else:
                    out1.write(",0")
            out1.write("\n")
            out2.write("\n")
    print("output files: {}, {}".format(outfile1, outfile2))

if __name__ == '__main__':
    main()
