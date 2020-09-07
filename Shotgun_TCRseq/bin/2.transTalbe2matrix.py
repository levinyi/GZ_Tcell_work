#!/usr/bin/pythn3 
import sys,os
from collections import Counter
import numpy as np
import argparse


def usage():
    print("""Usage:
    python3 {0} <reads_barcode_file>

example:
    python3 {0} Total.G87E2.TRAB.clone.reads.barcode.txt
    python3 {0} -i Total.G208E1.TRAB.clone.reads.barcode.txt --cloneCountThreshold 19 --overallThreshold 5 --column_median 10 --row_median 10

    
Updates:

    20200818: redesign the structure. output for 5 files. and add 4 kinds of filtering parmarters.
    20200607: optimized code. replace output file name "Total.TRAB.96wells.boole.matrix.csv" to "Total.TRAB.wells.boole.matrix.csv".
    20200315: output read counts to a csv file.
    20200304: fix format from .txt to .csv to save more memory.
    20200301: created script.
    """.format(os.path.basename(sys.argv[0])))


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-i', '--input', action='store', dest='input_file', help="input file")
    parser.add_argument('-a', '--cloneCountThreshold', action='store',dest='cloneCountThreshold', default=3, type=int, help='mixcr filtering clonotypes.')
    parser.add_argument('-b', '--overallThreshold', action='store',dest='overall_threshold',default=5,type=int,help='number of filtering clonotype by overall.')
    parser.add_argument('-c', '--column_median', action='store',dest='column_median',default=10, type=int,help='10 means 10%')
    parser.add_argument('-d', '--row_median', action='store',dest='row_median',default=10, type=int, help='10 means 10%')
    parser.add_argument('-f', '--filter_noise_wells', action='store_true',dest='filter_noise_wells', default=False, help='filter out all 1 rows.')
    return parser.parse_args()


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


def addtwodimdict(thedict, key_a, key_b, value):
    if key_a in thedict:
        thedict[key_a].update({key_b: value})
    else:
        thedict.update({key_a:{key_b:value}})
    return thedict


def deal_table(input_file):
    """docstring for deal_table"""
    clonotype_dict = {}
    barcode_list = []  # for check barcode number!
    cloneid_list = []  # should equal to the input file column 1 unique
    read_count_dict = {}
    tra_barcode_dict = {}
    trb_barcode_dict = {}
    with open(input_file, "r") as f:
        for line in f:
            cloneid, readid, barcode = line.split()
            well = barcode.rstrip("\n")[4:]
            two_dim_dict(read_count_dict, cloneid, well, 1)

            if cloneid.startswith("a"):
                two_dim_dict(tra_barcode_dict, well, cloneid, 1)
            elif cloneid.startswith("b"):
                two_dim_dict(trb_barcode_dict, well, cloneid, 1)
            else:
                sys.exit("Error: clone id error, check your input file!")

            if well not in barcode_list:
                barcode_list.append(well)
            if cloneid not in cloneid_list:
                cloneid_list.append(cloneid)
            if cloneid in clonotype_dict:
                if well not in clonotype_dict[cloneid]:
                    clonotype_dict.setdefault(cloneid, []).append(well)
            else:
                clonotype_dict.setdefault(cloneid, []).append(well)
    return clonotype_dict, barcode_list, cloneid_list, read_count_dict, tra_barcode_dict, trb_barcode_dict


def main():
    """docstring for main"""
    parser = _argparse()

    input_file = parser.input_file
    column_median = parser.column_median
    row_median = parser.row_median
    cloneCountThreshold = parser.cloneCountThreshold
    overall_threshold = parser.overall_threshold
    filterNoise = parser.filter_noise_wells

    #############################################################
    # deal with input file, generate dicts.
    clonotype_dict, barcode_list, cloneid_list, read_count_dict, tra_barcode_dict, trb_barcode_dict = deal_table(input_file)
    barcode_number = len(barcode_list)
    print("Barcode numbers: {}".format(len(barcode_list)))

    prefix = input_file.split(".")[1]
    outfile1 = "Total." + prefix + ".TRAB.wells.boole.matrix.Raw.csv"
    outfile2 = "Total." + prefix + ".TRAB.wells.count.matrix.Raw.csv"
    outfile3 = "Total." + prefix + ".TRAB.wells.boole.matrix.Filtered.csv"
    outfile4 = "Total." + prefix + ".TRAB.wells.count.matrix.Filtered.csv"
    outfile5 = "Total." + prefix + ".TRAB.clone.reads.barcode.Filtered.txt"

    ###############################################################
    # create three dict for calculate row median and column median and cloneCoun.
    row_median_dict = {}
    clone_count_dict = {}
    for cloneid in cloneid_list:
        count_list = [int(value) for value in read_count_dict[cloneid].values()]
        cloneCount = np.sum(count_list)
        cloneid_in_row_median = np.median(count_list) / row_median
        row_median_dict[cloneid] = cloneid_in_row_median
        clone_count_dict[cloneid] = cloneCount
    # print(clone_count_dict)
    column_median_dict = {}
    for barcode in barcode_list:
        alist = [int(value) for value in tra_barcode_dict[barcode].values()]
        blist = [int(value) for value in trb_barcode_dict[barcode].values()]
        tra_median = np.median(alist) / column_median
        trb_median = np.median(blist) / column_median
        addtwodimdict(column_median_dict, barcode, 'TRA', tra_median)
        addtwodimdict(column_median_dict, barcode, 'TRB', trb_median)

    ################################################################     
    out1 = open(outfile1, "w")
    out2 = open(outfile2, "w")
    out3 = open(outfile3, "w")
    out4 = open(outfile4, "w")
    out5 = open(outfile5, "w")
    # write header to each file.
    out1.write("clone_id,clone_type,{}\n".format(",".join(barcode_list)))
    out2.write("clone_id,clone_type,{}\n".format(",".join(barcode_list)))
    out3.write("clone_id,clone_type,{}\n".format(",".join(barcode_list)))
    out4.write("clone_id,clone_type,{}\n".format(",".join(barcode_list)))


    ###############################################################
    # print out1, out2, these two files are raw data. 
    rid_clone_list = []
    last_boole_filtered_dict = {}
    last_count_filtered_dict ={}

    for clone in cloneid_list:
        if clone.startswith('a'):
            clone_type = "TRA"
        elif clone.startswith('b'):
            clone_type = "TRB"
        
        out1.write("{},{}".format(clone, clone_type))
        out2.write("{},{}".format(clone, clone_type))

        if clone_count_dict[clone] >= cloneCountThreshold:
            for barcode in barcode_list:
                if barcode in read_count_dict[clone]:
                    read_count = read_count_dict[clone][barcode]
                    out1.write(",1")
                    out2.write(",{}".format(read_count))
                    if read_count >= row_median_dict[clone] and read_count >= column_median_dict[barcode][clone_type] and read_count >= overall_threshold:
                        addtwodimdict(last_boole_filtered_dict, clone, barcode, 1)
                        addtwodimdict(last_count_filtered_dict, clone, barcode, read_count)
                    else:
                        addtwodimdict(last_boole_filtered_dict, clone, barcode, 0)
                        addtwodimdict(last_count_filtered_dict, clone, barcode, 0)
                        rid_clone_list.append("{}_{}".format(clone, barcode))
                else:
                    out1.write(",0")
                    out2.write(",0")
                    addtwodimdict(last_boole_filtered_dict, clone, barcode, 0)
                    addtwodimdict(last_count_filtered_dict, clone, barcode, 0)
                # out1.write("\n")
                # out2.write("\n")
        else:
            for barcode in barcode_list:
                if barcode in read_count_dict[clone]:
                    read_count = read_count_dict[clone][barcode]
                    out1.write(",1")
                    out2.write(",{}".format(read_count))
                else:
                    out1.write(",0")
                    out2.write(",0")
                rid_clone_list.append("{}_{}".format(clone, barcode))
        out1.write("\n")
        out2.write("\n")

    ##############################################################################
    # print out3, out4 : These two files are filtered data: filter out 1) all the 0 rows and 2) all the 1 rows(if --filter_noise_wells selected).
    if filterNoise:
        outfile6 = prefix + ".Clones.in.All.Wells.txt"
        out6 = open(outfile6, "w")

    for clone in cloneid_list:
        if clone in last_boole_filtered_dict:
            if np.sum([value for value in last_boole_filtered_dict[clone].values()]) == 0:
                continue
            if filterNoise:
                if np.sum([value for value in last_boole_filtered_dict[clone].values()]) == barcode_number:
                    for barcode in last_boole_filtered_dict[clone].keys():
                        rid_clone_list.append("{}_{}".format(clone, barcode))
                    # print("Filtered: {} 1 in all rows.".format(clone))
                    out6.write("{}\n".format(clone))
                    continue
            if clone.startswith('a'):
                clone_type = "TRA"
            elif clone.startswith("b"):
                clone_type = "TRB"
            out3.write("{},{}".format(clone, clone_type))
            out4.write("{},{}".format(clone, clone_type))
            for barcode in barcode_list:
                out3.write(",{}".format(last_boole_filtered_dict[clone][barcode]))
                out4.write(",{}".format(last_count_filtered_dict[clone][barcode]))
            out3.write("\n")
            out4.write("\n")

    ###############################################################################
    # print filtered file.
    with open(input_file, "r") as f:
        for line in f:
            cloneid, readid, barcode = line.split()
            well = barcode.rstrip("\n")[4:]
            if cloneid + '_' + well in rid_clone_list:
                continue
            out5.write("{}\t{}\t{}\n".format(cloneid, readid, barcode))

    out5.close()
    out4.close()
    out3.close()
    out2.close()
    out1.close()
    print("output files: {}\n{}\n{}\n{}\n{}\n".format(outfile1, outfile2, outfile3, outfile4, outfile5))
    if filterNoise:
        out6.close()
        print("{}\n".format(outfile6))

if __name__ == '__main__':
    main()
