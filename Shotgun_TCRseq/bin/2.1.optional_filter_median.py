#!/usr/bin/pythn3
import sys
import os 
import numpy as np
from collections import Counter

def usage():
    print("""Usage:
    python3 {0} <reads_barcode_file>

example:
    python3 {0} Total.G87E2.TRAB.clone.reads.barcode.txtpython3 /cygene/work/00.test/pipeline/Shotgun_TCRseq/bin/2.1.optional_filter_median.py Total.G125E2.TRAB.clone.reads.barcode.txt

Updates:
    20200724: change median.
    20200723: created script.
    """.format(os.path.basename(sys.argv[0])))


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
            if cloneid.startswith("a"):
                two_dim_dict(tra_barcode_dict, well, cloneid, 1)
            elif cloneid.startswith("b"):
                two_dim_dict(trb_barcode_dict, well, cloneid, 1)
            two_dim_dict(read_count_dict, cloneid, well, 1)
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
    if len(sys.argv) != 2:
        usage()
        sys.exit("Error: Put your input file here")

    input_file = sys.argv[1] # Total.G112E1.TRAB.clone.reads.barcode.txt
    clonotype_dict, barcode_list, cloneid_list, read_count_dict, tra_barcode_dict, trb_barcode_dict = deal_table(input_file)
    barcode_number = len(barcode_list)
    print("Barcode numbers: {}".format(barcode_number))


    prefix = input_file.split(".")[1]
    outfile1 = "Total."+prefix+".TRAB.wells.boole.matrix.Filtered.csv"
    outfile2 = "Total."+prefix+".TRAB.wells.count.matrix.Filtered.csv"
    outfile3 = "Total."+prefix+".TRAB.clone.reads.barcode.Filtered.txt"
    

    barcode_median_dict = {}
    for barcode in barcode_list:
        alist = [int(value) for value in tra_barcode_dict[barcode].values()]
        blist = [int(value) for value in trb_barcode_dict[barcode].values()]
        # print("alist:{}".format(alist))
        # print("blist:{}".format(blist))
        tra_median = np.median(alist)/10
        trb_median = np.median(blist)/10
        # print("tra_median:{}".format(tra_median))
        # print("trb_median:{}".format(trb_median))
        addtwodimdict(barcode_median_dict, barcode,'TRA', tra_median)
        addtwodimdict(barcode_median_dict, barcode,'TRB', trb_median)
    
    # print(barcode_median_dict)
    rid_list = []
    with open(outfile1, "w") as out1, open(outfile2, "w") as out2:
        out1.write("clone_id,clone_type,{}\n".format(",".join(barcode_list)))
        out2.write("clone_id,clone_type,{}\n".format(",".join(barcode_list)))
        for clone in cloneid_list:
            if clone.startswith('a'):
                clone_type = "TRA"
            elif clone.startswith('b'):
                clone_type = "TRB"
            #########################
            # calculate median or sum or dict.items
            each_clone_reads_list = [int(value) for value in read_count_dict[clone].values()]
            
            # while len(each_clone_reads_list) < barcode_number:
            #     each_clone_reads_list.append(0)
            median = np.median(each_clone_reads_list)
            out1.write("{},{}".format(clone, clone_type))
            out2.write("{},{}".format(clone, clone_type))
            for barcode in barcode_list:
                if barcode in read_count_dict[clone]:
                    read_count = read_count_dict[clone][barcode]
                    if read_count >= (median/10) and read_count >= barcode_median_dict[barcode][clone_type] and read_count >= 5:
                        out1.write(",1")
                        out2.write(",{}".format(read_count))
                    else:
                        out1.write(",0")
                        out2.write(",0")
                        rid_list.append("{}_{}".format(clone,barcode))
                else:
                    out1.write(",0")
                    out2.write(",0")
                    rid_list.append("{}_{}".format(clone,barcode))
            out1.write("\n")
            out2.write("\n")
    print("output files: {}, {}, {}\n".format(outfile1, outfile2, outfile3))

    with open(input_file, "r") as f, open(outfile3, "w") as out3:
        for line in f:
            cloneid, readid, barcode = line.split()
            well = barcode.rstrip("\n")[4:]

            if cloneid + '_' + well in rid_list: 
                continue
            #######################################
            if cloneid == 'b0' or cloneid == 'b3':
                continue
            #######################################
            out3.write("{}\t{}\t{}\n".format(cloneid, readid, barcode))


if __name__ == '__main__':
    main()
