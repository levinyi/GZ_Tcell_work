import sys
from itertools import islice
import math


def usage():
    """docstring for usage"""
    """
    
    python  reshape2frequency.py G41E1L2.merged.umi.count.reshape.xls >G41E1L2.raw.freq.xls
    """


def dict_freqency(adict):
    total_sum = sum((int(k) for k in adict.values()))
    new_dict = {}
    for each in adict:
        new_dict[each] = float(adict[each])/total_sum
    return new_dict


input_file = sys.argv[1]  # G39E1L2.merged.umi.count.reshape.xls
a_dict = {}
molecule_dict = {}
with open(input_file, "r") as f:
    for line in islice(f, 1, None):  # skip header line.
        line = line.rstrip("\n")
        Clonotype, TRV, CDR3, TRJ, UMIcount, ReadsNumber = line.split("\t")
        a_dict[Clonotype] = ReadsNumber
        molecule_dict[Clonotype] = UMIcount


new_dict = dict_freqency(a_dict)
print("Clonotype\tSupporting_Reads\tReads_Frequency\tFrequency_Log10\tMolecules")
for k in new_dict:
    print("{}\t{}\t{}\t{}\t{}".format(k, a_dict[k], new_dict[k], math.log(new_dict[k], 10), molecule_dict[k]))
