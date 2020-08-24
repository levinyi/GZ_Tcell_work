import sys
from itertools import islice
import numpy as np


def deal_file(afile, bfile, clone_id_list):
    mixcr_dict ={}
    for each in [afile, bfile]:
        with open(each, "r") as f:
            for line in islice(f, 1, None):
                line = line.rstrip("\n")
                CloneId, Clonotype, TRV, CDR3, TRJ, ReadsNumber = line.split("\t")
                if CloneId.startswith("a"):
                    clone_type = 'TRA'
                elif CloneId.startswith("b"):
                    clone_type = "TRB"
                if CloneId not in clone_id_list:
                    clone_id_list.append(CloneId)
                mixcr_dict[CloneId] = ReadsNumber
                # print("{}\t{}\tReads-in-Mixcr\t{}".format(CloneId, clone_type, ReadsNumber))
    return mixcr_dict, clone_id_list

input1 = sys.argv[1] # Total.G208E1.TRAB.wells.count.matrix.Raw.csv
input2 = sys.argv[2] # G208E1L1.mixcr.out.clonotypes.TRA.count.txt
input3 = sys.argv[3] # G208E1L2.mixcr.out.clonotypes.TRB.count.txt

well_dict = {}
clone_id_list = []
with open(input1, "r") as f:
    for line in islice(f, 1, None):
        line = line.rstrip("\n")
        a = line.split(",")
        count_sum = np.sum([int(each) for each in a[2:]])
        if a[0] not in clone_id_list:
            clone_id_list.append(a[0])
        well_dict[a[0]] = count_sum
        # print("{}\t{}\tReads-in-Well\t{}".format(a[0], a[1], count_sum))

mixcr_dict, clone_id_dict = deal_file(input2, input3, clone_id_list)

print("cloneId,cloneType,Reads-in-Well,Reads-in-Mixcr")
for each in clone_id_list:
    if each.startswith("a"):
        clone_type = 'TRA'
    elif each.startswith("b"):
        clone_type = "TRB"

    print("{},{},{},{}".format(each, clone_type, well_dict.get(each, 0), mixcr_dict.get(each, 0)))

