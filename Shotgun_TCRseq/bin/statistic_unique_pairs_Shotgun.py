import sys
from itertools import islice
from collections import Counter

total_pairs_number = 0
pairs_list = []
cloneid_list = []
with open(sys.argv[1], "r") as f:
    for line in islice(f, 1, None):
        line = line.rstrip("\n")
        
        total_pairs_number += 1
        a = line.split("\t")
        if line.startswith("a"):
            tra_id = a[0]
            trb_id = a[3]
        elif line.startswith("NULL"):
            tra_id = a[18]
            trb_id = a[15]
        else:
            sys.exit("Error: wrong input file")
        cloneid_list.append(tra_id)
        cloneid_list.append(trb_id)
        pairs_list.append(tra_id+'_'+trb_id)

result = Counter(cloneid_list)

A1B1 = 0
A2B1 = 0
A1B2 = 0
A2B2 = 0
other = 0
for each_pair in pairs_list:
    a, b = each_pair.split("_")

    if result[a] == 1 and result[b] == 1 :
        # print("A1B1\t{}_{}".format(a,b))
        A1B1 += 1
    elif result[a] == 2 and result[b] == 1:
        # print("A2B1\t{}_{}".format(a,b))
        A2B1 += 1
    elif result[a] == 1 and result[b] == 2:
        # print("A1B2\t{}_{}".format(a,b))
        A1B2 += 1
    elif result[a] == 2 and result[b] == 2:
        # print("A2B2\t{}_{}".format(a,b))
        A2B2 += 1
    else :
        # print("other\t{}_{}".format(a,b))
        other += 1


print("A1B1\t{}\t{}".format(A1B1, A1B1/float(total_pairs_number)))
print("A2B1\t{}\t{}".format(A2B1, A2B1/float(total_pairs_number)))
print("A1B2\t{}\t{}".format(A1B2, A1B2/float(total_pairs_number)))
print("A2B2\t{}\t{}".format(A2B2, A2B2/float(total_pairs_number)))
print("others\t{}\t{}".format(other, other/float(total_pairs_number)))
print("total_pairs\t{}".format(total_pairs_number))
