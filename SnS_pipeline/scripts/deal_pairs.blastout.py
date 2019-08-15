import sys
import os
import numpy as np

blastout = sys.argv[1]

reads_dict = {}
with open(blastout, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        line = line.rstrip("\n")

        query, subject, identity, alignment_length, mismatches, gap_opens, q_start, q_end, s_start, s_end, evalue, bit_score = line.split()

        reads_dict.setdefault(query, []).append(subject)


no_use = 0
perfect_match = 0
mispairing = 0
total_reads = 0
pairs_name = []
mispairing_name = []
for each in reads_dict:
    # deal_list(reads_dict[each])
    total_reads += 1
    group_a = []
    group_b = []
    for subject in reads_dict[each]:
        if subject.endswith('a'):
            group_a.append(subject.rstrip("a"))
        elif subject.endswith('b'):
            group_b.append(subject.rstrip('b'))

    if len(group_a) == 0 and len(group_b) != 0:
        # all the group are B:
        no_use += 1
    elif len(group_a) != 0 and len(group_b) == 0:
        no_use += 1
    elif len(group_a) == 0 and len(group_b) == 0:
        print "no read find"
    elif len(group_a) != 0 and len(group_b) != 0:
        interact = list(set(group_a) & set(group_b))
        if interact:
            perfect_match +=1
            pairs_name.append(interact)
        else:
            mispairing += 1
            mispairing_name.append((group_a[0],group_b[0]))

print("no_use: %s\t%.2f%%" % (no_use, float(no_use) / total_reads*100))
print("perfect_match: %s\t%.2f%%" % (perfect_match, float(perfect_match) / total_reads*100))
print("mispairing : %s\t%.2f%%" % (mispairing, float(mispairing) / total_reads*100))
print("total_reads : %s" % total_reads)


with open(sys.argv[2],"w") as out_put:
    for each in pairs_name:
        out_put.write("{0}\n".format(''.join(each)))

with open(sys.argv[3], "w") as out_put2:
    for each in mispairing_name:
        out_put2.write("{0}\t{1}\n".format(each[0],each[1]))