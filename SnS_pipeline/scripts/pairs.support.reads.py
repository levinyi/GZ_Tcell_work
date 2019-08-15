import sys

pairs_file = sys.argv[1]
reference_file = sys.argv[2]


def deal_ref(reference_file):
    """docstring for deal_ref"""
    reference_name = []
    with open(reference_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            a = line.split()
            name = a[0].rstrip("a,b")
            reference_name.append(name)

    return reference_name

pairs_dict = {}
with open(pairs_file, "r") as f2:
    for line in f2:
        line = line.rstrip("\n")
        if line not in pairs_dict:
            pairs_dict[line] = 1
        else:
            pairs_dict[line] = pairs_dict[line] + 1

reference_list = deal_ref(reference_file)
for each in reference_list:
    print each, pairs_dict.get(each, 0)

