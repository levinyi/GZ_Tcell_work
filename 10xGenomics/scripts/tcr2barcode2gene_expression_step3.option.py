import sys

filtered_file = sys.argv[1]
raw_file = sys.argv[2]

header = []
with open(filtered_file, "r") as f:
    line = f.readline() # first line
    line = line.rstrip("\n")
    a = line.split()
    for each in a :
        header.append(each)
header = list(set(header))
with open(raw_file, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("barcode"): # first line
            a = line.split()
            new_header = []
            for each in a:
                if each in header:
                    new_header.append("in_filter")
                else:
                    new_header.append("only_in_raw")
            print("barcode_note\t{}".format("\t".join(new_header)))
        print("{}".format(line))
