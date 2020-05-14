import sys

input_file = sys.argv[1]

with open(input_file, "r") as f:
    for line in f:
        if len(line.split("\t")) == 15:
            # if "insert" in line.split("\t")[5] or "del" in line.split("\t")[5]:
            #     chrom = line.split("\t")[0]
            #     start = int(line.split("\t")[1]) -2
            #     end = int(line.split("\t")[2])-1
            # else:
            chrom = line.split("\t")[0]
            start = int(line.split("\t")[1]) -1
            end = int(line.split("\t")[2])
        else:
            chrom = line.split("\t")[6]
            start = int(line.split("\t")[7]) - 1
            end = int(line.split("\t")[8])
        print("{}\t{}\t{}".format(chrom, start, end))