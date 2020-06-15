import sys

input_file = "consensus_annotations.csv"

clonotype_list = []
with open(input_file, "r") as f:
    for line in f:
        if line.startswith("clonotype_id"):
            continue
        line = line.rstrip("\n")
        clonotype_id,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis = line.split(",")
        clonotype_list.append(v_gene+j_gene+cdr3)

with open(input_file, "r") as f:
    for line in f:
        if line.startswith("clonotype_id"):
            print("{},shared_YN,Number_of_sharing".format(line.rstrip("\n")))
            continue
        line = line.rstrip("\n")
        clonotype_id,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis = line.split(",")
        n = clonotype_list.count(v_gene+j_gene+cdr3)
        if n == 1:
            YN = "NO"
        else:
            YN = "YES"
        print("{},{},{}".format(line, YN, n))
