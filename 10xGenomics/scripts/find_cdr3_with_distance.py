import sys
import json
sys.path.append("/cygene/script/python_packages")
import edit_distance

input_file = "filtered_contig_annotations.csv"
input_cdr3list = sys.argv[1]

cdr3list = []
with open(input_cdr3list, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        cdr3list.append(line)

match_dict = {}
with open(input_file, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("barcode"):
            print("target_cdr3,edit_distance_less_than1,{}".format(line))
            continue
        barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,raw_consensus_id = line.split(",")
        if chain == "TRB" or chain == "TRA":
            for each in cdr3list:
                if edit_distance.minEditDist(cdr3, each) <= 1:
                    print("{},{},{}".format(each,cdr3,line))
                    match_dict.setdefault(each,[]).append(barcode)
#print(json.dumps(match_dict,indent=1))
