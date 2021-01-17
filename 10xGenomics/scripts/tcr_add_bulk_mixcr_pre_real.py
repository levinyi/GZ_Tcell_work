import sys

def get_tcr_dict(afile):
    adict = {}
    with open(afile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("TCR_id"):
                header=line.split(",")[1:]
                continue
            c = line.split(",")
            adict[c[0]] = ",".join(c[1:])
    return adict, header




tcr_file = sys.argv[1]
pre_real_tra_file = sys.argv[2]
pre_real_trb_file = sys.argv[3]

tra_dict, tra_header = get_tcr_dict(pre_real_tra_file)
trb_dict, trb_header = get_tcr_dict(pre_real_trb_file)


with open(tcr_file, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("clonotype_id"):
            print("{},{},{}".format(line,",".join(tra_header),",".join(trb_header)))
            continue
        c = line.split(",")
        TRA_clonotype = c[4].split("*")[0].replace("/","") + c[10] + c[6].split("*")[0].replace("/","")
        TRB_clonotype = c[17].split("*")[0].replace("/","") + c[23] + c[19].split("*")[0].replace("/","")
        #print(TRA_clonotype, TRB_clonotype)
        print("{},{},{}".format(",".join(c),tra_dict.get(TRA_clonotype,",".join(['0']*13)),trb_dict.get(TRB_clonotype, ",".join(['0']*13))))
