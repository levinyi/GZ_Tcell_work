import sys
import pandas as pd

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


def main():
    tcr_file = sys.argv[1]
    pre_real_tra_file = sys.argv[2]
    pre_real_trb_file = sys.argv[3]
    output = sys.argv[4]
    # tra_dict, tra_header = get_tcr_dict(pre_real_tra_file)
    # trb_dict, trb_header = get_tcr_dict(pre_real_trb_file)
    tra_table = pd.read_csv(pre_real_tra_file, sep=",")
    trb_table = pd.read_csv(pre_real_trb_file, sep=",")

    tcr_data = pd.read_csv(tcr_file, sep=",")
    tcr_data['TRA_clonotype'] = tcr_data.apply(lambda x : x['v_gene'].replace("/","") + x['cdr3'] + x['j_gene'].replace("/",""), axis=1)
    tcr_data['TRB_clonotype'] = tcr_data.apply(lambda x : x['v_gene.1'].replace("/","") + x['cdr3.1'] + x['j_gene.1'].replace("/",""), axis=1)
    print(tcr_data.head)
    # merge three files by clonotypes.
    tcr_data = pd.merge(tcr_data, tra_table, how='left', left_on=['TRA_clonotype'], right_on=['TCR_id(Real_samples_union)'])
    tcr_data = pd.merge(tcr_data, trb_table, how='left', left_on=['TRB_clonotype'], right_on=['TCR_id(Real_samples_union)'])
    print(tcr_data.head)
    tcr_data.to_csv(output, sep=",", index=False)
    #####################################
    '''
    with open(tcr_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("clonotype_id"):
                print("{},{},{}".format(line,",".join(tra_header),",".join(trb_header)))
                continue
            c = line.split(",")
            TRA_clonotype = c[4].split("*")[0].replace("/","") + c[10] + c[6].split("*")[0].replace("/","")
            TRB_clonotype = c[17].split("*")[0].replace("/","") + c[23] + c[19].split("*")[0].replace("/","")
            ## print(TRA_clonotype, TRB_clonotype)
            # print("{},{},{}".format(",".join(c),tra_dict.get(TRA_clonotype,",".join(['0']*13)),trb_dict.get(TRB_clonotype, ",".join(['0']*13))))
    '''
if __name__ == '__main__':
    main()
