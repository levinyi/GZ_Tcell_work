import sys
import os


def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict



def deal_10x(afile):
    tradict = {}
    trbdict = {}
    with open(afile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("clonotype_id"):
                continue
            c = line.split(",")
            tra_clonotype = c[4].split("*")[0].replace("/", "") + c[10] + c[6].split("*")[0]
            trb_clonotype = c[17].split("*")[0].replace("/", "") + c[23] + c[19].split("*")[0]
            tradict[tra_clonotype] = c[-1]
            trbdict[trb_clonotype] = c[-1]
    return tradict, trbdict



def main():
    mixcr_file = sys.argv[1]
    tcr_type = sys.argv[2] # TRA or TRB
    vdj_files = sys.argv[3:]
    
    tradict = {}
    trbdict = {}
    sample_name_list = []
    for vdj_file in vdj_files:
        sample_name = os.path.basename(vdj_file).split(".")[0]
        sample_name_list.append(sample_name)
        tradict1, trbdict1 = deal_10x(vdj_file)
        # merge dict
        for clonotype in tradict1:
            addtwodimdict(tradict, clonotype, sample_name, tradict1[clonotype])
        for clonotype in trbdict1:
            addtwodimdict(trbdict, clonotype, sample_name, trbdict1[clonotype])

    if tcr_type == "TRA":
        vdj_dict = tradict
    elif tcr_type == "TRB":
        vdj_dict = trbdict
    else:
        sys.exit("Error: Wrong tcr types, should be :TRA or TRB.")
    
    # print(vdj_dict)

    output = open(os.path.basename(mixcr_file).rstrip("xls") + "add_TCRseq.freq.xls", "w")
    with open(mixcr_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("Clonetype"):
                output.write("{}\t{}\n".format(line,"\t".join(sample_name_list)))
                continue
            c = line.split()
            output.write("{}".format(line))
            for sample_name in sample_name_list:
                if c[0] in vdj_dict:
                    output.write("\t{}".format(vdj_dict[c[0]].get(sample_name, "NULL")))
                else:
                    output.write("\tNULL")
            output.write("\n")
    output.close()
    print("output file : {}add_TCRseq.freq.xls".format(os.path.basename(mixcr_file).rstrip("xls") ))

if __name__ == '__main__':
    main()
