import sys
import os

def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict


def deal_freq_file(freq_file):
    adict = {}
    clonotype_list = []
    with open(freq_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("CloneId"):
                continue
            CloneId, Clonotype, Reads, Frequency, FreqLog10 = line.split("\t")
            if Clonotype not in adict:
                adict[Clonotype] = Frequency
                clonotype_list.append(Clonotype)
    return adict, clonotype_list


def main():
    afile = sys.argv[1]
    sample_name_list = sys.argv[2:]
    a_big_list = []
    a_big_dict = {}
    for sample_name in sample_name_list:
        print(sample_name)
        freq_file = sample_name + '.filtered2reads.freq.xls'
        print(freq_file)
        if not os.path.exists(freq_file):
            sys.exit("file not exists!")
        clonotype_dict, clonotype_list = deal_freq_file(freq_file)

        # merge dict 
        for clonotype in clonotype_dict:
            addtwodimdict(a_big_dict, clonotype, sample_name, clonotype_dict[clonotype] )
        # merge list
        for clonotype in clonotype_list:
            a_big_list.append(clonotype)
    
    # iterate list
    output_file = open("_VS_".join(sample_name_list) + '.add.freq.xls', "w")
    output_file.write("Clonetype\t{}\n".format("\t".join(sample_name_list)))
    
    with open(afile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            c = line.split("\t")
            output_file.write("{}".format(line))
            for sample in sample_name_list:
                if c[0] not in a_big_dict:
                    output_file.write("\t1e-10".format())
                else:
                    output_file.write("\t{}".format(a_big_dict[c[0]].get(sample, "1e-10")))
            output_file.write("\n")

    #for clonotype in set(a_big_list):
    #     output_file.write("{}".format(clonotype))
    #    for sample in sample_name_list:
     #       output_file.write("\t{}".format(a_big_dict[clonotype].get(sample, "1e-10")))
     #   output_file.write("\n")

if __name__ == '__main__':
    main()
