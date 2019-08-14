import sys
import os


def usage():
    '''
    step 1: python /cygene/work/00.test/pipeline/TCR_AB_type_pipeline/scripts/pairs2umi.py  G13E3L8.pairs > G13E3L8.pairs.umi.txt 
    step 2: python /cygene/work/00.test/pipeline/TCR_AB_type_pipeline/scripts/pairs2umi.merge.filter.py G13E3L8.pairs.umi.txt >G13E3L8.pairs.umi.total.txt
    '''


def two_dim_dict(thedict, key_a, key_b, value):
    if key_a in thedict:
        if key_b in thedict[key_a]:
            value = thedict[key_a][key_b] + value
            thedict[key_a].update({key_b: value})
        else:
            thedict[key_a].update({key_b: value})
    else:
        thedict.update({key_a: {key_b: value}})
    return thedict


def deal_umi_file(afile):
    adict = {}
    with open(afile,"r") as f:
        for line in f:
            line = line.rstrip("\n")
            pairs,readid,umi = line.split("\t")
            two_dim_dict(adict,pairs,umi,1)
    return adict

def main():
    input_file = sys.argv[1]

    clone2umi_dict = deal_umi_file(input_file)
    new_dict = {}
    with open(sys.argv[2],"w") as f1:
        for pairs in clone2umi_dict:
            for umi in clone2umi_dict[pairs]:
                f1.write("%s\t%s\t%s\n"% (pairs,umi,clone2umi_dict[pairs][umi]))
                if clone2umi_dict[pairs][umi] >=3:
                    two_dim_dict(new_dict,pairs,umi,clone2umi_dict[pairs][umi])
    # print new_dict
    with open(sys.argv[3],"w") as f2:
        for pairs in new_dict:
            umi_sup_reads = 0
            # if len(new_dict[pairs])>=3:
            for umi in new_dict[pairs]:
                umi_sup_reads  += new_dict[pairs][umi]

            f2.write("%s\t%s\t%s\n"%(pairs,len(new_dict[pairs]),umi_sup_reads))
            # print("%s\t%s\t%s"%(pairs,len(new_dict[pairs]),umi_sup_reads))

if __name__ == '__main__':
    main()
