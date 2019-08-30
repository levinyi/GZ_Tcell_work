import sys

def usage():
    """docstring for usage
    python $0 G33E1L1.umi.count.aa.all.xls >G33E1L1.umi.clonotype.matrix.reshape.3.xls
    """

    

def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a:{key_b: val}})
    return thedict

adict = {}
umi_list = []
clonotype_list = []
rank_dict = {}
umi_rank_dict = {}
inputfile = sys.argv[1]
with open(inputfile,"r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("unique"):
            continue
        Clonotype, UMIs, ReadsNumber = line.split("\t")
        umi_list.append(UMIs)
        clonotype_list.append(Clonotype)
        addtwodimdict(adict, UMIs, Clonotype, ReadsNumber)
        if Clonotype not in rank_dict:
            rank_dict[Clonotype] = int(ReadsNumber)
        else:
            rank_dict[Clonotype] += int(ReadsNumber)
        if UMIs not in umi_rank_dict:
            umi_rank_dict[UMIs] = int(ReadsNumber)
        else:
            umi_rank_dict[UMIs] += int(ReadsNumber)

umi_list = list(set(umi_list))
clonotype_list = list(set(clonotype_list))

'''
print "\t","\t".join(clonotype_list)
for each_umi in umi_list:
    print each_umi,"\t",
    for each_clone in clonotype_list:
        print adict[each_umi].get(each_clone, 0),"\t",
    print ""
'''
for each_umi in umi_list:
    for each_clone in clonotype_list:
        print("{}\t{}\t{}\t{}\t{}".format(each_umi, each_clone, adict[each_umi].get(each_clone, 0), umi_rank_dict[each_umi],rank_dict[each_clone]))
