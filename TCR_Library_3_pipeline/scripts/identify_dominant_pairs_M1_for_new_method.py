import sys
freq_file = sys.argv[1]
output_file = sys.argv[2]

clones_A_file = sys.argv[3]
clones_B_file = sys.argv[4]

def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict

def deal_clones_file(clones_file,flag):
    clones_dict = {}
    with open(clones_file,"r") as f:
        for line in f:
            line = line.rstrip("\n")
            if len(line.split("\t")) >15:
                cloneId = line.split("\t")[0] # cloneId
                bestVGene = line.split("\t")[5].split("*")[0] # allVHitsWithScore
                bestJGene = line.split("\t")[7].split("*")[0] # allJHitsWithScore
                aaSeqCDR3 = line.split("\t")[32] # 
            else:
                cloneId,cloneCount,bestVGene,bestJGene,nSeqCDR3,aaSeqCDR3,qualCDR3 = line.split("\t")
            cloneId = flag + cloneId
            addtwodimdict(clones_dict, cloneId,'TRV',bestVGene)
            addtwodimdict(clones_dict, cloneId,'TRJ',bestJGene)
            addtwodimdict(clones_dict, cloneId,'cdr3',aaSeqCDR3)
    return clones_dict

clones_A_dict = deal_clones_file(clones_A_file,'a')
clones_B_dict = deal_clones_file(clones_B_file,'b')

sum_dict = {}
with open(freq_file,"r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.endswith("freq"):
            continue
        a,b,x = line.split()
        if not a.startswith("a"): # for new method
            a = 'a' + str(a)
            b = 'b' + str(b)
        if int(x) >=3:
            if a not in sum_dict:
                sum_dict[a] = int(x)
            else:
                sum_dict[a] += int(x)
            if b not in sum_dict:
                sum_dict[b] = int(x)
            else:
                sum_dict[b] +=int(x)


with open(freq_file,"r") as f, open(output_file,"w") as op_file:
    op_file.write("Alpha Clone ID\tTRAV\tCDR3Alpha\tTRAJ\tBeta Clone Id\tTRBV\tCDR3Beta\tTRBJ\tRead Count\tA_total_reads\tB_total_reads\tM1A Score\tM1B Score\tM1 Score\n")
    for line in f:
        line = line.rstrip("\n")
        if line.endswith("freq"):
            continue
        a,b,x = line.split()
        if int(x) >=3:
            if a.startswith("a"):
                M1A = int(x)/float(sum_dict[a])
                M1B = int(x)/float(sum_dict[b])
                M1 = (int(x)/float(sum_dict[a]))*(int(x)/float(sum_dict[b]))
            else:
                a = 'a' + str(a)
                b = 'b' + str(b)
                M1A = int(x)/float(sum_dict[a])
                M1B = int(x)/float(sum_dict[b])
                M1 = (int(x)/float(sum_dict[a]))*(int(x)/float(sum_dict[b]))

            op_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(a,clones_A_dict[a]['TRV'], clones_A_dict[a]['cdr3'], clones_A_dict[a]['TRJ'],b, clones_B_dict[b]['TRV'], clones_B_dict[b]['cdr3'], clones_B_dict[b]['TRJ'],x,sum_dict[a],sum_dict[b],M1A,M1B,M1))

