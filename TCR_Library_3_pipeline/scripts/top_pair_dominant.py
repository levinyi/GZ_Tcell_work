import sys


def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict

clone_dict = {}
with open(sys.argv[1]) as f:
    for line in f:
        line = line.rstrip("\n")
        if line.endswith("freq"):
            continue
        a, b, freq = line.split()
        addtwodimdict(clone_dict, a, b, freq)

big_dict = {}
clone_reads_dict ={}
for k, v in clone_dict.items():
    if len(v) == 1:
        for i, j in v.items():
            y_top1 = float(j) / float(j)
            addtwodimdict(big_dict, 'top1', k, y_top1)
            addtwodimdict(clone_reads_dict, 'top1', k, float(j))
            addtwodimdict(big_dict, 'top2', k, y_top1)
            addtwodimdict(clone_reads_dict, 'top2', k, float(j))
            addtwodimdict(big_dict, 'top3', k, y_top1)
            addtwodimdict(clone_reads_dict, 'top3', k, float(j))
    elif len(v) > 1:
        # print k,v
        value = []
        for i, j in v.items():
            value.append(int(j))
        total_reads = sum(value)

        y_top1 = max(value)
        value.remove(y_top1)
        top1_reads = y_top1
        y_top2 = max(value)
        value.remove(y_top2)
        top2_reads = top1_reads+y_top2
        y_top1_p = top1_reads / float(total_reads)
        y_top2_p = top2_reads / float(total_reads)

        if len(value) == 0:
            top3_reads = top2_reads
            y_top3_p = y_top2_p
        else:
            y_top3 = max(value)
            top3_reads = top2_reads+ y_top3
            y_top3_p = top3_reads / float(total_reads)
        addtwodimdict(big_dict, 'top1', k, y_top1_p)
        addtwodimdict(clone_reads_dict, 'top1', k, top1_reads)
        addtwodimdict(big_dict, 'top2', k, y_top2_p)
        addtwodimdict(clone_reads_dict, 'top2', k, top2_reads)
        addtwodimdict(big_dict, 'top3', k, y_top3_p)
        addtwodimdict(clone_reads_dict, 'top3', k, top3_reads)
        
        #print "%s\t%s\t%s\t%.2f\t%s\t%.2f\t%s\t%.2f" % (index, k, top1_reads, y_top1_p, top2_reads, y_top2_p, top3_reads, y_top3_p)

print("CloneId\tCAT\tRank\tReads\tValue")

cloneid_rank_dict = {}
for cat, v in big_dict.items():
    index = 1
    for each in sorted(v.items(),key=lambda item:item[1],reverse=True):
        addtwodimdict(cloneid_rank_dict,cat,each[0],index)
        index +=1

for cat, v in big_dict.items():
    index = 1
    for cloneid,value in v.items():
        print("%s\t%s\t%s\t%s\t%s"%(cloneid,cat,cloneid_rank_dict[cat][cloneid],clone_reads_dict[cat][cloneid],value))


