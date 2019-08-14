import sys

M1_file = sys.argv[1]
flag = sys.argv[2] # 1 or 2; 1 means filter M1=1,2 means not filter.
all_gene_list = ["TRAV1-2","TRAV1-1","TRAV2","TRAV3","TRAV4","TRAV5","TRAV6","TRAV7","TRAV8-6","TRAV8-4","TRAV8-7","TRAV8-1","TRAV8-3","TRAV8-2","TRAV9-1","TRAV9-2","TRAV10","TRAV11","TRAV12-1","TRAV12-2", "TRAV12-3","TRAV13-1","TRAV13-2","TRAV14DV4","TRAV16","TRAV17","TRAV18","TRAV19","TRAV20","TRAV21","TRAV22","TRAV23DV6","TRAV24","TRAV25","TRAV26-2","TRAV26-1","TRAV27","TRAV29DV5","TRAV30", "TRAV34","TRAV35","TRAV36DV7","TRAV38-2DV8","TRAV38-1","TRAV39","TRAV40","TRAV41","TRBV1","TRBV2","TRBV3-1","TRBV3-2","TRBV4-1","TRBV4-3","TRBV4-2","TRBV5-5","TRBV5-3","TRBV5-1","TRBV5-4", "TRBV5-7","TRBV5-8","TRBV5-6","TRBV6-4","TRBV6-2","TRBV6-9","TRBV6-7","TRBV6-8","TRBV6-1","TRBV6-6","TRBV6-3","TRBV6-5","TRBV7-9","TRBV7-4","TRBV7-1","TRBV7-2","TRBV7-8","TRBV7-7","TRBV7-3", "TRBV7-6","TRBV9","TRBV10-1","TRBV10-2","TRBV10-3","TRBV11-1","TRBV11-2","TRBV11-3","TRBV12-4","TRBV12-5","TRBV12-3","TRBV12-1","TRBV12-2","TRBV13","TRBV14","TRBV15","TRBV16","TRBV17","TRBV18", "TRBV19","TRBV20-1","TRBV21-1","TRBV23-1","TRBV24-1","TRBV25-1","TRBV26","TRBV27","TRBV28","TRBV29-1","TRBV30"]

def store_dict(adict,k,v):
    if k not in adict:
        adict[k] = v
    else:
        adict[k] = adict[k] + v
    return adict

gene_dict = {}
with open(M1_file,"r") as f:
    for line in f:
        if line.startswith("Alpha"):
            continue
        line = line.rstrip("\n")
        AlphaCloneID,TRAV,CDR3Alpha,TRAJ,BetaCloneId,TRBV,CDR3Beta,TRBJ,ReadCount,A_total_reads,B_total_reads,M1AScore,M1BScore,M1Score = line.split("\t")
        #deal_TRAV
        #print TRAV,TRBV,ReadCount
        TRAV = TRAV.split("*")[0].replace(" ","")
        TRBV = TRBV.split("*")[0].replace(" ","")
        if flag == '1':
            if float(M1Score) == 1:
                store_dict(gene_dict,TRAV,int(ReadCount))
                store_dict(gene_dict,TRBV,int(ReadCount))
        else:
            store_dict(gene_dict,TRAV,int(ReadCount))
            store_dict(gene_dict,TRBV,int(ReadCount))

for each in all_gene_list:
    if each.startswith("TRA"):
        print each,gene_dict.get(each,0),"TRA"
    else:
        print each,gene_dict.get(each,0),"TRB"
