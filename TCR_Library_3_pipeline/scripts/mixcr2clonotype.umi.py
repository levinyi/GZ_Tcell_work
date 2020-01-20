import sys
import argparse

def usage():
    """
    try me:
    # this is for Library3 method. input file should be reads file.
    python mixcr2clonotype.umi.py mixcr.result <read dir>
    
    20200107: modify format.
    xxx: created.
    """
    

def _argparse():
    parser = argparse.ArgumentParser(description="A Python-Venn")
    parser.add_argument('-i', '--input', action='store',      dest='mixcr_file', required=True, help='mixcr result file,usually like so: G22_R2.mixcr.out.clonotypes.TRA.txt.')
    parser.add_argument('-r', '--reads', action='store',      dest='reads_file', required=True, help='you need to provide reads file')
    parser.add_argument('-umi', '--umi', action='store_true', dest='umi', default=True, help='reads with umi or without umi, default with umi.')
    parser.add_argument('-d', '--dir',   action='store',      dest='fq_dir', required=True, help='fastq file dir. a lot of fastq file in this dir.')
    parser.add_argument('-m','--merged', action='store_true', dest='merged', help='if mixcr result is composed by multi sample.')
    parser.add_argument('-l', '--list',  action='store',      dest='list_file', default='list', help='if -m paramter is used, you need to privide a sample list.')


    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
    return parser.parse_args()


def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict

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

def deal_clones_file(clones_file):
    clones_dict = {}
    with open(clones_file,"r") as f:
        for line in f:
            line = line.rstrip("\n")
            if len(line.split("\t")) >15:
                cloneId = line.split("\t")[0] # cloneId
                bestVGene = line.split("\t")[5] # allVHitsWithScore
                bestJGene = line.split("\t")[7] # allJHitsWithScore
                aaSeqCDR3 = line.split("\t")[32] #
            else:
                cloneId,cloneCount,bestVGene,bestJGene,nSeqCDR3,aaSeqCDR3,qualCDR3 = line.split("\t")
            
            addtwodimdict(clones_dict, cloneId,'TRV',bestVGene)
            addtwodimdict(clones_dict, cloneId,'TRJ',bestJGene)
            addtwodimdict(clones_dict, cloneId,'cdr3',aaSeqCDR3)

    umi_dict = {}
    cloneType_support_reads_dict = {}
    p = re.compile(r'^(TR.*)\*00\(.*')
    with open(clones_file, "r") as f:
        project_name = os.path.basename(clonetype_file).split(".")[0]
        for line in f:
            if line.startswith("cloneId"):
                continue
            c = line.split("\t")
            cloneId, cloneCount, cloneFraction, targetSequences, targetQualities, allVHitsWithScore, allDHitsWithScore, allJHitsWithScore,\
                allCHitsWithScore, allVAlignments, allDAlignments, allJAlignments, allCAlignments, nSeqFR1, minQualFR1, nSeqCDR1, minQualCDR1,\
                nSeqFR2, minQualFR2, nSeqCDR2, minQualCDR2, nSeqFR3, minQualFR3, nSeqCDR3, minQualCDR3, nSeqFR4, minQualFR4, aaSeqFR1, aaSeqCDR1,\
                aaSeqFR2, aaSeqCDR2, aaSeqFR3, aaSeqCDR3, aaSeqFR4, refPoints = line.split("\t")
            allVHitsWithScore = allVHitsWithScore.split(",")[0]
            allJHitsWithScore = allJHitsWithScore.split(",")[0]
            V = p.match(allVHitsWithScore)
            J = p.match(allJHitsWithScore)
            try:
                fastq_file = project_name + '/' + project_name + '.' + cloneId + '_R2.fastq.gz'
                fq = gzip.open(fastq_file, "r")
            except IOError:
                fastq_file = project_name + '/' + project_name + '.' + cloneId + '.fastq.gz'
                fq = gzip.open(fastq_file, "r")

            try:
                if float(cloneCount) >= 3 and V and J:
                    for record in SeqIO.parse(fq, "fastq"):
                        umi = str(record.seq[:8])
                        #two_dim_dict(umi_dict, V.group(1) + aaSeqCDR3 + J.group(1), umi, 1)
                        two_dim_dict(umi_dict, V.group(1) + targetSequences + J.group(1), umi, 1)
                        #two_dim_dict(cloneType_support_reads_dict, V.group(1) + aaSeqCDR3 + J.group(1), 'reads_number', 1)
                        two_dim_dict(cloneType_support_reads_dict, V.group(1) + targetSequences + J.group(1), 'reads_number', 1)
            finally:
                fq.close()

    output = open(sys.argv[2], "w")
    output.write("uniqueCloneType\tUMIs\treadsNumber\n")
    # print("uniqueCloneType\treadsNumber\tumiNumber")
    return clones_dict


def main():
    """docstring for main"""
    parser = _argparse()
    # check arguments dependency.
    if parser.merged  and not parser.list:
        parser.error("The -m argument requires the -l\n you need privide a list file that contains sample names.")

    with open(mixcr_file,"w") as f:
        if 


    

    output = open(sys.argv[2], "w")
    output.write("uniqueCloneType\tUMIs\treadsNumber\n")
    # print("uniqueCloneType\treadsNumber\tumiNumber")



if __name__ == '__main__' :
    main()
