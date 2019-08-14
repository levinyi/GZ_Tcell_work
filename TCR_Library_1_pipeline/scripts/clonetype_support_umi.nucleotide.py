import sys
import os
from Bio import SeqIO
import gzip
import re


def usage():
    """docstring for usage
    python clonetype_support_umi.py G17E3L1.mixcr.out.clonotypes.TRA.txt G17E3L1.umi.count.aa.all.xls > G17E3L1.umi.count.aa.xls
    """


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


def main():

    clonetype_file = sys.argv[1]

    umi_dict = {}
    cloneType_support_reads_dict = {}
    p = re.compile(r'^(TR.*)\*00\(.*')
    with open(clonetype_file, "r") as f:
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

    # write to file 
    for clonotype in umi_dict:
        for umi in umi_dict[clonotype]:
            output.write("%s\t%s\t%s\n" % (clonotype, umi, umi_dict[clonotype][umi]))
    output.close()

    '''
    for clonotype in umi_dict:
        if len(umi_dict[clonotype]) >= 3:  # we just want to have that each clonotype has more than 3 umis.
            umi_count = 0
            supporting_reads = 0
            for umi in umi_dict[clonotype]:
                if umi_dict[clonotype][umi] >= 3:  # we just want each umi has more than 3 reads.
                    output.write("%s\t%s\t%s\n" % (clonotype, umi, umi_dict[clonotype][umi]))
                    umi_count += 1
                    supporting_reads += 1
                else:
                    output.write("%s\t%s\t%s\n" % (clonotype, umi, umi_dict[clonotype][umi]))
            if umi_count >= 3:
                print("%s\t%s\t%s" % (clonotype, supporting_reads, umi_count))

        else:
            # all to detail
            for umi in umi_dict[clonotype]:
                output.write("%s\t%s\t%s\n" % (clonotype, umi, umi_dict[clonotype][umi]))
    output.close()
    '''
if __name__ == '__main__':
    main()
