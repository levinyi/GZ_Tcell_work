import sys
import os
from Bio import SeqIO
import gzip
import re


def usage():
    """docstring for usage
    python clonetype_support_umi.TSO.py G17E3L1.mixcr.out.clonotypes.TRA.txt <umi_fastq> G17E3L1.umi.count.aa.all.xls > G17E3L1.umi.count.aa.xls
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


def deal_umi_file(fastq):
    umi_dict = {}
    try:
        fq = gzip.open(fastq,"r")
        for record in SeqIO.parse(fq,"fastq"):
            umi_dict[record.id] = str(record.seq)[:8]
    except IOError:
        sys.exit("error with fastq file during deal umi function")
    finally:
        fq.close()
    return umi_dict

def main():

    clonetype_file = sys.argv[1]
    project_name = os.path.basename(clonetype_file).split(".")[0]
    
    fastq_with_umi_file = sys.argv[2]

    fastq_umi_dict = deal_umi_file(fastq_with_umi_file)
    print "finished deal_umi_fastq"

    umi_dict = {}
    cloneType_support_reads_dict = {}
    p = re.compile(r'^(TR.*)\*00\(.*')
    with open(clonetype_file, "r") as f:
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
            except IOError:
                fastq_file = project_name + '/' + project_name + '.' + cloneId + '.fastq.gz'
            else:
                fq = gzip.open(fastq_file, "r")
                if float(cloneCount) >= 3 and V and J:
                    for record in SeqIO.parse(fq, "fastq"):
                        umi = fastq_umi_dict[record.id]
                        two_dim_dict(umi_dict, V.group(1) + aaSeqCDR3 + J.group(1), umi, 1)
                        two_dim_dict(cloneType_support_reads_dict, V.group(1) + aaSeqCDR3 + J.group(1), 'reads_number', 1)
            finally:
                fq.close()

    # writ to output.
    with open(sys.argv[3], "w") as output:
        output.write("uniqueCloneType\tUMIs\treadsNumber\n")
        for clonotype in umi_dict:
            for umi in umi_dict[clonotype]:
                output.write("{0}\t{1}\t{2}\n".format(clonotype, umi, umi_dict[clonotype][umi]))

if __name__ == '__main__':
    main()
