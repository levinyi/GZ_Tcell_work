import sys
import gzip
import os
from Bio import SeqIO


def usage():
    print("""Usage:
    python {0} <barcode file> <mixcr file> <TRA/TRB> 

    Example:
    python {0} G87E2L1.barcode.read.txt G87E2L1.mixcr.out.clonotypes.TRA.txt TRA
    
    updates:
    20200607: updates.
    """.format(os.path.basename(sys.argv[0])))


def check_input_file():
    """docstring for check_input_file"""
    barcode_filename = os.path.basename(sys.argv[1]).split(".")[0]
    mixcr_filename = os.path.basename(sys.argv[2]).split(".")[0]
    if barcode_filename != mixcr_filename:
        sys.stderr("warning: input {} not match with {}".format(sys.argv[1],sys.argv[2]))


def deal_bc_readid(afile):
    """docstring for deal_bc_readid"""
    bc_read_dict = {}
    with open(afile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            bc_name, read_id = line.split()
            bc_read_dict[read_id] = bc_name
    return bc_read_dict


def get_TCR_symbol(tcr_type):
    """docstring for get_TCR_symbol"""
    if tcr_type == 'TRA':
        tcr_symbol = 'a'
    elif tcr_type == 'TRB':
        tcr_symbol = 'b'
    else:
        sys.exit("wrong TCR_type sys.argv[3]\n")
    return tcr_symbol


def main():
    """docstring for main"""
    if len(sys.argv) != 4:
        usage()
        sys.exit("Error: You must put 3 input arguments.")
    bc_readsId_file = sys.argv[1]
    mixcr_result_file = sys.argv[2]  # G84E2L2.mixcr.out.TR[A|B].txt
    tcr_type = sys.argv[3]    # TRA or TRB
    basename = os.path.basename(mixcr_result_file).split(".")[0]

    # check_input_file()
    bc_read_dict = deal_bc_readid(bc_readsId_file)

    tcr_symbol = get_TCR_symbol(tcr_type)

    # deal with mixcr file and output the result.
    with open(mixcr_result_file, "r") as f1, open(basename+"."+tcr_type+'.clone.reads.barcode.txt', "w") as output:
        for line in f1:
            if line.startswith("cloneId"):
                continue
            cloneId, cloneCount, cloneFraction, targetSequences, targetQualities, allVHitsWithScore, allDHitsWithScore, allJHitsWithScore,\
            allCHitsWithScore, allVAlignments, allDAlignments, allJAlignments, allCAlignments, nSeqFR1, minQualFR1, nSeqCDR1, minQualCDR1,\
            nSeqFR2, minQualFR2, nSeqCDR2, minQualCDR2, nSeqFR3, minQualFR3, nSeqCDR3, minQualCDR3, nSeqFR4, minQualFR4, aaSeqFR1, aaSeqCDR1,\
            aaSeqFR2, aaSeqCDR2, aaSeqFR3, aaSeqCDR3, aaSeqFR4, refPoints = line.split("\t")
            if eval(cloneCount) >= 3:
                fastq = basename+"/"+basename+'.'+cloneId+'_R1.fastq.gz'
                fastq_file = gzip.open(fastq, "r")
                for record in SeqIO.parse(fastq_file, "fastq"):
                    if bc_read_dict.get(str(record.id)):
                        output.write("{}{}\t{}\t{}\n".format(tcr_symbol, cloneId, str(record.id), bc_read_dict[str(record.id)]))
                fastq_file.close()


if __name__ == '__main__':
    main()

