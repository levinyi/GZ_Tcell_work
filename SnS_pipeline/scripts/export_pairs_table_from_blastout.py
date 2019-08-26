import sys
from Bio import SeqIO
import threading
import Queue
import time


def usage():
    """
    python ../scripts/export_pairs_table_from_blastout.py G34E3L1.TRA.blast.out G34E3L1.TRB.blast.out ../database/SplitCDR3JOligos.zip.fa
    """

def deal_zip_file(zip_file):
    zip_dict = {}
    with open(zip_file, "r") as f:
        for record in SeqIO.parse(f,"fasta"):
            name = str(record.id).rstrip(".zip")
            zip_dict[name] = str(record.seq)
    return zip_dict

def deal_blast_out(blast_file):
    blast_out_dict = {}

    with open(blast_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            query, subject, identity, alignmentlength, mismatches, gap_opens, q_start, q_end, s_start, s_end, evalue, bit_score = line.split("\t")
            blast_out_dict[query] = subject
    return blast_out_dict

def main():
    blast_out_TRA = sys.argv[1]
    blast_out_TRB = sys.argv[2]
    reference_zip = sys.argv[3] 
    output_file = open("export_pairs_table.out.xls", "w")
    log_file = open("blastn.pairs.log", "w")

    log_file.write("start time: {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ))
    blast_out_dict_A = deal_blast_out(blast_out_TRA)
    blast_out_dict_B = deal_blast_out(blast_out_TRB)
    zip_dict = deal_zip_file(reference_zip)
    log_file.write("end time :{}\n".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    log_file.write("time consuming: XXX\n")
    
    # multi process:
    '''
    log_file.write("start time : {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ))
    q = Queue.Queue()
    result = list()
    t1 = threading.Thread(target=deal_blast_out, name='thread1', args=(blast_out_TRA,))
    t2 = threading.Thread(target=deal_blast_out, name='thread2', args=(blast_out_TRB,))
    t3 = threading.Thread(target=deal_zip_file,  name='thread3', args=(reference_zip,))

    t1.start()
    print("Start deal with blast TRA file 1")
    t2.start()
    print("Start deal with blast TRB file 2")
    t3.start()
    print("Start deal with reference file")

    t1.join()
    t2.join()
    t3.join()

    while not q.empty():
        print "waiting ..."
        result.append(q.get())
    print "there are %s results in this threading." % len(result)

    for item in result:
        if len(item) == 2:
            if item[1] == blast_out_TRA:
                print "successfully get blast TRA dict "
                blast_out_dict_A = item[0]
            elif item[1] == blast_out_TRB:
                print "successfully get blast TRB dict "
                blast_out_dict_B = item[0]
            elif item[1] == reference_zip:
                print "successfully get reference dict"
                zip_dict = item[0]
    log_file.write("finished multi process at : {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    log_file.write("total time consuming:XXX ")
    '''
    # merge 2 dicts
    total_mapped_reads = set(blast_out_dict_A.keys()) | set(blast_out_dict_B.keys())
    log_file.write("total read id is : {0}\n".format(len(total_mapped_reads)))
    for each in total_mapped_reads:
        # print each
        if each in blast_out_dict_A:
            name_A = blast_out_dict_A[each].rstrip("a")
        else:
            name_A = "NULL"

        if each in blast_out_dict_B:
            name_B = blast_out_dict_B[each].rstrip("b")
        else:
            name_B = "NULL"

        if name_A != 'NULL' and name_B != 'NULL' and name_A == name_B:
            flag = 1
        elif name_A != 'NULL' and name_B != 'NULL' and name_A != name_B:
            flag = 2
        else:
            flag = 0

        output_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(each, name_A, name_B, zip_dict.get(name_A,"NULL"), zip_dict.get(name_B,"NULL"),flag))

    log_file.write("All finished at : {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ))
    output_file.close()
    log_file.close()

if __name__ == '__main__':
    main()