import sys,os
import gzip
from Bio import SeqIO
import multiprocessing


def usage():
    print('''
    python {0} <sample_fastq> <fastqs_path> <tcr_type>

    python {0} G128E1L1.fastq.gz G128_G114_G115_TRA TRA

    '''.format(os.path.basename(sys.argv[0])))


def deal_fastq(fastq_file):
    if fastq_file.endswith(".gz"):
        handle = gzip.open(fastq_file, "r")
    else:
        handle = open(fastq_file, "r")

    adict = {}
    for record in SeqIO.parse(handle, "fastq"):
        adict[str(record.id)] = str(record.seq)
    handle.close()
    return adict


def mycallback(delta_list):
    print("I'm writing once")
    with open(sys.argv[4], 'a+') as f:
        f.write("{}\n".format("\n".join(map(str, delta_list))))


def main():
    sample_fastq = sys.argv[1]  # r/cygene/work/00.test/pipeline/Shotgun_TCRseq/bin/extract_read_from_merged_mixcr_result.pyawdata/xxx.fastq.gz.
    merged_fastq_path = sys.argv[2]  # G128_G114_G115_TRA
    tcr_type = sys.argv[3]  # TRA or TRB 
    output_file = sys.argv[4]

    if tcr_type == 'TRA':
        t_type = 'a'
    elif tcr_type == 'TRB':
        t_type = 'b'
    else:
        usage()
        sys.exit("Error: wrong tcr_type.")
    
    # PE or SE read ?
    fastq_dict = deal_fastq(sample_fastq)
    fastq_list = [merged_fastq_path+'/'+f for f in os.listdir(merged_fastq_path) if f.endswith("R1.fastq.gz")]
    output = open(output_file, "w")
    for each_fastq in fastq_list:
        cloneid = each_fastq.split(".")[1].split("_")[0]
        if each_fastq.endswith(".gz"):
            handle = gzip.open(each_fastq, "r")
        else:
            handle = open(each_fastq, "r")

        for record in SeqIO.parse(handle, "fastq"):
            if str(record.id) in fastq_dict:
                output.write("{}{}\t{}\n".format(t_type, cloneid, str(record.id)))
    output.close()
    '''
    p = multiprocessing.Pool()
    for each_fastq in fastq_list:
        p.apply_async(extract_read_id, args=(fastq_dict, each_fastq, tcr_type,), callback=mycallback)
        print("doing file :{}".format(each_fastq))
    p.close()
    p.join()
    '''

if __name__ == '__main__':
    main()

