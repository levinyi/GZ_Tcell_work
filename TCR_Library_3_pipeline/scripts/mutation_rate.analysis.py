""" this is the main room"""
import sys
import os


def main():
    """docstring for main"""
    work_dir = os.path.basename('./')
    data_dir = work_dir
    if os.path.exists(data_dir):
        sys.exit("no such file")
    analysis_dir = ''
    sub_anls_dir = analysis_dir + '/' +''
    if not os.path.exists('muta')


    # write to the work.sh for each 
    for each in data_file:
        with open(path + "/work.sh","w") as f:
            f.write("python /cygene/work/27.G13/G13AB/deal_umi.py ../../../data/G22E5L1_R2.fq.gz  ../../G22_R1/*.fastq.gz")
            f.write("perl /home/xiaofan/for_shiyi/UTR/02.merge_umi.pl")
            f.write("ls *.umi.cons |while read line;do name=`ls $line |awk -F '/' '{print$NF}'`;echo \"less $line |awk '{print\">\"\$2\"\n\"\$4}' >$name.fa\";done |sh")
            f.write("")


if __name__ == '__main__':
    main()
