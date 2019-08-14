import sys
import os
from Bio import SeqIO
import gzip


def usage():
    '''
    python deal_cutadapt_log_PE.py 11_R1.cutadapt.log /cygene/work/20.RP1G11/data/11_R2.fq.gz /cygene/work/20.RP1G11/analysis
    '''


def deal_fastq(fastq):
    adict = {}
    for record in SeqIO.parse(fastq, "fastq"):
        adict[record.id] = record
    return adict


logfile = sys.argv[1]
PE_file = sys.argv[2]
result_dir = sys.argv[3]

if PE_file.endswith("gz"):
    PE_header = gzip.open(PE_file)
else:
    PE_header = open(PE_file)

PE_dict = deal_fastq(PE_header)


file_name = os.path.basename(logfile).split(".")[0]
sample, R = file_name.split("_")
print sample, R

with open(logfile, "r") as f, open(result_dir + '/' + file_name + ".p1.fq", "w") as OUT1, open(result_dir + '/' + file_name + '.p2.fq', "w") as OUT2:
    if R == 'R1':
        with open(result_dir + '/' + sample + "_R2.p1.fq", "w") as PE_OUT1, open(result_dir + '/' + sample + "_R2.p2.fq", "w") as PE_OUT2:
            for line in f:
                line = line.rstrip("\n")
                c = line.split()
                if c[2] == '-1':
                    continue
                elif c[3] == '0': # adapter from start.
                    OUT2.write("@%s\n%s\n+\n%s\n" % (c[0], c[6], c[9]))
                    SeqIO.write(PE_dict[c[0]], PE_OUT2, "fastq")
                else:
                    if c[4] == '451':
                        OUT1.write("@%s\n%s\n+\n%s\n" % (c[0], c[5], c[8]))
                        SeqIO.write(PE_dict[c[0]], PE_OUT1, "fastq")
                    else:
                        OUT1.write("@%s\n%s\n+\n%s\n" % (c[0], c[5], c[9]))
                        SeqIO.write(PE_dict[c[0]], PE_OUT1, "fastq")
                        OUT2.write("@%s\n%s\n+\n%s\n" % (c[0], c[7], c[11]))
                        SeqIO.write(PE_dict[c[0]], PE_OUT2, "fastq")
    elif R == 'R2':
        with open(result_dir + '/' + sample + "_R1.p1.fq", "w") as PE_OUT1, open(result_dir + '/' + sample + "_R1.p2.fq", "w") as PE_OUT2:
            for line in f:
                line = line.rstrip("\n")
                c = line.split()
                if c[2] == '-1':
                    continue
                elif c[3] == '0':
                    OUT2.write("@%s\n%s\n+\n%s\n" % (c[0], c[6], c[9]))
                    SeqIO.write(PE_dict[c[0]], PE_OUT2, "fastq")
                else:
                    if len(c[5]) >= 30:
                        OUT1.write("@%s\n%s\n+\n%s\n" % (c[0], c[5], c[9]))
                        SeqIO.write(PE_dict[c[0]], PE_OUT1, "fastq")
                    OUT2.write("@%s\n%s\n+\n%s\n" % (c[0], c[7], c[11]))
                    SeqIO.write(PE_dict[c[0]], PE_OUT2, "fastq")
print("finished")
