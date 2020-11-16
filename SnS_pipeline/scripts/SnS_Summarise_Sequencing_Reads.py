import sys
import os

if len(sys.argv) == 1:
    workpath = './'
else:
    workpath = sys.argv[1]

output1 = open("Total.Count.UMI.info.summary.csv","w")
#output1.write("sample,total,match_A,match_B,six_total,match_perfect,percent\n")
adict = {}
content = []
for each in os.listdir(workpath):
    if "pair.count" in each:
        sample_name = each.split(".")[0]
        with open(each, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("total"):
                    header = [each for each in line.split()]
                    continue
                c = line.split()
                adict[sample_name] = c[0]
                content.append("{},{}".format(sample_name, ",".join(c)))
                #output1.write("{},{}\n".format(sample_name, ",".join(c)))
output1.write("sample,{}\n".format(",".join(header)))
for each in content:
    output1.write("{}\n".format(each))
output1.close()
if os.path.exists("sequencing_info.txt"):
    output2 = open("NGS_Reads.summary.csv", "w")
    output2.write("Sample,PF Reads,Cell number,read per cell,MiSeq,Targeted Reads before sequencing\n")
    with open("sequencing_info.txt", "r") as f:
        for line in f:
            line = line.rstrip("\n")
            c = line.split("\t")
            sample_name = c[0].replace(".", "")
            reads_number = adict[sample_name]
            cell_number = c[2].split()[0]
            reads_per_cell = float(reads_number)/float(cell_number)
            targeted_reads = c[20]
            output2.write("{},{},{},{},MiSeq run-1,{}\n".format(sample_name, adict[sample_name], cell_number, reads_per_cell, targeted_reads))
    output2.close()
