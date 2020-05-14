import sys
import re
from Bio import SeqIO


def deal_genome(genome_file):
    genome_dict = {}
    for record in SeqIO.parse(genome_file, "fasta"):
        genome_dict[record.id] = str(record.seq).upper()
    return genome_dict

genome_file = sys.argv[1]
myfile = sys.argv[2]  # data1.xls
indel_info = sys.argv[3]
genome_dict = deal_genome(genome_file)


# pattern = re.compile(r'chr(.*):g\.(\d+_?\d+)(.*)')
pattern = re.compile(r'(chr.*):g\.(\d+_?\d+)(.*)')  # chr9:g.139564754_139564755insC
with open(myfile, "r") as f, open(indel_info, "w") as f2:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("OM_ID"):
            continue
        transcript_ref = line.split("\t")[2]
        mutation_info = line.split("\t")[9]
        CodingDNAchange = line.split("\t")[10]
        match = pattern.match(mutation_info)
        if match:
            chromsome, position, mutation_type = match.groups()
            if ">" in mutation_type:
                position = int(position)
                if genome_dict[chromsome][position-1] == mutation_type[0]:
                    oldseq = genome_dict[chromsome][position-5 :position+10]
                    truenewseq = genome_dict[chromsome][:position-1] + mutation_type[-1] + genome_dict[chromsome][position:]
                    testnewseq = truenewseq[position-5: position+10]
                    print(mutation_info, chromsome, position, mutation_type,"yes",genome_dict[chromsome][position-1], oldseq, testnewseq)
                    genome_dict[chromsome] = truenewseq
                else:
                    print(mutation_info, chromsome, position, mutation_type,"NO",genome_dict[chromsome][position-1], mutation_type[0])
            elif "delins" in mutation_type:  # chr17:g.7578476_7578504del
                pos1, pos2 = position.split("_")
                pos1 = int(pos1)
                pos2 = int(pos2)

                oldseq = genome_dict[chromsome][pos1-5: pos2+10]
                
                del_num = pos2 - pos1 + 1
                insseq = re.match(r'.*delins(.*)', mutation_type)

                if del_num == len(insseq.group(1)):
                    truenewseq = genome_dict[chromsome][:pos1-1] + insseq.group(1) + genome_dict[chromsome][pos2:]
                    testnewseq = truenewseq[pos1-5: pos2+10]
                    genome_dict[chromsome] = truenewseq
                    print("this is a deletion and insertion",mutation_info, chromsome, pos1,pos2,del_num,insseq.group(1),oldseq, testnewseq)
            elif mutation_type.endswith("del"):  
                pos1, pos2 = position.split("_")
                pos1 = int(pos1)  # 25623958
                pos2 = int(pos2)  # 25623958

                oldseq = genome_dict[chromsome][pos1-5: pos2+10]

                del_num = int(pos2)-int(pos1) + 1
                truenewseq = genome_dict[chromsome][:pos1-1] + del_num*"N" + genome_dict[chromsome][pos2:]
                testnewseq = truenewseq[pos1-5: pos2+10]
                genome_dict[chromsome] = truenewseq
                print("this is in deletion loop",mutation_info, chromsome, pos1, pos2, del_num, oldseq, testnewseq)
            elif "ins" in mutation_type:  # chr9:g.139564754_139564755insC
                pos1, pos2 = position.split("_")
                pos1 = int(pos1)
                pos2 = int(pos2)
                
                oldseq = genome_dict[chromsome][pos1-5: pos2+10]

                ins_num = int(pos2) - int(pos1)
                insseq = re.match(r'.*ins(.*)', mutation_type)
                truenewseq = genome_dict[chromsome][:pos1-1] + genome_dict[chromsome][pos1-1:pos1].lower() + genome_dict[chromsome][pos1:]
                testnewseq = truenewseq[pos1-5: pos1+10]
                genome_dict[chromsome] = truenewseq
                print("this is in insertion loop: ",transcript_ref, mutation_info, chromsome, pos1, pos2, mutation_type, insseq.group(1),oldseq, testnewseq)
                f2.write("{}\tinsert\t{}\t{}\t{}\n".format(transcript_ref,CodingDNAchange,insseq.group(1),ins_num))
            elif mutation_type.endswith("dup"):  # chr7:g.103014422_103015786dup
                pos1, pos2 = position.split("_")
                pos1 = int(pos1)
                pos2 = int(pos2)
                
                oldseq = genome_dict[chromsome][pos1-5: pos2+10]

                dup_num = int(pos2) - int(pos1)
                # insseq = re.match(r'.*ins(.*)', mutation_type)
                truenewseq = genome_dict[chromsome][:pos1-1] + genome_dict[chromsome][pos1-1:pos2+1].lower() + genome_dict[chromsome][pos2+1:]
                testnewseq = truenewseq[pos1-5: pos2+10]
                genome_dict[chromsome] = truenewseq
                print("this is in dup loop:",transcript_ref, mutation_info, chromsome, pos1, pos2, mutation_type, oldseq, testnewseq)
                f2.write("{}\tdup\t{}\t-\t{}\n".format(transcript_ref, CodingDNAchange,dup_num))

        else:
            print ("No match!!")

with open("genome_replaced_mutation.fa","w") as output:
    for each in genome_dict:
        # if each in useful_chr_list:
        output.write(">{}\n{}\n".format(each,genome_dict[each]))
        # output.write(">{}\n{}\n".format(each,genome_dict[each].replace("D","")))
