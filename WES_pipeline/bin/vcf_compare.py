import sys

def usage():
    """
    compare two vcf file . ignore line where filter is not PASS.

    python diff_vcf.py vcf1 vcf2
    
    """

def deal_vcf_file(avcf):
    adict = {}
    alist = []
    if avcf.endswith(".gz"):
        import gzip
        handle = gzip.open(avcf, "rb")
        print("yes gzip")
    else:
        handle = open(avcf, "r")

    for line in handle:
        if line.startswith("#"):
            continue
        a = line.split("\t")
        identity = a[0] + '_' + a[1]
        alist.append(identity)
        adict[identity] = a
    return alist, adict



file0 = sys.argv[1]
# file0 = "CR001001.filtered.liftOver.vcf"
# file1 = "/cygene/work/P0000-Blackbird/CR001-GA001_ByNovogene/CR001001_TU_gDNA.merged.xls"
file1 = sys.argv[2]

# output1 = open("shared_but_not_pass.xls", "w")
# output2 = open("shared_and_pass.xls", "w")
# output3 = open("overdetected.xls", "w")
# output4 = open("not_detected.xls", "w")

vcf1_list, vcf1_dict = deal_vcf_file(file0)
vcf2_list, vcf2_dict = deal_vcf_file(file1)

# new_list = list(set(vcf1_list).union(vcf2_list))
new_list = list(set(vcf1_list).intersection(vcf2_list))

for each in new_list:
    print("{}\t{}".format(vcf1_dict[each], vcf2_dict[each]))
    '''
    if each in vcf1_list and each in vcf2_list:
        if vcf1_dict[each][6] == "PASS":
            output2.write("{}\t{}\t{}\n".format(each,vcf1_dict[each][9].split(":")[2], vcf2_dict[each][41].split(":")[4]))
        else:
            output1.write("{}\t{}\t{}\n".format(each,vcf1_dict[each][9].split(":")[2], vcf2_dict[each][41].split(":")[4]))
    if each in vcf1_list and each not in vcf2_list:
        output3.write("{}\t{}\n".format(each,vcf1_dict[each][9].split(":")[2] ))
    if each not in vcf1_list and each in vcf2_list:
        output4.write("{}\t{}\n".format(each,vcf2_dict[each][41].split(":")[4]))
output1.close()
output2.close()
output3.close()
output4.close()
'''