import sys

def usage():
    """
    compare  vcf file and a xls file.  ignore line where filter is not PASS.

    python diff_vcf.py vcf1 vcf2
    
    """

def deal_with_xls(afile):
    adict = {}
    alist = []
    with open(afile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("CHROM"):
                continue
            a = line.split("\t")
            identity = a[0]+'_'+a[1]
            alist.append(identity)
            adict[identity] = a
    # print(len(a))
    return alist, adict


def deal_vcf_file(avcf):
    adict = {}
    alist = []
    with open(avcf, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            a = line.split("\t")
            identity = a[0].replace("chr", "") + '_' + a[1]
            alist.append(identity)
            adict[identity] = a
    return alist, adict



file0 = sys.argv[1]
# file0 = "CR001001.filtered.liftOver.vcf"
file1 = "/cygene/work/P0000-Blackbird/CR001-GA001_ByNovogene/CR001001_TU_gDNA.merged.xls"

output1 = open("shared_but_not_pass.xls", "w")
output2 = open("shared_and_pass.xls", "w")
output3 = open("overdetected.xls", "w")
output4 = open("not_detected.xls", "w")

vcf_list, vcf_dict = deal_vcf_file(file0)
xls_list, xls_dict = deal_with_xls(file1)

new_list = list(set(vcf_list).union(xls_list))

for each in new_list:
    if each in vcf_list and each in xls_list:
        if vcf_dict[each][6] == "PASS":
            output2.write("{}\t{}\t{}\n".format(each,vcf_dict[each][9].split(":")[2], xls_dict[each][41].split(":")[4]))
        else:
            output1.write("{}\t{}\t{}\n".format(each,vcf_dict[each][9].split(":")[2], xls_dict[each][41].split(":")[4]))
    if each in vcf_list and each not in xls_list:
        output3.write("{}\t{}\n".format(each,vcf_dict[each][9].split(":")[2] ))
    if each not in vcf_list and each in xls_list:
        output4.write("{}\t{}\n".format(each,xls_dict[each][41].split(":")[4]))
output1.close()
output2.close()
output3.close()
output4.close()
