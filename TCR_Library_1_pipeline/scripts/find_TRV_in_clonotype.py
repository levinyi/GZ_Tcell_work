import sys
import re
import os


def usage():
    pass
    '''
    python find_TRV_in_clonotype.py G21E6L2.merged.umi.count.xls TRA >outpt.xls
    python find_TRV_in_clonotype.py G21E6L2.merged.umi.count.xls TRB >outpt.xlsa
    The second parameter: alternative for [TRA] or [TRB].  only two parameters.

    '''


def deal_Vgene(afile, TRV_GENE):
    adict = {}
    c = re.compile(r'(TR[AB]V[0-9]+-?[0-9]?(DV\d)?)(.*)(TR[ABD]J.*)')
    with open(afile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("unique"):
                continue
            a = line.split("\t")
            clonotype = a[0]
            p = c.match(clonotype)
            if p:
                v_gene = p.group(1)
                if v_gene in TRV_GENE:
                    if v_gene not in adict:
                        adict[v_gene] = 1
                    else:
                        adict[v_gene] = adict[v_gene] + 1
                else:
                    # print("not found",v_gene)
                    sys.stderr.write('not find %s\n' % v_gene)
            else:
                sys.stderr.write('one not match! %s\n' % line)
                # print("one not match!",line)
    return adict


def main():
    '''
    this is discription.
    '''
    TRAV_GENE = {"TRAV1-1": "TRAV1-1", "TRAV1-2": "TRAV1-2", "TRAV10": "TRAV10", "TRAV11": "TRAV11_Pseudo", "TRAV12-1": "TRAV12-1", "TRAV12-2": "TRAV12-2",
                 "TRAV12-3": "TRAV12-3", "TRAV13-1": "TRAV13-1", "TRAV13-2": "TRAV13-2", "TRAV14DV4": "TRAV14DV4", "TRAV16": "TRAV16", "TRAV17": "TRAV17", "TRAV18": "TRAV18",
                 "TRAV19": "TRAV19", "TRAV2": "TRAV2", "TRAV20": "TRAV20", "TRAV21": "TRAV21", "TRAV22": "TRAV22", "TRAV23DV6": "TRAV23DV6", "TRAV24": "TRAV24",
                 "TRAV25": "TRAV25", "TRAV26-1": "TRAV26-1", "TRAV26-2": "TRAV26-2", "TRAV27": "TRAV27", "TRAV29DV5": "TRAV29DV5", "TRAV3": "TRAV3", "TRAV30": "TRAV30",
                 "TRAV34": "TRAV34", "TRAV35": "TRAV35", "TRAV36DV7": "TRAV36DV7", "TRAV38-1": "TRAV38-1", "TRAV38-2DV8": "TRAV38-2DV8", "TRAV39": "TRAV39", "TRAV4": "TRAV4",
                 "TRAV40": "TRAV40", "TRAV41": "TRAV41", "TRAV5": "TRAV5", "TRAV6": "TRAV6", "TRAV7": "TRAV7", "TRAV8-1": "TRAV8-1", "TRAV8-2": "TRAV8-2", "TRAV8-3": "TRAV8-3",
                 "TRAV8-4": "TRAV8-4", "TRAV8-6": "TRAV8-6", "TRAV8-7": "TRAV8-7_ORF", "TRAV9-1": "TRAV9-1", "TRAV9-2": "TRAV9-2", "TRAV8-5": "TRAV8-5_Pseudo_NULL", "TRAV28": "TRAV28_Pseudo_NULL",
                 "TRAV33": "TRAV33_Pseudo_NULL", "TRAV31": "TRAV31_Pseudo_NULL"}

    TRBV_GENE = {
        "TRBV1": "TRBV1_Pseudo", "TRBV10-1": "TRBV10-1",        "TRBV10-2": "TRBV10-2",       "TRBV10-3": "TRBV10-3", "TRBV11-1": "TRBV11-1", "TRBV11-2": "TRBV11-2",
        "TRBV11-3": "TRBV11-3",     "TRBV12-1": "TRBV12-1_Pseudo", "TRBV12-2": "TRBV12-2_Pseudo", "TRBV12-3": "TRBV12-3", "TRBV12-4": "TRBV12-4", "TRBV12-5": "TRBV12-5",
        "TRBV13": "TRBV13", "TRBV14": "TRBV14", "TRBV15": "TRBV15",
        "TRBV16": "TRBV16", "TRBV17": "TRBV17_ORF", "TRBV18": "TRBV18", "TRBV19": "TRBV19", "TRBV2": "TRBV2", "TRBV20-1": "TRBV20-1", "TRBV21-1": "TRBV21-1",
        "TRBV23-1": "TRBV23-1_ORF", "TRBV24-1": "TRBV24-1", "TRBV25-1": "TRBV25-1", "TRBV26": "TRBV26_Pseudo_NULL", "TRBV27": "TRBV27", "TRBV28": "TRBV28",
        "TRBV29-1": "TRBV29-1", "TRBV3-1": "TRBV3-1", "TRBV3-2": "TRBV3-2", "TRBV30": "TRBV30", "TRBV4-1": "TRBV4-1", "TRBV4-2": "TRBV4-2", "TRBV4-3": "TRBV4-3",
        "TRBV5-1": "TRBV5-1", "TRBV5-3": "TRBV5-3_ORF", "TRBV5-4": "TRBV5-4", "TRBV5-5": "TRBV5-5", "TRBV5-6": "TRBV5-6", "TRBV5-7": "TRBV5-7_ORF", "TRBV5-8": "TRBV5-8",
        "TRBV6-1": "TRBV6-1", "TRBV6-2": "TRBV6-2", "TRBV6-3": "TRBV6-3", "TRBV6-4": "TRBV6-4", "TRBV6-5": "TRBV6-5", "TRBV6-6": "TRBV6-6", "TRBV6-7": "TRBV6-7_ORF",
        "TRBV6-8": "TRBV6-8", "TRBV6-9": "TRBV6-9", "TRBV7-1": "TRBV7-1_ORF", "TRBV7-2": "TRBV7-2", "TRBV7-3": "TRBV7-3", "TRBV7-4": "TRBV7-4", "TRBV7-6": "TRBV7-6",
        "TRBV7-7": "TRBV7-7", "TRBV7-8": "TRBV7-8", "TRBV7-9": "TRBV7-9", "TRBV9": "TRBV9", "TRBV5-2": "TRBV5-2_Pseudo_NULL", "TRBV7-5": "TRBV7-5_Pseudo_NULL", "TRBV22-1": "TRBV22-1_Pseudo_NULL"}

    input_file = sys.argv[1]  # G21E6L2.merged.umi.count.xls
    symbol = sys.argv[2]  # alternative for [TRA] or [TRB].  only two parameters.

    # check file exists
    if not os.path.exists(input_file):
        sys.exit("Error: %s is not exits, please check your input file!" % input_file)

    if symbol == 'TRA':
        adict = deal_Vgene(input_file, TRAV_GENE)
        for each in TRAV_GENE:
            print("{}\t{}".format(TRAV_GENE[each], adict.get(each, 0)))
    elif symbol == 'TRB':
        adict = deal_Vgene(input_file, TRBV_GENE)
        for each in TRBV_GENE:
            print("{}\t{}".format(TRBV_GENE[each], adict.get(each, 0)))
    else:
        sys.exit("Error: input symbol is not correct! please check your input symbol: alternative for [TRA] or [TRB].")

if __name__ == '__main__':
    main()
