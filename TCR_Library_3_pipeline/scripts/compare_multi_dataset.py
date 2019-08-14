import sys
import os


def deal_M1_file(afile):
    TRA_list = []
    TRB_list = []
    with open(afile, "r") as f:
        for line in f:
            alpha, AV, Acdr3, AJ, beta, BV, Bcdr3, BJ, a, b, c, d, e, f = line.split("\t")
            TRA_clonotype = AV + Acdr3 + AJ
            TRB_clonotype = BV + Bcdr3 + BJ
            TRA_list.append(TRA_clonotype)
            TRB_list.append(TRB_clonotype)
    return list(set(TRA_list)), list(set(TRB_list))


def deal_readcount_file(afile):
    alist = []
    adict = {}
    with open(afile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            clonotype, count = line.split("\t")
            alist.append(clonotype)
            adict[clonotype] = count
    return list(set(alist)),adict


def main():
    dest_mixcr_M1 = sys.argv[1]
    dest_name = os.path.basename(dest_mixcr_M1).split(".")[0]

    set1_mixcr_M1 = sys.argv[2]
    set1_TRA_readcount_file = sys.argv[3]
    set1_TRB_readcount_file = sys.argv[4]
    column1 = os.path.basename(set1_mixcr_M1).split(".")[0]

    set2_mixcr_M1 = sys.argv[5]
    set2_TRA_readcount_file = sys.argv[6]
    set2_TRB_readcount_file = sys.argv[7]
    column2 = os.path.basename(set2_mixcr_M1).split(".")[0]

    set1_M1_TRA_list, set1_M1_TRB_list = deal_M1_file(set1_mixcr_M1)
    set1_TRA_readcount_list,set1_TRA_readcount_dict = deal_readcount_file(set1_TRA_readcount_file)
    set1_TRB_readcount_list,set1_TRB_readcount_dict = deal_readcount_file(set1_TRB_readcount_file)

    set2_M1_TRA_list, set2_M1_TRB_list = deal_M1_file(set2_mixcr_M1)
    set2_TRA_readcount_list,set2_TRA_readcount_dict = deal_readcount_file(set2_TRA_readcount_file)
    set2_TRB_readcount_list,set2_TRB_readcount_dict = deal_readcount_file(set2_TRB_readcount_file)

    print("%(dest_name)s_pairs_TRA\t%(dest_name)s_pairs_TRB\tsupporting_reads\tM1\talpha_in_%(column1)s_pairs\talpha_in_%(column1)s_total\tbeta_in_%(column1)s_pairs\tbeta_in_%(column1)s_total\talpha_in_%(column2)s_Pairs\talpha_in_%(column2)s_total\tbeta_in_%(column2)s_pairs\tbeta_in_%(column2)s_total\t%(column1)s_A_support_reads\t%(column1)s_B_support_reads\t%(column2)s_A_support_reads\t%(column2)s_B_support_reads" % {'dest_name': dest_name, 'column1': column1, 'column2': column2})
    with open(dest_mixcr_M1, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            alpha, AV, Acdr3, AJ, beta, BV, Bcdr3, BJ, a, b, c, d, e, f = line.split("\t")
            TRA_clonotype = AV + Acdr3 + AJ
            TRB_clonotype = BV + Bcdr3 + BJ
            # print TRB_clonotype
            if TRA_clonotype in set1_M1_TRA_list:
                A1 = "Yes"
            else:
                A1 = "No"
            if TRA_clonotype in set1_TRA_readcount_list:
                B1 = "Yes"
            else:
                B1 = "No"
            if TRB_clonotype in set1_M1_TRB_list:
                C1 = "Yes"
            else:
                C1 = "No"
            if TRB_clonotype in set1_TRB_readcount_list:
                D1 = "Yes"
            else:
                D1 = "No"
##############################################################
            if TRA_clonotype in set2_M1_TRA_list:
                A2 = "Yes"
            else:
                A2 = "No"
            if TRA_clonotype in set2_TRA_readcount_list:
                B2 = "Yes"
            else:
                B2 = "No"
            if TRB_clonotype in set2_M1_TRB_list:
                C2 = "Yes"
            else:
                C2 = "No"
            if TRB_clonotype in set2_TRB_readcount_list:
                D2 = "Yes"
            else:
                D2 = "No"
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (TRA_clonotype, TRB_clonotype, a, f, A1, B1, C1, D1, A2, B2, C2, D2, set1_TRA_readcount_dict.get(TRA_clonotype,"NULL"), set1_TRB_readcount_dict.get(TRB_clonotype,"NULL"), set2_TRA_readcount_dict.get(TRA_clonotype,"NULL"), set2_TRB_readcount_dict.get(TRB_clonotype,"NULL")))

if __name__ == '__main__':
    main()
