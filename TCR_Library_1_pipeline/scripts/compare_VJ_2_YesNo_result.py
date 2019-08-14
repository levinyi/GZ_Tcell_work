import sys
import os
import reshape_clonotype

def usage():
    '''
    try me:
    cd /cygene/work/29.G22_G23_G27_G28_G29_G30/compare_26.G24-G25.VS.G22_G23
    python compare_VJ_2_YesNo_result.py G24E1L1.pairs.freq.relaxed.M1.xls G22E2L1.clonotype.umi.reads.xls G22E2L2.clonotype.umi.reads.xls G23E2L1.clonotype.umi.reads.xls G23E2L2.clonotype.umi.reads.xls
    '''
    pass


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
            if not line.startswith("TR"):
                continue
            line = line.rstrip("\n")
            # clonotype, count = line.split("\t")
            a = line.split("\t")  # TRAV26-1CIVRVVDSGSWSGDMRFTRAJ43 TRAV26-1        CIVRVVDSGSWSGDMRF       TRAJ43  2       15
            alist.append(a[0])
            adict[a[0]] = a[5]
    return list(set(alist)), adict


def main():
    dest_mixcr_M1 = sys.argv[1]
    dest_name = os.path.basename(dest_mixcr_M1).split(".")[0]

    set1_TRA_readcount_file = sys.argv[2]
    set1_TRB_readcount_file = sys.argv[3]
    column1 = os.path.basename(set1_TRA_readcount_file).split(".")[0]

    set2_TRA_readcount_file = sys.argv[4]
    set2_TRB_readcount_file = sys.argv[5]
    column2 = os.path.basename(set2_TRA_readcount_file).split(".")[0]

    set1_TRA_readcount_list, set1_TRA_readcount_dict = deal_readcount_file(set1_TRA_readcount_file)
    set1_TRB_readcount_list, set1_TRB_readcount_dict = deal_readcount_file(set1_TRB_readcount_file)

    # set2_M1_TRA_list, set2_M1_TRB_list = deal_M1_file(set2_mixcr_M1)
    set2_TRA_readcount_list, set2_TRA_readcount_dict = deal_readcount_file(set2_TRA_readcount_file)
    set2_TRB_readcount_list, set2_TRB_readcount_dict = deal_readcount_file(set2_TRB_readcount_file)

    print("%(dest_name)s_pairs_TRA\t%(dest_name)s_TRAV\t%(dest_name)s_TRACDR3\t%(dest_name)s_TRAJ\t%(dest_name)s_pairs_TRB\t%(dest_name)s_TRBV\t%(dest_name)s_TRBCDR3\t%(dest_name)s_TRBJ\tsupporting_reads\tM1\talpha_in_%(column1)s_pairs\tbeta_in_%(column1)s_pairs\talpha_in_%(column2)s_Pairs\tbeta_in_%(column2)s_pairs\t%(column1)s_A_support_reads\t%(column1)s_B_support_reads\t%(column2)s_A_support_reads\t%(column2)s_B_support_reads" % {'dest_name': dest_name, 'column1': column1, 'column2': column2})
    with open(dest_mixcr_M1, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            alpha, AV, Acdr3, AJ, beta, BV, Bcdr3, BJ, a, b, c, d, e, f = line.split("\t")
            TRA_clonotype = AV + Acdr3 + AJ
            TRB_clonotype = BV + Bcdr3 + BJ
            if TRA_clonotype in set1_TRA_readcount_list:
                B1 = "Yes"
            else:
                B1 = "No"
            if TRB_clonotype in set1_TRB_readcount_list:
                D1 = "Yes"
            else:
                D1 = "No"
##############################################################
            if TRA_clonotype in set2_TRA_readcount_list:
                B2 = "Yes"
            else:
                B2 = "No"
            if TRB_clonotype in set2_TRB_readcount_list:
                D2 = "Yes"
            else:
                D2 = "No"
            TRAV,TRA_cdr3,TRAJ =reshape_clonotype.split_clonotype(TRA_clonotype)
            TRBV,TRB_cdr3,TRBJ = reshape_clonotype.split_clonotype(TRB_clonotype) 
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (TRA_clonotype, TRAV,TRA_cdr3,TRAJ,TRB_clonotype, TRBV,TRB_cdr3,TRBJ,a, f, B1, D1, B2, D2, set1_TRA_readcount_dict.get(TRA_clonotype, "NULL"), set1_TRB_readcount_dict.get(TRB_clonotype, "NULL"), set2_TRA_readcount_dict.get(TRA_clonotype, "NULL"), set2_TRB_readcount_dict.get(TRB_clonotype, "NULL")))

if __name__ == '__main__':
    main()
