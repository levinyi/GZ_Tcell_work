#!--coding:utf-8
import scipy.special
import random
import sys,os
from itertools import islice
from multiprocessing import Pool


def usage():
    print("""Usage:
    python3 {0} <bools matrix file>  <output file>

Example:
    python3  {0} Total.G87E2.TRAB.96wells.matrix.csv permutation_test.10000.out.txt

Updates:
    20200607: optimized code.
    20200316: updated.
    """.format(os.path.basename(sys.argv[0])))


def combin(a, b):
    return scipy.special.comb(a, b)


def cal_P(w, wa, wb, wab):
    p = (combin(w, wab)*combin(w-wab, wa-wab)*combin(w-wa, wb-wab))/float(combin(w, wa)*combin(w, wb))
    return p


def deal_raw_data(afile):
    TRA_level_dict = {}
    TRB_level_dict = {}
    TRA_level_list = []
    TRB_level_list = []

    with open(afile, "r") as f:
        for line in islice(f, 1, None):
            line = line.rstrip("\n")
            # line = line.rstrip(",")
            clone_types = line.split(",")[1]
            level = sum([int(each) for each in line.split(",")[2:]])
            if clone_types == 'TRA':
                TRA_level_list.append(int(level))
            elif clone_types == 'TRB':
                TRB_level_list.append(int(level))
            else:
                sys.exit("clone type don't match 'TRA/TRB'!")
    column_length = len(line.split(","))-2
    # print(column_length)
    # fix by shiyi.
    for each in range(1, column_length+1):
        TRA_level_dict[each] = TRA_level_list.count(each)
        TRB_level_dict[each] = TRB_level_list.count(each)
    return TRA_level_dict, TRB_level_dict, column_length


def mycallback(delta_list):
    with open(sys.argv[2], 'a+') as f:
        f.writelines("{}\n".format("\t".join(map(str, delta_list))))


def main_work(w, TRB_level_dict):
    delta_list = []  # a empty list store values.
    for j in range(1, w+1):
        p_list = []
        for i in range(1, w+1):
            Ti = TRB_level_dict[i]  # observed value 
            if Ti != 0:
                gamma_i_list = [random.random() for each in range(0, Ti)]
                gamma_i = max(gamma_i_list)
                cumulate_p = 0
                for wab in range(0, min(i, j)+1):
                    if wab >= i+j-w:
                        p = cal_P(w, i, j, wab)
                        if cumulate_p < gamma_i:
                            cumulate_p += p
                            Nij = wab
                        else:
                            Nij = wab - 1
                            break
                p_true = 0
                for each_Nij in range(Nij, w+1):
                    p_true += cal_P(w, i, j, each_Nij)     
                # print("j={},i={},Ti={},max(gamma_i)={},Nij({}-{}),Nij={},cumulate_p={},P={}".format(j, i, Ti, gamma_i, max(0,(i+j-96)), min(i,j), Nij,cumulate_p, p_true))
                p_list.append(p_true)
        pj = min(p_list)
        delta_list.append(pj)
    return delta_list
    # print("{}".format("\t".join(map(str, delta_list))))


def main():
    """docstring for main"""
    if len(sys.argv) != 3:
        usage()
        sys.exit("Error: Wrong input file")

    input_file = sys.argv[1]
    TRA_level_dict, TRB_level_dict, w = deal_raw_data(input_file)
    # w = 96
    # main_work(w, TRB_level_dict)
    
    pool = Pool(w)
    for i in range(1, 10001):
        pool.apply_async(main_work, args=(w, TRB_level_dict,),callback=mycallback)
    pool.close()
    pool.join()


if __name__ == '__main__':
    main()

