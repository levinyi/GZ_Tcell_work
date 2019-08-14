import sys
sys.path.append("/cygene/script/python_packages")
import edit_distance

def usage():
    '''
    python merge_umi.py .umi.count.aa.all.xls >xxx.xls
    '''
    pass

def two_dim_dict(thedict, key_a, key_b, value):
    if key_a in thedict:
        if key_b in thedict[key_a]:
            value = thedict[key_a][key_b] + value
            thedict[key_a].update({key_b: value})
        else:
            thedict[key_a].update({key_b: value})
    else:
        thedict.update({key_a: {key_b: value}})
    return thedict


def main():
    # step 01 :read data and store to a dict. like below:
    
    '''
    # umi_dict = {
      'TRAV1-1CAVRVLKGNNNARLMFTRAJ31': {'TATCACGA': 15}, 
        'TRAV1-1CAVRVSGSARQLTFTRAJ22': {'CATGGAGC': 64, 
                                        'AATGGAGC': 1, 
                                        'TATGGAGC': 1, 
                                        'CTTGGAGC': 1
                                        }
      }
    '''
    umi_dict = {}
    with open(sys.argv[1], "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("unique"):
                print("uniqueCloneType\tUMI\treadsNumber")
                continue
            clone, umi, readCount = line.split("\t")
            umi_dict.setdefault(clone, {})[umi] = int(readCount)
    # print umi_dict

    # step 02 : merge umi.
    umi_merge_dict = {}
    for clone in umi_dict:
        # frist of all, sort dict by number.
        sorted_umi_list = sorted(umi_dict[clone].items(), key=lambda item: item[1], reverse=True)
        # sorted by value. example: [('ATGC', 80), ('ATCG', 70), ('TCGA', 60)]
        # print sorted_umi_list
        # index = 1
        while len(sorted_umi_list) > 0:
            # print("start while --------------\n")
            a = sorted_umi_list.pop(0)  # pop with first,default is last. a is a tuple : ('ATCG',80)
            # print a
            a_umi, a_count = a[:]
            # print("a_umi,a_count = a[:] : %s,%s = a[:]" %(a_umi,a_count))
            two_dim_dict(umi_merge_dict, clone, a_umi, a_count)

            # print "now big dict has: ", umi_merge_dict
            # print "list remain : ", sorted_umi_list
            # print("start for loop sorted_umi_list --------------")
            need_remove_list = []
            for each in sorted_umi_list:
                if edit_distance.minEditDist(a_umi, each[0]) <= 2:
                    # print("compare two: %s and %s and <=2 "%(a_umi,each[0]))
                    two_dim_dict(umi_merge_dict, clone, a_umi, each[1])
                    # print "print again :", umi_merge_dict
                    # sorted_umi_list.remove(each)
                    need_remove_list.append(each)
                # print("sorted_umi_list  remain: ",sorted_umi_list)
            # index += 1
            sorted_umi_list = [x for x in sorted_umi_list if x not in need_remove_list]
            # print "list length: ",len(sorted_umi_list)

    for each_clone in umi_merge_dict:
        for umi in umi_merge_dict[each_clone]:
            print("%s\t%s\t%s"%(each_clone,umi,umi_merge_dict[each_clone][umi]))

if __name__ == '__main__':
    main()