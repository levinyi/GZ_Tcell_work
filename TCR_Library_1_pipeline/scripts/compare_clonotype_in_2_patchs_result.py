import sys


def uasge():
    '''
    test by type:
    python compare_clonotype_in_2_patchs_result.py /cygene/work/29.G22_G23_G27_G28_G29_G30/compare_G29_G30_G26/G29E1L1.clonotype.umi.reads.xls /cygene/work/29.G22_G23_G27_G28_G29_G30/compare_G29_G30_G26/G0026E6L1_R1.split.reads.xls
    '''
    pass


class Clonotype:

    def __init__(self, TRV, CDR3, TRJ):
        self.TRV = TRV
        self.TRJ = TRJ
        self.CDR3 = CDR3

    def conbine_name(self):
        print("%s|%s|%s" % (self.TRV, self.CDR3, self.TRJ))

    def grade(self, grade):
        self.grade = grade

class cells(object):
    def __init__(self,clone_list,TRV_list,TRJ_list,CDR3_list):
        self.clonotypes = clone_list
        self.TRV_list = TRV_list
        self.TRJ_list = TRJ_list
        self.CDR3_list = CDR3_list

    def add_to_list(self,TRV):
        TRV_list.append(TRV)

def deal_file(f):
    a_dict = {}
    clone_list = []
    TRV_list = []
    TRJ_list = []
    CDR3_list = []
    # instance_list = []
    with open(f, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("Clonotype"):
                continue
            a = line.split("\t")
            # aninstance = Clonotype(a[1], a[2], a[3])
            a_dict[a[0]] = line
            # instance_list.append(aninstance)
            clone_list.append(a[0])
            TRV_list.append(a[1])
            CDR3_list.append(a[2])
            TRJ_list.append(a[3])
    return a_dict, clone_list, TRV_list,TRJ_list,CDR3_list


def main():
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    f1_dict, f1_clone_list,f1_TRV_list,f1_TRJ_list,f1_CDR3_list = deal_file(file1)
    f2_dict, f2_clone_list,f2_TRV_list,f2_TRJ_list,f2_CDR3_list = deal_file(file2)
    union = list(set(f1_clone_list).union(set(f2_clone_list)))

    for each in union:
        # if each not in file1_dict:
        #     print each.TRV.counter() 
        #     print each.TRJ
        print("%s\t%s" % (f1_dict.get(each, each + "\tNULL\tNULL\tNULL\tNULL\tNULL"), f2_dict.get(each, "\tNULL\tNULL\tNULL\tNULL\tNULL")))

if __name__ == '__main__':
    main()
