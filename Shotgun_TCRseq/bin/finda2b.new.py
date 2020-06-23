"""docstring of module docstring"""
import sys


def usage():
    """
    python G115red.xls G114green.xls  
    red file must first.

    """

def deal_file(afile, sorted_by):
    """docstring for main"""
    if sorted_by == 'trb':
        column = 4
    elif sorted_by == 'tra':
        column = 1

    name_list = []
    adict = {}
    with open(afile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("Note"):
                continue
            a = line.split()
            name_list.append(a[column])
            adict.setdefault(a[column], []).append(line)
    return name_list, adict


def main():
    """docstring for main"""
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    sorted_by = sys.argv[3]

    name_list_1, file1_dict = deal_file(file1, sorted_by)
    name_list_2, file2_dict = deal_file(file2, sorted_by)
    name_list = list(set(name_list_1 + name_list_2))
    O1 = open("red." + sorted_by + ".xls", "w")
    O2 = open("green." + sorted_by + ".xls", "w")

    for each in name_list:
        if each in file1_dict and each not in file2_dict:
            for each_a in file1_dict[each]:
                O1.write("1\t{}\n".format(each_a))
        elif each in file1_dict and each in file2_dict:
            a_max = 0
            b_max = 0
            for each_c in file1_dict[each]:
                shared_wells_a = each_c.split("\t")[12]
                # print("AAAA",shared_wells_a)
                if int(shared_wells_a) > a_max:
                    a_max = int(shared_wells_a)
                O1.write("2\t{}\n".format(each_c))
            for each_d in file2_dict[each]:
                shared_wells_d = each_d.split("\t")[12]
                # print("shared_wells _d ",shared_wells_d)
                if int(shared_wells_d) > b_max:
                    b_max = int(shared_wells_d)
                O2.write("2\t{}\n".format(each_d))
            # print("{}\t{}\t{}\t".format(each, a_max, b_max))
        elif each not in file1_dict and each in file2_dict:
            for each_b in file2_dict[each]:
                O2.write("3\t{}\n".format(each_b))
    O1.close()
    O2.close()

if __name__ == '__main__':
    main()
