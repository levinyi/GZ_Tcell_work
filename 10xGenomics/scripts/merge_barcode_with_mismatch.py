import sys
sys.path.append("/cygene/script/python_packages")
import edit_distance


def main():
    barcode_dict = {}
    with open(sys.argv[1], "r") as f:
        for line in f:
            line = line.rstrip("\n")
            a, b = line.split("\t")
            barcode_dict[a] = int(b)

    merged_dict = {}
    sorted_barcode_list = sorted(barcode_dict.items(), key=lambda item: item[1], reverse=True)

    while len(sorted_barcode_list) > 0:
        a = sorted_barcode_list.pop(0)
        barcode, reads = a[:]
        merged_dict[barcode] = reads

        need_remove_list = []
        for each in sorted_barcode_list:
            if edit_distance.minEditDist(barcode, each[0]) == 1:
                merged_dict[barcode] = merged_dict[barcode] + each[1]
                need_remove_list.append(each)
        sorted_barcode_list = [x for x in sorted_barcode_list if x not in need_remove_list] 


    with open("barcode.merged.txt", "w") as f2:
        for barcode in merged_dict:
            f2.write("{}\t{}\n".format(barcode, merged_dict[barcode]))
        # f.write("Total_Reads\t{}\n".format(total_reads))


if __name__ == '__main__':
    main()
