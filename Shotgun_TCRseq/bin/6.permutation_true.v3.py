import numpy as np
import sys,os


def usage():
    print('''Usage:
    python3 {0} <permutation file> <total pairs> <5000> > <output file>
Example:
    python3 {0} \\
            permutation_test_10000.out.txt \\
            Total.pairs.FromA2B.0317.threshold.50.xls 5000 > Total.pairs.FromA2B.0317.threshold.50.add.xls
Updates:
    20200317: filtered data: p2<0.1 && shared_wells >= 5.
    20200317: update input file format. more 8 columns added.
    20200316: use sys to import input file.
    20200312: updated.
    20200309: created.
    '''.format(os.path.basename(sys.argv[0])))


def main():
    # data = np.loadtxt("permutation_test.10000.out.txt")
    data = np.loadtxt(sys.argv[1])
    print("TRA_clone\tTRA_count\tTRA_clone_wells\tTRB_clone\tTRB_count\tTRB_clone_wells\tshared_wells\tp1\tratio\tp2\tCell.Num.Freq\tTheoretical.Well.Number\tTRA-dropout\tTRB-dropout\tTheoretical.Cell.Num.Freq")
    input_file = sys.argv[2]  # Total.G87E2.pairs.unique.4.xls
    cell_number = eval(sys.argv[3])  # 5000
    with open(input_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("TR"):
                continue
            TRA_cloneId, acount, TRA_wells, TRB_cloneId, bcount, TRB_wells, shared_wells, p, ratio = line.split("\t")
            bigger = []
            for each in data[:, int(shared_wells) - 1]:
                if float(p) > each:
                    bigger.append(each)
            if len(bigger) == 0:
                p_min = 0
            else:
                p_min = len(bigger) / float(10000)
            cell_num_Freq = eval(shared_wells) / float(cell_number)
            Theore_WellNum = eval(TRA_wells) + eval(TRB_wells) - eval(shared_wells)
            TRA_dropout = (Theore_WellNum - int(TRA_wells)) / float(Theore_WellNum)
            TRB_dropout = (Theore_WellNum - int(TRB_wells)) / float(Theore_WellNum)
            Theore_WellNum_freq = Theore_WellNum / float(cell_number)
            
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{}\t{:.2f}\t{:.2f}\t{:.2f}".format(
                TRA_cloneId, acount, TRA_wells, TRB_cloneId,bcount,
                TRB_wells, shared_wells, p, ratio, p_min,
                cell_num_Freq, Theore_WellNum, TRA_dropout, TRB_dropout,
                Theore_WellNum_freq,)
            )


if __name__ == "__main__":
    main()
