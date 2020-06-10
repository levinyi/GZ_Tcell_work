import numpy as np
import sys,os


def usage():
    print('''Usage:
    python3 {0} <permutation_result> <xls> > output_file

Example:
    python3 6.permutation_true.py \\
            permutation_test_10000.out.txt \\
            Total.G87E2.pairs.unique.3.update.merged.xls  > 6.result.xls
Update:
    20200607: optimized code.
    20200316: use sys to import input file.
    20200312: updated.
    20200309: created.
    '''.format(os.path.basename(sys.argv[0])))


def main():
    if len(sys.argv) != 3:
        usage()
        sys.exit("Error: Wrong input files.")

    # data = np.loadtxt("permutation_test.10000.out.txt")
    data = np.loadtxt(sys.argv[1])
    print("TRA_clone\tTRA_clone_wells\tTRB_clone\tTRB_clone_wells\tshared_wells\tp1\tp2")
    input_file = sys.argv[2]  # Total.G87E2.pairs.unique.4.xls
    with open(input_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("TR"):
                continue
            TRA_cloneId, TRA_wells, TRB_cloneId, TRB_wells, shared_wells, p = line.split("\t")
            bigger = []
            for each in data[:, int(shared_wells) - 1]:
                if float(p) > each:
                    bigger.append(each)
            if len(bigger) == 0:
                p_min = 0
            else:
                p_min = len(bigger) / float(10000)

            print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                TRA_cloneId, TRA_wells, TRB_cloneId,
                TRB_wells, shared_wells, p, p_min)
            )


if __name__ == "__main__":
    main()
