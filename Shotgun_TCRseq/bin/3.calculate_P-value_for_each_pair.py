import sys,os
import scipy.special


def usage():
    print("""Usage:
    python3 {0} <Total clone barcode file> > <output file>

Example:
    time python3 {0} Total.G87E2.TRAB.clone.reads.barcode.txt >Total.G87E2.pairs.unique.3.xls

Updates:
    20200607    change wells as a parameter.(96wells,or 48wells.)
    20200312    optimimized speed at line 77.
    20200312    fix bugs.(calculate all min(b).see: line 87.)
    20200309    created
    """.format(os.path.basename(sys.argv[0])))

class Clonotype(object):
    """docstring for Clonotype"""
    def __init__(self, barcodes, total_counts, wells):
        self.barcodes = barcodes
        self.total_counts = total_counts
        self.wells = wells
        
def combin(a, b):
    return scipy.special.comb(a, b)


def cal_P(w, wa, wb, wab):
    p = (combin(w, wab)*combin(w-wab, wa-wab)*combin(w-wa, wb-wab))/float(combin(w, wa)*combin(w, wb))
    return p


def add2dict(adict, key, value):
    """docstring for add2dict"""
    if key in adict:
        if value not in adict[key]:
            adict.setdefault(key, []).append(value)
    else:
        adict.setdefault(key, []).append(value)
    return adict


def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict


def main():
    """docstring for main"""
    if len(sys.argv) != 2:
        usage()
        sys.exit("Error: Wrong input file")
    input_file = sys.argv[1]  # Total.G87E2.TRAB.clone.reads.barcode.txt
    TRA_dict = {}
    TRB_dict = {}
    count_dict = {}
    with open(input_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            clone_id, read, barcode = line.split()
            if clone_id not in count_dict:
                count_dict[clone_id] = 1
            else:
                count_dict[clone_id] = count_dict[clone_id] + 1
            if barcode[2:3] == 'A':
                add2dict(TRA_dict, clone_id, barcode[3:])
            elif barcode[2:3] == 'B':
                add2dict(TRB_dict, clone_id, barcode[3:])
            else:
                sys.exit("barcode error.")

    clonotype2p_dict = {}
    clonotype2wells_dict = {}
    shared_wells_dict = {}
    for each_A_clone in TRA_dict:
        barcode_of_A = TRA_dict[each_A_clone]
        Wa = len(barcode_of_A)
        clonotype2wells_dict[each_A_clone] = Wa
        for each_B_clone in TRB_dict:
            barcode_of_B = TRB_dict[each_B_clone]
            Wb = len(barcode_of_B)
            clonotype2wells_dict[each_B_clone] = Wb
            Wab = len(list(set(barcode_of_A).intersection(set(barcode_of_B))))
            shared_wells_dict.setdefault(each_A_clone, {})[each_B_clone] = Wab
            W = 48  # W = 96
            if Wab == 0:
                p = 1
            elif W >= Wab > 0:
                p = 0
                for each_Wab in range(Wab, min(Wa, Wb) + 1):
                    p += cal_P(W, Wa, Wb, each_Wab)
                clonotype2p_dict.setdefault(each_A_clone, {})[each_B_clone] = p
            else:
                sys.exit("Wrong Wab!")

    print("TRA_cloneId\tTRA_clone_wells\tTRB_cloneId\tTRB_clone_wells\tShared_Wells\tmin_P-value")
    for each_A_clone in clonotype2p_dict:
        # modified by author at 20200312.
        min_value = min(clonotype2p_dict[each_A_clone].values())
        min_keys = [k for k in clonotype2p_dict[each_A_clone] if clonotype2p_dict[each_A_clone][k] == min_value]

        for each_trb in min_keys:
            print("{}\t{}\t{}\t{}\t{}\t{}".format(
                each_A_clone,
                clonotype2wells_dict[each_A_clone],
                each_trb,
                clonotype2wells_dict[each_trb],
                shared_wells_dict[each_A_clone][each_trb],
                min_value,)
            )


if __name__ == '__main__':
    main()

