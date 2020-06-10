import sys
import scipy.special


def usage():
    print("""Usage
    python {0} <reads barcode file>  <ratio> > <output file>
Example:
    time python3 {0} Total.G87E2.TRAB.clone.reads.barcode.txt 10 >Total.G87E2.pairs.unique.3.xls
Updates:
    20200312    fix bugs.(calculate all min(b).see: line 87.)
    20200309    created
    """.format(os.path.basename(sys.argv[0])))


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


def threetwodimdict(thedict, key_a, key_b, key_c, value):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        if key_b in thedict[key_a]:
            thedict[key_a][key_b].update({key_c: value})
        else:
            thedict[key_a].update({key_b: {key_c: value}})
    else:
        thedict.update({key_a:{key_b: {key_c: value}}})
    return thedict


def calculate_p_value(w):
    prob_dict = {}
    for i in range(1, w+1):
        for j in range(1, w+1):
            for wab in range(max(0, i + j - w), min(i, j) + 1):
                prop = cal_P(w, i, j, wab)
                threetwodimdict(prob_dict, i, j, wab, prop)
    p_value_dict = {}
    for i in range(1, w+1):
        for j in range(1, w+1):
            for wab in range(max(0, i + j - w), min(i, j)+1):
                pvalue = 0
                for each_wab in range(wab, min(i, j)+1):
                    pvalue += prob_dict[i][j][each_wab]
                threetwodimdict(p_value_dict,i,j,wab,pvalue)
    return p_value_dict


def main():
    """docstring for main"""
    if len(sys.argv) != 3:
        usage()
        sys.exit("Error: Wrong input files")

    input_file = sys.argv[1]  # Total.G87E2.TRAB.clone.reads.barcode.txt
    threshold = eval(sys.argv[2]) # should be an integer.
    W = 48  # may be modified every time.
    p_value_dict = calculate_p_value(W)

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
            if W >= Wab > 0:
                p = p_value_dict[Wa][Wb][Wab]
                clonotype2p_dict.setdefault(each_A_clone, {})[each_B_clone] = p

    print("TRA_cloneId\tTRA_count\tTRA_clone_wells\tTRB_cloneId\tTRB_count\tTRB_clone_wells\tShared_Wells\tmin_P-value\tratio")
    for each_A_clone in clonotype2p_dict:
        # modified by author at 20200317.
        candidate_keys = [] 
        for each in sorted(clonotype2p_dict[each_A_clone].items(),key=lambda item : item[1], reverse=False):
            ratio = count_dict[each_A_clone] / float(count_dict[each[0]])
            if float(1/threshold) < ratio < threshold:
                if len(candidate_keys) != 0:
                    if each[1] == candidate_keys[0][1]:
                        candidate_keys.append(each)
                    else:
                        break
                else:
                    candidate_keys.append(each)
        # modified by author at 20200312.
        # min_value = min(clonotype2p_dict[each_A_clone].values())
        # min_keys = [k for k in clonotype2p_dict[each_A_clone] if clonotype2p_dict[each_A_clone][k] == min_value]

        # for each_trb in min_keys:
        for each_can in candidate_keys:
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                each_A_clone,
                count_dict[each_A_clone],
                clonotype2wells_dict[each_A_clone],
                each_can[0],
                count_dict[each_can[0]],
                clonotype2wells_dict[each_can[0]],
                shared_wells_dict[each_A_clone][each_can[0]],
                each_can[1],
                count_dict[each_A_clone]/float(count_dict[each_can[0]]),)
            )


if __name__ == '__main__':
    main()
