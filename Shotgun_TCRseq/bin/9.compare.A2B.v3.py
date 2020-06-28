import sys

a2b_file = sys.argv[1]  # 7.result.20200312.filtered.xls
b2a_file = sys.argv[2]  # 7.result.20200312.filtered.xls


    
adict = {}
alist = []
with open(a2b_file, "r") as f1:
    for line in f1:
        line = line.rstrip("\n")
        colum_1 = line.split()[0]
        colum_3 = line.split()[3]
        keya = colum_1 + '_' + colum_3
        alist.append(keya)
        adict[keya] = line

bdict = {}
with open(b2a_file, "r") as f2:
    for line in f2:
        line = line.rstrip("\n")
        colum_1 = line.split()[0]
        colum_3 = line.split()[3]
        keyb = colum_3 + '_' + colum_1
        alist.append(keyb)
        bdict[keyb] = line

alist = list(set(alist))
print("TRA_clone\tTRA_count\tTRA_clone_wells\tTRB_clone\tTRB_count\tTRB_clone_wells\tshared_wells\tp1\tratio\tp2\tCell.Num.Freq\tTheoretical.Well.Number\tTRA-dropout\tTRB-dropout\tTheoretical.Cell.Num.Freq\tTRB_clone\tTRB_count\tTRB_clone_wells\tTRA_clone\tTRA_count\tTRA_clone_wells\tshared_wells\tp1\tratio\tp2\tCell.Num.Freq\tTheoretical.Well.Number\tTRA-dropout\tTRB-dropout\tTheoretical.Cell.Num.Freq")
for each in alist:
    print("{}\t{}".format(
        adict.get(each,"NULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL"),
        bdict.get(each, "NULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL"))
    )


