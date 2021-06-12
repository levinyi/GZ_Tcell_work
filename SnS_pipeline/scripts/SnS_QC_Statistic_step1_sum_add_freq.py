import sys
import re
sum_file = sys.argv[1]
freq_file = sys.argv[2]

adict = {}
with open(freq_file, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("pair"):
            continue
        # pair,CDR3_acc, CDR3_freq, perfect_count, perfect_freq = line.split()
        # adict[pair] = perfect_freq
        a = line.split()
        adict[a[0]] = a[4]


# pattern = re.compile(r'.*\.(clonotype.*)\.[ATCG]{8}.*')
# pattern = re.compile(r'.*\.(.*)\.[ATCG]{12}.*')
with open(sum_file, "r") as f:
    print("Molecule\tCDR3aJ.Identity\tCDR3aJ.Coverage%\tCDR3aJ.Mutation\tCDR3aJ.InDel\tconA.Identity\tconA.Coverage%\tconA.Mutation\tconA.InDel\tTRAV.Identity\tTRAV.Coverage%\tTRAV.Mutation\tTRAV.InDel\tCDR3bJ.Identity\tCDR3bJ.Coverage%\tCDR3bJ.Mutation\tCDR3bJ.InDel\tconB.Identity\tconB.Coverage%\tconB.Mutation\tconB.InDel\tTRBV.Identity\tTRBV.Coverage%\tTRBV.Mutation\tTRBV.InDel\tSum.of.5.TRA\tSum.of.5.TRB\tTCR_freq")
    for line in f:
        if line.startswith(("Molecule","molecule")):
            continue
        line = line.rstrip("\n")
        c = line.split()
        sum5_TRA = int(c[3])+int(c[4])+int(c[7])+int(c[8])+int(c[11])+int(c[12])+int(c[15])+int(c[16])+int(c[19])+int(c[20])
        sum5_TRB = int(c[3])+int(c[4])+int(c[7])+int(c[8])+int(c[23])+int(c[24])+int(c[15])+int(c[16])+int(c[19])+int(c[20])
        

        m = c[0].split(".")
        # print(m)

        # 20210526 new batch: P0036-SnS-Gen3p3_BkBdWv5_WouterLiAl/20210402_SnS_Gen3p3_wouterOVA5LiliAlena_hsTRBCRQR8/
        # for example molecule name: G382E2L1.194.CD4.1000.1.b|SP08.AAAGAATTTTCG.1
        # split with ".": G382E2L1 194 CD4 1000 1 b|SP08 AAAGAATTTTCG 1
        # pair in acc_freq file:   190.CD4.1000.1|SP05
        #m[-3] = m[-3][1:]
        # print(m)
        #match = "{}{}".format(".".join(m[1:-3]), m[-3])
        
        # 20210611 new batch: P0037-SnS-Gen3p3_NCIWv1_CR13/20210604_SnS_NCIWv1_CR13all_hsTRBCRQR8_N12_ExoI/G412-G419/UMI5_Choose20210609
        # for example molecule name : G412E2L1.NCI4323.Clone1.1.e|SP1.ACAGCGGGGTGG.1
        # split woth "." G412E2L1 NCI4323 Clone1 1 e|SP1 ACAGCGGGGTGG 1
        # pair in acc_freq file: NCI4323.Clone252.1
        #m = c[0].split(".")
        # print(m)
        #match = ".".join(m[1:-3])
         
        # for: G421E1L1.a0_b0|SP2.TGGAGCAACCAC.1
        # split with "." G421E1L1 a0_b0|SP2 TGGAGCAACCAC 1
        # pair in acc_freq file: a1341_b1714|SP2
        match = m[1]
        # print("match :", match)
        # match = ".".join(c[0].split(".")[1:-2])
        tcr_freq = adict[match]

        print("{}\t{}\t{}\t{}".format("\t".join(c),sum5_TRA,sum5_TRB,tcr_freq))
