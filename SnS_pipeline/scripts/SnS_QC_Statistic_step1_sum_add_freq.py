import sys
import re
import pandas as pd
sum_file = sys.argv[1]
freq_file = sys.argv[2]
out_file = sys.argv[3]

adict = {}
with open(freq_file, "r") as f:
    """pair    CDR3_acc  CDR3_freq  match_AB6_count  match_AB6_freq"""
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("pair"): # skip header
            continue
        # pair,CDR3_acc, CDR3_freq, perfect_count, perfect_freq = line.split()
        # adict[pair] = perfect_freq
        a = line.split()
        # adict[a[0]] = a[4] # match_AB6_freq
        if len(a) == 5:
            adict[a[0]] = a[4] # match_AB6_freq
        else:
            adict[a[0]] = "NA"

# pattern = re.compile(r'.*\.(clonotype.*)\.[ATCG]{8}.*')
# pattern = re.compile(r'.*\.(.*)\.[ATCG]{12}.*')


data = pd.read_csv(sum_file, sep="\t", dtype=str)
header = data.columns.tolist()
# print("old header: ", header)
header_list = ["Sum.of.5.TRA","Sum.of.5.TRB","TCR_freq"]
header.extend(header_list)
# print("new header:", header)

for index, row in data.iterrows():
    sum5_TRA = int(row['CDR3aJ.Mutation'])+int(row['CDR3aJ.InDel'])+int(row['conA.Mutation'])+int(row['conA.InDel'])+int(row['TRAV.Mutation'])+int(row['TRAV.InDel'])+int(row['CDR3bJ.Mutation'])+int(row['CDR3bJ.InDel'])+int(row['conB.Mutation'])+int(row['conB.InDel'])
    sum5_TRB = int(row['CDR3aJ.Mutation'])+int(row['CDR3aJ.InDel'])+int(row['conA.Mutation'])+int(row['conA.InDel'])+int(row['TRBV.Mutation'])+int(row['TRBV.InDel'])+int(row['CDR3bJ.Mutation'])+int(row['CDR3bJ.InDel'])+int(row['conB.Mutation'])+int(row['conB.InDel'])
    
    # clonotype11123.1.e|SP54
    # clonotype11123.1|SP54
    # G0413E1L2.NCI4323.Clone1000.1.b|SP2.ACCCTTTGACAG.1
    # G454E1L1.clonotype1001|SP5.AACAAGGAGGCT.1
    # G456E1L1.RobPrinsCollab.clonotype11123.1.e|SP54.AACACGACTGGG.1 : RobPrinsCollab.clonotype11123.1|SP54
    # G382E2L1.194.CD4.1000.1.b|SP08.AAAGAATTTTCG.1 : 190.CD4.1000.1|SP05
    # G412E2L1.NCI4323.Clone1.1.e|SP1.ACAGCGGGGTGG.1: NCI4323.Clone252.1
    m = row['Molecule'].split(".")
    
    # try:
    m[-3] = m[-3][1:] # G456E1L1.RobPrinsCollab.clonotype11123.1.e|SP54.AACACGACTGGG.1 : RobPrinsCollab.clonotype11123.1|SP54
    tcr_id = "{}{}".format(".".join(m[1:-3]), m[-3])

    # or 
    # tcr_id = ".".join(m[1:-2]) # G454E1L1.clonotype942|SP5.ACCCTCGCTGAC.1
    
    if tcr_id in adict:
        tcr_freq = adict[tcr_id]
    else:
        print("Error: ", tcr_id)
        sys.exit(1)
    data.loc[index, 'Sum.of.5.TRA'] = sum5_TRA
    data.loc[index, 'Sum.of.5.TRB'] = sum5_TRB
    data.loc[index, 'TCR_freq'] = tcr_freq
# data = data.round(decimals=0).astype(object)
# data = data.astype(str)
data['Sum.of.5.TRA'] = data['Sum.of.5.TRA'].astype(int)
data['Sum.of.5.TRB'] = data['Sum.of.5.TRB'].astype(int)
data.to_csv(out_file, sep="\t", index=False, header=header)
    
'''
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
        print(m)
        # G454E1L1.clonotype1001|SP5.AACAAGGAGGCT.1

        # 20210906 new batch: P0037-SnS-Gen3p3_NCIWv1_CR13/20210902_SnS_NCIWv1_mmTRBC_N12_ExoI_PE300_reSeq/G0412-G0413/UMI_10choose
        # for example molecule name : G0413E1L2.NCI4323.Clone1000.1.b|SP2.ACCCTTTGACAG.1
        # 20210526 new batch: P0036-SnS-Gen3p3_BkBdWv5_WouterLiAl/20210402_SnS_Gen3p3_wouterOVA5LiliAlena_hsTRBCRQR8/
        # for example molecule name: G382E2L1.194.CD4.1000.1.b|SP08.AAAGAATTTTCG.1
        # split with ".": G382E2L1 194 CD4 1000 1 b|SP08 AAAGAATTTTCG 1
        # pair in acc_freq file:   190.CD4.1000.1|SP05
        m[-3] = m[-3][1:] # bug fix: remove the first character "b" in the molecule name
        # print(m)
        match = "{}{}".format(".".join(m[1:-3]), m[-3])
        # print(match)

        # 20210611 new batch: P0037-SnS-Gen3p3_NCIWv1_CR13/20210604_SnS_NCIWv1_CR13all_hsTRBCRQR8_N12_ExoI/G412-G419/UMI5_Choose20210609
        # for example molecule name : G412E2L1.NCI4323.Clone1.1.e|SP1.ACAGCGGGGTGG.1
        # split with "." G412E2L1 NCI4323 Clone1 1 e|SP1 ACAGCGGGGTGG 1
        # pair in acc_freq file: NCI4323.Clone252.1
        # print(m)
        # match = ".".join(m[1:-3])
         
        # for example molecule name: G421E1L1.a0_b0|SP2.TGGAGCAACCAC.1
        # split with "." G421E1L1 a0_b0|SP2 TGGAGCAACCAC 1
        # pair in acc_freq file: a1341_b1714|SP2
        # match = m[1]
        # print("match :", match)
        # match = ".".join(c[0].split(".")[1:-2])
        tcr_freq = adict[match]

        print("{}\t{}\t{}\t{}".format("\t".join(c), sum5_TRA, sum5_TRB, tcr_freq))
'''
