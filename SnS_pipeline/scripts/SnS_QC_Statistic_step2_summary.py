import sys
import os
import argparse


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-f1', '--countfile', action='store', dest='count_file',   help="count file, usually named G260E2L1.pair.count.txt")
    parser.add_argument('-f2', '--freqfile',  action='store', dest="frequ_file",   help="frequency file, usually named G260E2L1.pair.acc_freq.txt")
    parser.add_argument('-f3', '--sumfile',   action='store', dest="summa_file",   help="sum file which generated by SnS_QC_Statistic_step1_sum_add_freq.py")
    parser.add_argument('-n', '--total_num', action='store', type=int, dest="total_number", help="total TCR number")
    parser.add_argument('-f', '--freq',      action='store', type=float, dest="frequency",    help="tcr frequency")
    parser.add_argument('-u', '--umi',      action='store', type=int, dest="umi", help="umi number")
    return parser.parse_args()


def deal_count_file(count_file):
    with open(count_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("total"):
                continue
            counts_sum = line.split()
    return counts_sum


def deal_frequency_file(frequ_file, frequency, total_number, counts_sum):
    molecule_number = 0
    tcr_freq_cutoff = 0
    acc90_freqcutoff = 0
    acc80_freqcutoff = 0
    match_perfect_number = 0
    with open(frequ_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("pair"):
                continue
            c = line.split()
            if c[1] != '0':
                molecule_number += 1
            if float(c[4]) >= frequency:
                tcr_freq_cutoff += 1
                match_perfect_number += float(c[4])
            if float(c[1]) > 0.9 and float(c[4]) >= frequency:
                acc90_freqcutoff += 1
            if float(c[1]) > 0.8 and float(c[4]) >= frequency:
                acc80_freqcutoff += 1
    print("total(reads),match_A(reads),match_B(reads),six_total(Mol),match_AB6,percent,Mol_per_TCR,~Reads/Molecule,,TCR by at least 1 molecule,TCRs with freq>{cutoff},TCRs: acc>0.90 & Freq>{cutoff},TCRs: acc>0.80 & Freq>{cutoff},match_perfect: TCR Freq > {cutoff}".format(**{"cutoff":frequency}))
    print("{},{},{},,{}/{},{}/{},{}/{},{}/{},{}".format(",".join(counts_sum),float(counts_sum[4])/float(total_number),float(counts_sum[0])/float(counts_sum[3]),
        molecule_number,total_number, 
        tcr_freq_cutoff,total_number, 
        acc90_freqcutoff,total_number, 
        acc80_freqcutoff,total_number, 
        match_perfect_number))
    print("")
    return match_perfect_number


def main():
    parser = _argparse()
    count_file = parser.count_file
    frequ_file = parser.frequ_file
    summa_file = parser.summa_file
    total_number = parser.total_number
    frequency = parser.frequency
    if parser.umi:
        umi = parser.umi
        print("(1) For all following samples. we used reads/UMI >= {} as the cutoff. I have looked through all reads/UMI distribution histograms. the choice of 10 works well.".format(umi))
        print("(2) PE300 on BA-type libraries can generally cover 56% of the entire TRAV region while TRBV is impossible to be estimated by NGS. So we assign 50% coverage for TRAV and use TRAV error rate to estimate TRBV regions. Theoretical TRAV was calculated by multiplying its errors twice.")
        print("")

    file_name = frequ_file.split(".")[0]
    print(file_name)
    # 
    counts_sum = deal_count_file(count_file)

    # 
    match_perfect_number = deal_frequency_file(frequ_file,frequency,total_number, counts_sum)

    if summa_file:
        mut_in_CDR3aJ = 0
        Indel_in_CDR3aJ = 0
        mut_in_conA = 0
        Indel_in_conA = 0
        mut_in_TRAV = 0
        Indel_in_TRAV = 0
        mut_in_CDR3bJ = 0
        Indel_in_CDR3bJ = 0
        mut_in_conB	= 0
        Indel_in_conB = 0

        total_molecules_in_file = 0
        overall_wo_TRBV = 0
        with open(summa_file, "r") as f:
            for line in f:
                if line.startswith("Molecule"):
                    continue
                line = line.rstrip("\n")
                c = line.split()
                if int(c[3]) == 0:
                    mut_in_CDR3aJ += 1
                if int(c[4]) == 0:
                    Indel_in_CDR3aJ += 1
                if int(c[7]) == 0:
                    mut_in_conA += 1
                if int(c[8]) == 0:
                    Indel_in_conA +=1
                if int(c[11]) == 0:
                    mut_in_TRAV += 1
                if int(c[12]) == 0:
                    Indel_in_TRAV += 1
                if int(c[15]) == 0:
                    mut_in_CDR3bJ += 1
                if int(c[16]) == 0:
                    Indel_in_CDR3bJ += 1
                if int(c[19]) == 0:
                    mut_in_conB += 1
                if int(c[20]) == 0:
                    Indel_in_conB += 1
                if int(c[25]) == 0:
                    overall_wo_TRBV += 1
                total_molecules_in_file += 1
        print("Match_perfect molecules in the file \".sum\": ,{}".format(total_molecules_in_file))
        print("mut in CDR3aJ,Indel in CDR3aJ,mut in conA,Indel in conA,mut in CDR3bJ,Indel in CDR3bJ,mut in conB,Indel in conB,mut in TRAV,Indel in TRAV,Overall w/o TRBV,Error-free% in TRAV,Overall w/o TRBV,Overall Err-free w/ TRBV,Error free % at Freq>{}".format(frequency))
        print("{},{},{},{},{},{},{},{},{},{},{},,Theoretical(TRAV was multiplied twice)".format(
            mut_in_CDR3aJ, Indel_in_CDR3aJ, mut_in_conA, Indel_in_conA, 
            mut_in_CDR3bJ, Indel_in_CDR3bJ, mut_in_conB, Indel_in_conB, 
            mut_in_TRAV, Indel_in_TRAV, overall_wo_TRBV))
        Error_free_in_TRAV = (mut_in_TRAV**2)*(Indel_in_TRAV**2)/(float(total_molecules_in_file))**4
        Theoretical_pct = (mut_in_CDR3aJ*Indel_in_CDR3aJ*mut_in_conA*Indel_in_conA*mut_in_TRAV*Indel_in_TRAV*mut_in_CDR3bJ*Indel_in_CDR3bJ*mut_in_conB*Indel_in_conB)/float(pow(total_molecules_in_file,10))
        print("{:%},{:%},{:%},{:%},{:%},{:%},{:%},{:%},{:%},{:%},{:%},{:%},{:%},{:%},{:%}".format(
            mut_in_CDR3aJ/float(total_molecules_in_file), Indel_in_CDR3aJ/float(total_molecules_in_file), 
            mut_in_conA/float(total_molecules_in_file),   Indel_in_conA/float(total_molecules_in_file), 
            mut_in_CDR3bJ/float(total_molecules_in_file), Indel_in_CDR3bJ/float(total_molecules_in_file), 
            mut_in_conB/float(total_molecules_in_file),   Indel_in_conB/float(total_molecules_in_file),
            mut_in_TRAV/float(total_molecules_in_file),   Indel_in_TRAV/float(total_molecules_in_file), 
            overall_wo_TRBV/float(total_molecules_in_file),
            Error_free_in_TRAV,
            Theoretical_pct,
            Theoretical_pct*(mut_in_TRAV/float(total_molecules_in_file))*(Indel_in_TRAV/float(total_molecules_in_file)),
            match_perfect_number*Theoretical_pct*(mut_in_TRAV/float(total_molecules_in_file))*(Indel_in_TRAV/float(total_molecules_in_file)),
            ))
        print("")
    print("")

if __name__ =='__main__':
    main()
