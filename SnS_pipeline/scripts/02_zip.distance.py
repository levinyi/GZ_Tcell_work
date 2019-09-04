import sys
sys.path.append("/cygene/script/python_packages")
import edit_distance
import revers


def usage():
    ''' python ../scripts/zip.distance.py export_pairs_table.out.xls >G34E3L1.pairs.total.step.2.xls '''

print("Read_Id\tAlpha_Ref\tBeta_Ref\tAlpha_Ref_Zip\tBeta_Ref_Zip\tTag\tDistance_of_Alpha_and_Beta_zip\tReads_zip\tdistance_of_Alpha_and_Reads_zip\tdistance_of_Beta_and_Reads_zip\t")
with open(sys.argv[1], "r") as f:
    for line in f:
        line = line.rstrip("\n")
        reads, refa, refb, zipa, zipb, flag, zipr = line.split("\t")  # zipr means zip of reads.
        zipr = revers.DNA_reverse(revers.DNA_complement(zipr))  # added by dsy on:2019-09-03
        distance1 = edit_distance.minEditDist(zipa, zipb)
        distance2 = edit_distance.minEditDist(zipa, zipr)
        distance3 = edit_distance.minEditDist(zipb, zipr)
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(reads, refa, refb, zipa, zipb, flag, distance1, zipr, distance2, distance3))

