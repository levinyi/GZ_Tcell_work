import sys
sys.path.append("/cygene/script/python_packages")
import edit_distance

with open(sys.argv[1], "r") as f:
    for line in f:
        line = line.rstrip("\n")
        reads, refa, refb, zipa, zipb, flag, zipr = line.split("\t")
        distance1 = edit_distance.minEditDist(zipa, zipb)
        distance2 = edit_distance.minEditDist(zipa, zipr)
        distance3 = edit_distance.minEditDist(zipb, zipr)
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(reads, refa, refb, zipa, zipb, flag, distance1, zipr, distance2, distance3))

