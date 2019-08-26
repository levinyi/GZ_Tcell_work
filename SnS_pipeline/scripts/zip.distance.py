import sys
sys.path.append("/cygene/script/python_packages")
import edit_distance

with open(sys.argv[1],"r") as f:
    for line in f:
        line = line.rstrip("\n")
        reads,refa,refb,zipa,zipb,flag = line.split("\t")
        distance = edit_distance.minEditDist(zipa,zipb)
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(reads,refa,refb,zipa,zipb,flag,distance))
