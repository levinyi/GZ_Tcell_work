import sys
import numpy as np


alist = []
with open(sys.argv[1],"r") as f:
    for line in f:
        line = line.rstrip("\n")
        alist.append(int(line))

print("mean: {}\nmedian: {}\nmax:{}\nmin: {}".format(np.mean(alist),np.median(alist), np.max(alist),np.min(alist)))

