import csv
import gzip
import os
import scipy.io
from itertools import islice 
 
matrix_dir = "./"
# mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx.gz"))


features_path = os.path.join(matrix_dir, "features.tsv.gz")
feature_ids = [row[0] for row in csv.reader(gzip.open(features_path), delimiter="\t")]
gene_names = [row[1] for row in csv.reader(gzip.open(features_path), delimiter="\t")]
feature_types = [row[2] for row in csv.reader(gzip.open(features_path), delimiter="\t")]

barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path), delimiter="\t")]

with gzip.open("matrix.mtx.gz","r") as f:
    for line in islice(f,3,None):
        a,b,c = line.split()
        print("{}\t{}\t{}\t{}\t{}".format(feature_ids[int(a)-1], gene_names[int(a)-1], feature_types[int(a)-1], barcodes[int(b)-1], c))
        