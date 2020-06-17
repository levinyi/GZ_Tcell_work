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

def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a:{key_b: val}})
    return thedict

big_dict = {}
gene_list = []
barcode_list = []
with gzip.open("matrix.mtx.gz","r") as f:
    for line in islice(f,3,None):
        a,b,c = line.split()
        addtwodimdict(big_dict, feature_ids[int(a)-1], barcodes[int(b)-1], c)
        gene_list.append(feature_ids[int(a)-1])
        barcode_list.append(barcodes[int(b)-1])
        # print("{}\t{}\t{}\t{}\t{}".format(feature_ids[int(a)-1], gene_names[int(a)-1], feature_types[int(a)-1], barcodes[int(b)-1], c))


print "\t","\t".join(barcode_list)
for i in xrange(len(gene_list)):
    print gene_list[i],
    for j in xrange(len(barcode_list)):
        print big_dict[gene_list[i]].get(barcode_list[j],0),"\t",
    print ""
