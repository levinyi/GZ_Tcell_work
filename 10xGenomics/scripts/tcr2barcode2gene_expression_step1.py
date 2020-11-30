import sys
import gzip

"""
    first, you should have a tcr list file.
    secend, you should have all_contig_annotations.csv which is generated by cellranger vdj 
    third, you should have barcodes.tsv.gz file which is generated by cellranger counts under filtered_feature_bc_matrix or raw_feature_bc_matrix folder.
"""
tcr_list = sys.argv[1]
contig_annotation_file = sys.argv[2]
barcodes_file = sys.argv[3]

clonotype2barcode_dict = {}
with open(contig_annotation_file, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        c = line.split(",")
        clonotype2barcode_dict.setdefault(c[-2], []).append(c[0])

print("barcode,clonotype,CSF,expansion")
total_tcr = 0
total_barcode = 0
not_found = []
barcode2clonotype_dict = {}
tcr_dict = {}
with open(tcr_list, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("TCR_id"):
            continue
        tcr_id, CFS = line.split()
        tcr_dict[tcr_id] = CFS
        tcr_id_fix = tcr_id.split(".")[0]
        total_tcr += 1
        if tcr_id_fix in clonotype2barcode_dict:
            for barcode in set(clonotype2barcode_dict[tcr_id_fix]):
                total_barcode +=1
                barcode2clonotype_dict[barcode] = tcr_id

with gzip.open(barcodes_file, "rb") as f:
    for line in f:
        line = line.decode("utf-8")
        barcode = line.rstrip("\n")
        if barcode in barcode2clonotype_dict:
            tcr_id = barcode2clonotype_dict[barcode]
            print("{},{},{},Known".format(barcode,tcr_id,tcr_dict[tcr_id]))
        else:
            print("{},,,Not_detected".format(barcode))
'''
sys.stderr.write("total tcr list = {}\n".format(total_tcr))
sys.stderr.write("total barcode  = {}\n".format(total_barcode))
sys.stderr.write("{}\n".format("\n".join(not_found)))
'''
