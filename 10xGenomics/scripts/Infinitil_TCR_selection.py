import sys

import pandas as pd

sys.path.append("/cygene/script/python_packages")
from MultiDict import three_dim

tcr_stats = pd.read_csv(sys.argv[1])
Treg_cluster_list = str(sys.argv[2]).split(",")  # 2,3,4

header = list(tcr_stats.columns)
# print(header)
unwanted = ['barcode','renamed_clone','original_clone','UMAP_1','UMAP_2',
"SCT_snn_res.0.4","SCT_snn_res.0.8","SCT_snn_res.1.2","SCT_snn_res.1.6","SCT_snn_res.2"]
header = [ele for ele in header if ele not in unwanted]
# print(header)

adict = {}

for index, row in tcr_stats.iterrows():
    for each in header:
        three_dim(adict, row['renamed_clone'], row['barcode'], each, row[each])

# import json
# print(json.dumps(adict, indent=4))

'''
{
    "group1_clonotype211": {
        "SA501001-G517E3-CD4_AAACCTGAGTACTTGC": {
            "v_genes": "TRBV5-1,TRAV13-2",
            "j_genes": "TRBJ2-1,TRAJ53",
            "cdr3s_aa": "CASSWQGASEGEQFF;CAENSGGSNYKLTF",
            "TCR_clonality": 1,
            "YostEx_AUCscore": 0.0,
            "TCRx_sig_AUCscore": 0.0339481555333998,
            "TR_1.0_AUCscore": 0.0369200394866732,
            "NeoTCR4_AUCscore": 0.0350947158524427,
            "NeoTCR8_AUCscore": 0.0427255187214329,
            "CD4": 0.0,
            "FOXP3": 1.0,
            "CD8A": 0.0,
            "CD8B": 0.0,
            "seurat_clusters": 20.0
        }
    }
}
'''
for clone, barcodes in adict.items():
    total_barcode = len(barcodes)

    CB_w_MooduleScores = 0    
    Treg_clusters_CB_count = 0
    Non_Treg_cluster_CB_count = total_barcode - Treg_clusters_CB_count

    foxp3_pos_count = 0

    # FACS
    CD4_FACS_CB_count = 0
    CD8_FACS_CB_count = 0
    
    # GEXP
    cd4_CB_count = 0
    CD4SP_GXP_CB_count = 0
    CD8SP_GXP_CB_count	= 0
    DP_GXP_CB_count	= 0
    DN_GXP_CB_count	= 0
    
    for barcode, v in barcodes.items():
        if v['seurat_clusters'] in Treg_cluster_list:
            Treg_clusters_CB_count +=1
        if v['FOXP3'] >0:
            foxp3_pos_count +=1
        # for FACS
        if "CD4" in barcode:
            CD4_FACS_CB_count += 1
        elif "CD8" in barcode:
            CD8_FACS_CB_count += 1
        # for GXP
        if v['CD4'] >0 and v['CD8A'] + v['CD8B'] >0:
            CD4SP_GXP_CB_count +=1
        if v['CD4']==0 and v['CD8A'] + v['CD8B'] >0 :
            CD8SP_GXP_CB_count +=1
        if v['CD4'] >0 and v['CD8A'] + v['CD8B'] >0:
            DP_GXP_CB_count +=0
        if v['CD4']==0 and v['CD8A'] + v['CD8B'] ==0:
            DN_GXP_CB_count +=0
print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
    clone,total_barcode, CB_w_MooduleScores,
    Treg_clusters_CB_count,Treg_clusters_CB_count/total_barcode,
    Non_Treg_cluster_CB_count,Non_Treg_cluster_CB_count/total_barcode,
    foxp3_pos_count, foxp3_pos_count/CB_w_MooduleScores,
    CD4_FACS_CB_count,
    CD8_FACS_CB_count,
))