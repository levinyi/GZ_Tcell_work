import sys

import pandas as pd
import numpy as np

sys.path.append("/cygene/script/python_packages")
from MultiDict import three_dim

tcr_stats = pd.read_csv(sys.argv[1])
tcr_stats = tcr_stats.fillna('NA')

Treg_cluster_list = str(sys.argv[2]).split(",")  # 2,3,4

# header
header = list(tcr_stats.columns)
unwanted = ['barcode','renamed_clone','original_clone','UMAP_1','UMAP_2',
            "SCT_snn_res.0.4","SCT_snn_res.0.8","SCT_snn_res.1.2",
            "SCT_snn_res.1.6","SCT_snn_res.2"]
header = [ele for ele in header if ele not in unwanted]

# dict
adict = {}
for index, row in tcr_stats.iterrows():
    for each in header:
        if row[each] != 'NA':
            three_dim(adict, row['renamed_clone'], row['barcode'], each, row[each])

import json
# print(json.dumps(adict, indent=4))

'''
{
    "group1_clonotype211": {
        "SA501001-G517E3-CD4_AAACCTGAGTACTTGC": {
            "v_genes": "TRBV5-1,TRAV13-2",
            "j_genes": "TRBJ2-1,TRAJ53",
            "cdr3s_aa": "CASSWQGASEGEQFF;CAENSGGSNYKLTF",
            "TCR_clonality": 1,
            "YostEx_AUCscore": nan,  # ???
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
# print(json.dumps(adict,indent=4))
# print(type(str(adict['group1_clonotype878']["SA501001-G517E1-CD4_AGGTCATTCGCGTAGC"]['YostEx_AUCscore'])))
print("renamed_clone\tTCR_clonality\t#CBs_w/_ModuleScores\tTreg_clusters_CB_count\t%Treg_clusters_CB\tNon-Treg_clusters_CB_count\t%Non-Treg_clusters_CB\tFOXP3_pos_CB_count\t%FOXP3_pos_CB\tCD4_FACS_CB_count\t%CD4_FACS_CB\tCD8_FACS_CB_count\t%CD8_FACS_CB\tCD4SP_GXP_CB_count\t%CD4SP_CB\tCD8SP_GXP_CB_count\t%CD8SP_CB\tDP_GXP_CB_count\t%DP_CB\tDN_GXP_CB_count\t%DN_CB")
for clone, barcodes in adict.items():
    total_barcode = len(barcodes)

    CBs_w_ModuleScores = 0
    Treg_clusters_CB_count = 0

    foxp3_pos_count = 0

    # FACS
    CD4_FACS_CB_count = 0
    CD8_FACS_CB_count = 0
    
    # GEXP
    CD4SP_GXP_CB_count = 0
    CD8SP_GXP_CB_count	= 0
    DP_GXP_CB_count	= 0
    DN_GXP_CB_count	= 0
    
    for barcode, v in barcodes.items():
        if 'YostEx_AUCscore' in v:
            CBs_w_ModuleScores += 1
            if v['seurat_clusters'] in Treg_cluster_list:
                Treg_clusters_CB_count += 1
            if v['FOXP3'] > 0:
                foxp3_pos_count += 1
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

    Non_Treg_cluster_CB_count = total_barcode - Treg_clusters_CB_count
    
    if CBs_w_ModuleScores == 0:
        print("{}\t{}\t{}\t{}\t0\t{}\t0\t{}\t0\t{}\t{}\t{}\t{}\t{}\t0\t{}\t0\t{}\t0\t{}\t0".format(
            clone, total_barcode, CBs_w_ModuleScores,
            Treg_clusters_CB_count,
            Non_Treg_cluster_CB_count, 
            foxp3_pos_count, 
            CD4_FACS_CB_count, CD4_FACS_CB_count/total_barcode,
            CD8_FACS_CB_count, CD8_FACS_CB_count/total_barcode,
            CD4SP_GXP_CB_count,
            CD8SP_GXP_CB_count,
            DP_GXP_CB_count, 
            DN_GXP_CB_count
        ))
    else:
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
            clone, total_barcode, CBs_w_ModuleScores,
            Treg_clusters_CB_count, Treg_clusters_CB_count/CBs_w_ModuleScores,
            Non_Treg_cluster_CB_count, Non_Treg_cluster_CB_count/CBs_w_ModuleScores,
            foxp3_pos_count, foxp3_pos_count/CBs_w_ModuleScores,
            CD4_FACS_CB_count, CD4_FACS_CB_count/total_barcode,
            CD8_FACS_CB_count, CD8_FACS_CB_count/total_barcode,
            CD4SP_GXP_CB_count, CD4SP_GXP_CB_count/CBs_w_ModuleScores,
            CD8SP_GXP_CB_count, CD8SP_GXP_CB_count/CBs_w_ModuleScores,
            DP_GXP_CB_count, DP_GXP_CB_count/CBs_w_ModuleScores,
            DN_GXP_CB_count, DN_GXP_CB_count/CBs_w_ModuleScores,
        ))
