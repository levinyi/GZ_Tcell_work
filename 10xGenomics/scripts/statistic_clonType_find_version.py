import sys
import json

# cellranger file

# clonotype file
VDJ_path = sys.argv[1]
GEX_path = sys.argv[2]
analysis_path = sys.argv[3]

version_file = VDJ_path+"/_versions"
GEX_csv_file = GEX_path+"/outs/metrics_summary.csv"
vdj_csv_file = VDJ_path+"/outs/metrics_summary.csv"
analysis_file = analysis_path+"/merged_tcr_table.csv"
with open(version_file, "r") as f:
    versions = json.load(f)
    print("cellranger: {}".format(versions['pipelines']))

# sequence depth
# before-QC total Cell number 
# before-QC total TCRs number
# Post-QC total Cell number
# Post-QC total TCR pair number
import pandas as pd
GEX_data = pd.read_csv(GEX_csv_file)
# print("sequence depth: {}".format(data['Estimated Number of Cells'][0]))
print("before-QC total Cell number GEX: {}".format(GEX_data['Estimated Number of Cells'][0]))
VDJ_data = pd.read_csv(vdj_csv_file)
print("before-QC total Cell number VDJ: {}".format(VDJ_data['Estimated Number of Cells'][0]))

# "raw_clonotype_id" used to find the unique clonotype
# read the analysis file
analysis_data = pd.read_csv(analysis_file)
print("Paired TCR number: {}".format(len(pd.unique(analysis_data['raw_clonotype_id']))))

CD4_file = analysis_path+"/CD4_Cluster_neoTCR.csv"
CD8_file = analysis_path+"/CD8_Cluster_neoTCR.csv"
CD4Treg_file = analysis_path+"/CD4_Treg_Cluster_neoTCR.csv"

cd4_data = pd.read_csv(CD4_file,sep="\t")
# print("cd4_data :", len(cd4_data))
# print("cd4_data :", cd4_data.head())
print("CD4 cell number: {}".format(len(cd4_data['clonotype'])))
cd8_data = pd.read_csv(CD8_file, sep="\t")
print("CD8 cell number: {}".format(len(cd8_data['clonotype'])))
cd4treg_data = pd.read_csv(CD4Treg_file, sep="\t")
print("CD4Treg cell number: {}".format(len(cd4treg_data['clonotype'])))
