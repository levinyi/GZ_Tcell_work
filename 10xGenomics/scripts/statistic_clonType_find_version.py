import sys
import os
import json


project_path = sys.argv[1]
VDJ_path_list = [each for each in os.listdir(project_path) if each.endswith('VDJ')]
GEX_path_list = [each for each in os.listdir(project_path) if each.endswith('GEX')]
analysis_path = project_path + '/analysis'
# print("project_path:", project_path)
# print("analysis_path:", analysis_path)
# print('VDJ_path_list:', VDJ_path_list)
# print('GEX_path_list:', GEX_path_list)

version_file = project_path + '/' + VDJ_path_list[0] +"/_versions"
GEX_csv_files = [project_path + '/'+ each +"/outs/metrics_summary.csv" for each in GEX_path_list]
vdj_csv_files = [project_path + '/'+ each +"/outs/metrics_summary.csv" for each in VDJ_path_list]
# print("version_file:", version_file)
# print('GEX_csv_files:', GEX_csv_files)
# print('vdj_csv_files:', vdj_csv_files)


with open(version_file, "r") as f:
    versions = json.load(f)
    print("# cellranger: {}".format(versions['pipelines']))

analysis_file = analysis_path + "/merged_tcr_table.csv"
# sequence depth
# before-QC total Cell number 
# before-QC total TCRs number
# Post-QC total Cell number
# Post-QC total TCR pair number
import pandas as pd

output_dict = {}
for gex_csv_file in GEX_csv_files:
    name = gex_csv_file.split('/')[-3].split('_')[0] # /cygene2/pipeline/10X/data/LC23/LC23-LN8-CD8neg_GEX/outs/metrics_summary.csv
    GEX_data = pd.read_csv(gex_csv_file)
    output_dict.setdefault(name, {})['before-QC_Total_Cell_number_GEX'] = GEX_data['Estimated Number of Cells'][0]
for vdj_csv_file in vdj_csv_files:
    name = vdj_csv_file.split('/')[-3].split('_')[0] # /cygene2/pipeline/10X/data/LC23/LC23-LN8-CD8neg_VDJ/outs/metrics_summary.csv
    VDJ_data = pd.read_csv(vdj_csv_file)
    # print("before-QC total Cell number VDJ: {}".format(VDJ_data['Estimated Number of Cells'][0]))
    output_dict.setdefault(name,{})['before-QC_Total_Cell_number_VDJ'] = VDJ_data['Estimated Number of Cells'][0]

# print("output_dict:", output_dict)
print("sample\tbefore-QC_Total_Cell_number_GEX\tbefore-QC_Total_Cell_number_VDJ")
for each in output_dict:
    print("{}\t{}\t{}".format(each, output_dict[each].get("before-QC_Total_Cell_number_GEX","null"), output_dict[each].get("before-QC_Total_Cell_number_VDJ","null")))

# read the analysis file
analysis_data = pd.read_csv(analysis_file)
print("Paired TCR number: {}".format(len(pd.unique(analysis_data['raw_clonotype_id']))))

CD4_file = analysis_path+"/CD4_Cluster_neoTCR.csv"
CD8_file = analysis_path+"/CD8_Cluster_neoTCR.csv"
CD4Treg_file = analysis_path+"/CD4_Treg_Cluster_neoTCR.csv"

if os.path.exists(CD4_file):
    cd4_data = pd.read_csv(CD4_file,sep="\t")
    print("CD4 cell number: {}".format(len(cd4_data['clonotype'])))
    cd8_data = pd.read_csv(CD8_file, sep="\t")
    print("CD8 cell number: {}".format(len(cd8_data['clonotype'])))
    cd4treg_data = pd.read_csv(CD4Treg_file, sep="\t")
    print("CD4Treg cell number: {}".format(len(cd4treg_data['clonotype'])))
