import sys
import json

# cellranger file

# clonotype file

version_file = "/cygene2/pipeline/10X/data/HC29/G480E1_VDJ/_versions"
with open(version_file, "r") as f:
    versions = json.load(f)
    print("cellranger: {}".format(versions['pipelines']))

with open() as f:
    versions = json.load(f)
    print("clonotype: {}".format(versions['clonotype']))

# sequence depth
# before-QC total Cell number 
# before-QC total TCRs number
# Post-QC total Cell number
# Post-QC total TCR pair number
import pandas as pd
csv_file = "/cygene2/pipeline/10X/data/HC29/G480E1_GEX/outs/metrics_summary.csv"
data = pd.read_csv(csv_file)
print("sequence depth: {}".format(data['sequence_depth'][0]))
print("before-QC total Cell number: {}".format(data['before_qc_total_cell_number'][0]))
