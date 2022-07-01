# prepare index, and put in ./index/INDEX.seq
# e.g., ./index/TXTO_Condensed_NameZipHvLtFull_NoFiller_CoVAbDab_SP1_20220223_235316.seq

# step 1: convert ab1 files to fastq files
01-ab1_to_fastq.pl

# step 2: trim low-quality bases from 3' end
02-qual_trim.pl

# step 3: assembly
03-analysis.pl ./index/INDEX.seq

# step 4: compare assembled sequence with index
04_check_output.pl ./index/INDEX.seq > sum.txt 
