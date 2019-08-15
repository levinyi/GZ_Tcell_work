#!/usr/bin/bash;
name="$1";
echo "blastn -db /cygene/database/SplitCDR3JOligos/plan2/SplitCDR3JOligos.A.B.fa  -query /cygene/work/40.G34.E2.E2.SnS/data/${name}_S1_L001_R1_001.fa -out ${name}.blast.out -outfmt 7 -num_threads 20 -max_hsps 3 -max_target_seqs 3
echo \"finished step 1\"

python /cygene/work/40.G34.E2.E2.SnS/scripts/deal_pairs.blastout.py ${name}.blast.out  ${name}.blast.out.perfect.pairs.txt ${name}.blast.out.mispairing.txt
echo \"finished step 2\"

python /cygene/work/40.G34.E2.E2.SnS/scripts/pairs.support.reads.py ${name}.blast.out.perfect.pairs.txt /cygene/work/40.G34.E2.E2.SnS/analysis/plan2/SplitCDR3JOligosA.reference.txt  >${name}.each.pairs.reads.xls
echo \"finished step 3\"

less ${name}.each.pairs.reads.xls |awk '\$2>0{print\$1\"\\t\"\$1\"\\t\"\$2}' >${name}.each.pairs.freq
echo \"finished step 4\"

less ${name}.blast.out.mispairing.txt |sort |uniq -c | awk '{print\$2\"\\t\"\$3\"\\t\"\$1}' > ${name}.blast.out.mispairing.freq
echo \"finished step 5\"

cat ${name}.each.pairs.freq ${name}.blast.out.mispairing.freq > ${name}.total.pairs.freq
echo \"finished step 6\"

Rscript /cygene/work/00.test/pipeline/TCR_Library_3_pipeline/scripts/matrix_raster.pairing.R ./ ${name}.total.pairs.filter.freq
echo \"finished step 7\"
" >work_SnS_dsy.sh

echo "Finished Create work_SnS_dsy.sh"