ls /cygene/data/191012_M06509_0014_000000000-CP7V6/G*.fq.gz | while read line; do name=`ls $line |awk -F '/' '{print$NF}' |awk -F '_' '{print$1}'`; \
# ls /cygene/data/191012_M06509_0014_000000000-CP7V6/G*.fq.gz | while read line; do name=${$(basename $line)//\./)[0]}; \
dir_name=$(dirname $line); \
echo "mixcr analyze shotgun --align -OsaveOriginalReads=true --species hs --starting-material rna --receptor-type tcr  -r $name.report $line $name.mixcr.out"; \
echo "mkdir $name";\
echo "mixcr exportReadsForClones -s $name.mixcr.out.clna $name/$name."; \
echo "python /cygene/work/00.test/pipeline/TCR_Library_1_pipeline/scripts/clonetype_support_umi.TSO.py $name.mixcr.out.clonotypes.TRX.txt /xxx/data/umi.fastq.gz  $name.umi.count.aa.all.xls"; \
echo "python /cygene/work/00.test/pipeline/TCR_Library_3_pipeline/scripts/merge_umi.py $name.umi.count.aa.all.xls >$name.merged.umi.count.xls" ; \
echo "Rscript /cygene/work/00.test/pipeline/TCR_Library_1_pipeline/scripts/histogram.R *.merged.umi.count.xls"; \
echo "python /cygene/work/00.test/pipeline/TCR_Library_3_pipeline/scripts/reshape.umi.distance.py $name.merged.umi.count.xls > $name.merged.umi.count.reshape.xls"; \
echo "python /cygene/work/00.test/pipeline/TCR_Library_3_pipeline/scripts/reshape2frequency.py $name.merged.umi.count.reshape.xls > $name.raw.freq.xls"; \
done >work.example.sh
