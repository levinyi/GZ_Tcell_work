# step 1 : split rd1 and rd2.
raw_reads=$1;
name=`ls $raw_reads | awk -F '.' '{print$1}' | awk -F '_' '{print$1}'`;
database='/cygene/work/00.test/pipeline/SnS_pipeline/database'
scripts='/cygene/work/00.test/pipeline/SnS_pipeline/scripts'


# analysis 
# step 1 : deal raw read TRB,
cutadapt -a adapter=GCATTGCACTTGTACGTACG -O 10  -o  $name.p1.useless.fq  -r  $name.p2.useless.fq  --info-file=$name.cutadapt.log  $raw_reads > $name.cutadapt.stats;
python /cygene/work/00.test/pipeline/TCR_Library_3_pipeline/scripts/deal_cutadapt_log.py -l $name.cutadapt.log -d ./ ;
mv $name.p1.fq $name.p1.TRB.fq ;

# deal raw read TRA
cutadapt -a adapter=GGCAGGGTCAGGGTTCTGGAT -O 10  -o  $name.p1.useless.fq  -r  $name.p2.useless.fq  --info-file=${name}_p2.cutadapt.log  $name.p2.fq > $name.p2.cutadapt.stats;
python /cygene/work/00.test/pipeline/TCR_Library_3_pipeline/scripts/deal_cutadapt_log.py -l ${name}_p2.cutadapt.log -d ./ ;
mv ${name}_p2.p2.fq $name.p2.TRA.fq ;
mv ${name}_p2.p1.fq $name.p1.zip.fq ;

# prepare for mapping.
seqkit fq2fa $name.p1.TRB.fq >$name.p1.TRB.fa
seqkit fq2fa $name.p2.TRA.fq >$name.p2.TRA.fa
seqkit fq2fa $name.p1.zip.fq >$name.p1.zip.fa

# blastn

blastn -db ${database}/SplitCDR3JOligosA.CDR3.fa  -query $name.p2.TRA.fa  -out $name.TRA.blast.out -outfmt 7 -num_threads 20 -max_hsps 1 -max_target_seqs 1
blastn -db ${database}/SplitCDR3JOligosB.CDR3.fa  -query $name.p1.TRB.fa  -out $name.TRB.blast.out -outfmt 7 -num_threads 20 -max_hsps 1 -max_target_seqs 1

#statistic
python ${scripts}/export_pairs_table_from_blastout.py $name.TRA.blast.out $name.TRB.blast.out ${database}/SplitCDR3JOligos.zip.fa
