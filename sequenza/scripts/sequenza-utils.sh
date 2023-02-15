# step 1.
sequenza-utils gc_wiggle -w 50 --fasta Homo_sapiens_assembly38.fasta -o hg38.gc50Base.wig.gz

# step 2.
sequenza-utils bam2seqz -n Normal.bam -t Tumor.bam --fasta Homo_sapiens_assembly38.fasta -gc hg38.gc50Base.wig.gz -o sample.seqz.gz

# step 3.
sequenza-utils seqz_binning --seqz sample.seqz.gz -w 50 -o small.seqz.gz

# or combine step2 and step3 in one line.
sequenza-utils bam2seqz -gc /reference/GRCh38.gc50Base.txt.gz --fasta /reference/GRCh38.d1.vd1.fa -n /data/normal.bam --tumor /data/tumor.bam -C chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chr24 chrX | sequenza-utils seqz_binning -w 50 -s - | gzip > /results/tumor_small.seqz.gz


