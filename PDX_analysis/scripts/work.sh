# raw reads align to mm10 and collect unmapped reads
/cygene/software/biosoftware/STAR/STAR-2.7.6a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/mm10/STAR-2.7.6a  --readFilesCommand zcat --readFilesIn ./data/R210628P01-G406E1L2_L2_1.fq.gz ./data/R210628P01-G406E1L2_L2_2.fq.gz  --outFileNamePrefix G406E1L2.mm10. --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outSAMattributes All
/cygene/software/biosoftware/STAR/STAR-2.7.6a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/mm10/STAR-2.7.6a  --readFilesCommand zcat --readFilesIn ./data/R210628P01-G406E2L2_L2_1.fq.gz ./data/R210628P01-G406E2L2_L2_2.fq.gz  --outFileNamePrefix G406E2L2.mm10. --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outSAMattributes All
/cygene/software/biosoftware/STAR/STAR-2.7.6a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/mm10/STAR-2.7.6a  --readFilesCommand zcat --readFilesIn ./data/R210628P01-G406E3L2_L2_1.fq.gz ./data/R210628P01-G406E3L2_L2_2.fq.gz  --outFileNamePrefix G406E3L2.mm10. --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outSAMattributes All
/cygene/software/biosoftware/STAR/STAR-2.7.6a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/mm10/STAR-2.7.6a  --readFilesCommand zcat --readFilesIn ./data/R210628P01-G406E4L2_L2_1.fq.gz ./data/R210628P01-G406E4L2_L2_2.fq.gz  --outFileNamePrefix G406E4L2.mm10. --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outSAMattributes All
echo "finish raw reads align to mm10 and collect unmapped reads"
# raw reads align to hg38 only without collect unmapped reads.
/cygene/software/biosoftware/STAR/STAR-2.7.6a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/hg38/STAR_2.7.6a  --readFilesCommand zcat --readFilesIn ./data/R210628P01-G406E1L2_L2_1.fq.gz ./data/R210628P01-G406E1L2_L2_2.fq.gz  --outFileNamePrefix G406E1L2.hg38. --outSAMtype BAM SortedByCoordinate --outSAMattributes All
/cygene/software/biosoftware/STAR/STAR-2.7.6a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/hg38/STAR_2.7.6a  --readFilesCommand zcat --readFilesIn ./data/R210628P01-G406E2L2_L2_1.fq.gz ./data/R210628P01-G406E2L2_L2_2.fq.gz  --outFileNamePrefix G406E2L2.hg38. --outSAMtype BAM SortedByCoordinate --outSAMattributes All
/cygene/software/biosoftware/STAR/STAR-2.7.6a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/hg38/STAR_2.7.6a  --readFilesCommand zcat --readFilesIn ./data/R210628P01-G406E3L2_L2_1.fq.gz ./data/R210628P01-G406E3L2_L2_2.fq.gz  --outFileNamePrefix G406E3L2.hg38. --outSAMtype BAM SortedByCoordinate --outSAMattributes All
/cygene/software/biosoftware/STAR/STAR-2.7.6a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/hg38/STAR_2.7.6a  --readFilesCommand zcat --readFilesIn ./data/R210628P01-G406E4L2_L2_1.fq.gz ./data/R210628P01-G406E4L2_L2_2.fq.gz  --outFileNamePrefix G406E4L2.hg38. --outSAMtype BAM SortedByCoordinate --outSAMattributes All
echo "echo finish raw reads align to  hg38 only without collect unmapped reads."
# unmapped reads align to hg38 and collect unmapped reads.
/cygene/software/biosoftware/STAR/STAR-2.7.6a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/hg38/STAR_2.7.6a   --readFilesIn G406E1L2.mm10.Unmapped.out.mate1 G406E1L2.mm10.Unmapped.out.mate2  --outFileNamePrefix G406E1L2.mm10.unmapped.to.hg38. --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx  --outSAMattributes All
/cygene/software/biosoftware/STAR/STAR-2.7.6a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/hg38/STAR_2.7.6a   --readFilesIn G406E2L2.mm10.Unmapped.out.mate1 G406E2L2.mm10.Unmapped.out.mate2  --outFileNamePrefix G406E2L2.mm10.unmapped.to.hg38. --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outSAMattributes All
/cygene/software/biosoftware/STAR/STAR-2.7.6a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/hg38/STAR_2.7.6a   --readFilesIn G406E3L2.mm10.Unmapped.out.mate1 G406E3L2.mm10.Unmapped.out.mate2  --outFileNamePrefix G406E3L2.mm10.unmapped.to.hg38. --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outSAMattributes All
/cygene/software/biosoftware/STAR/STAR-2.7.6a/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/hg38/STAR_2.7.6a   --readFilesIn G406E4L2.mm10.Unmapped.out.mate1 G406E4L2.mm10.Unmapped.out.mate2  --outFileNamePrefix G406E4L2.mm10.unmapped.to.hg38. --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outSAMattributes All
echo  "finish unmapped reads align to hg38 and collect unmapped reads."

/usr/local/bin/featureCounts -O -T 20 -t exon -g gene_name -a /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/hg38/gencode.v33.annotation.gtf -o G406E1L2.exon.gene_name.counts.txt  G406E1L2.mm10.unmapped.to.hg38.Aligned.sortedByCoord.out.bam
/usr/local/bin/featureCounts -O -T 20 -t exon -g gene_name -a /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/hg38/gencode.v33.annotation.gtf -o G406E2L2.exon.gene_name.counts.txt  G406E2L2.mm10.unmapped.to.hg38.Aligned.sortedByCoord.out.bam
/usr/local/bin/featureCounts -O -T 20 -t exon -g gene_name -a /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/hg38/gencode.v33.annotation.gtf -o G406E3L2.exon.gene_name.counts.txt  G406E3L2.mm10.unmapped.to.hg38.Aligned.sortedByCoord.out.bam
/usr/local/bin/featureCounts -O -T 20 -t exon -g gene_name -a /cygene/work/00.test/pipeline/RNAseq_pipeline/database/star/hg38/gencode.v33.annotation.gtf -o G406E4L2.exon.gene_name.counts.txt  G406E4L2.mm10.unmapped.to.hg38.Aligned.sortedByCoord.out.bam
echo "finish featureCounts"
python3 /cygene/work/00.test/pipeline/RNAseq_pipeline/bin/featureCounts2TPM.py -a G406E1L2.exon.gene_name.counts.txt -o G406E1L2.exon.gene_name.counts.TPM.txt
python3 /cygene/work/00.test/pipeline/RNAseq_pipeline/bin/featureCounts2TPM.py -a G406E2L2.exon.gene_name.counts.txt -o G406E2L2.exon.gene_name.counts.TPM.txt
python3 /cygene/work/00.test/pipeline/RNAseq_pipeline/bin/featureCounts2TPM.py -a G406E3L2.exon.gene_name.counts.txt -o G406E3L2.exon.gene_name.counts.TPM.txt
python3 /cygene/work/00.test/pipeline/RNAseq_pipeline/bin/featureCounts2TPM.py -a G406E4L2.exon.gene_name.counts.txt -o G406E4L2.exon.gene_name.counts.TPM.txt
echo "finish feature Counts 2 TPM"
##################################################### unique mapped bam 
# for EBV
STAR --runThreadN 10 --genomeDir /cygene/database/EBV/EBV_STAR --readFilesCommand zcat --readFilesIn ./data/R210628P01-G406E1L2_L2_1.fq.gz ./data/R210628P01-G406E1L2_L2_2.fq.gz  --outFileNamePrefix G406E1L2.EBV. --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /cygene/database/EBV/EBV_STAR --readFilesCommand zcat --readFilesIn ./data/R210628P01-G406E2L2_L2_1.fq.gz ./data/R210628P01-G406E2L2_L2_2.fq.gz  --outFileNamePrefix G406E2L2.EBV. --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /cygene/database/EBV/EBV_STAR --readFilesCommand zcat --readFilesIn ./data/R210628P01-G406E3L2_L2_1.fq.gz ./data/R210628P01-G406E3L2_L2_2.fq.gz  --outFileNamePrefix G406E3L2.EBV. --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /cygene/database/EBV/EBV_STAR --readFilesCommand zcat --readFilesIn ./data/R210628P01-G406E4L2_L2_1.fq.gz ./data/R210628P01-G406E4L2_L2_2.fq.gz  --outFileNamePrefix G406E4L2.EBV. --outSAMtype BAM SortedByCoordinate
#
echo "finish EBV"
# for ribosomal
/cygene/software/biosoftware/bbmap/bbduk.sh in1=G406E1L2.mm10.Unmapped.out.mate1 in2=G406E1L2.mm10.Unmapped.out.mate2 outu=G406E1L2.nonribo.fa outm=G406E1L2.ribo.fa ref=/cygene/database/human_ribosomal_DNA/human_ribosomal.fa 2>G406E1L2.bbduk.out
/cygene/software/biosoftware/bbmap/bbduk.sh in1=G406E2L2.mm10.Unmapped.out.mate1 in2=G406E2L2.mm10.Unmapped.out.mate2 outu=G406E2L2.nonribo.fa outm=G406E2L2.ribo.fa ref=/cygene/database/human_ribosomal_DNA/human_ribosomal.fa 2>G406E2L2.bbduk.out
/cygene/software/biosoftware/bbmap/bbduk.sh in1=G406E3L2.mm10.Unmapped.out.mate1 in2=G406E3L2.mm10.Unmapped.out.mate2 outu=G406E3L2.nonribo.fa outm=G406E3L2.ribo.fa ref=/cygene/database/human_ribosomal_DNA/human_ribosomal.fa 2>G406E3L2.bbduk.out
/cygene/software/biosoftware/bbmap/bbduk.sh in1=G406E4L2.mm10.Unmapped.out.mate1 in2=G406E4L2.mm10.Unmapped.out.mate2 outu=G406E4L2.nonribo.fa outm=G406E4L2.ribo.fa ref=/cygene/database/human_ribosomal_DNA/human_ribosomal.fa 2>G406E4L2.bbduk.out
echo "finish Ribosomal"


# pdx mapping summary
python3 ./PDX_analysis/scripts/pdx_summary.py  >pdx_summary.xls
