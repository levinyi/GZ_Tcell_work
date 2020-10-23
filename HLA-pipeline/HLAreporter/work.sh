./HLAreporter.sh 1MT000037WB1 HLA_DRB1 1MT000037WB1_R1.fastq.gz 1MT000037WB1_R2.fastq.gz

$1=1MT000037WB1_R1.fastq.gz
$2=1MT000037WB1_R2.fastq.gz
$3=1MT000037WB1

bwa aln exon23_high_resolution_multi_ref.fa 1MT000037WB1_R1.fastq.gz > 1MT000037WB1_1_exon23_high_resolution_multi_ref.sai
bwa aln exon23_high_resolution_multi_ref.fa 1MT000037WB1_R2.fastq.gz > 1MT000037WB1_2_exon23_high_resolution_multi_ref.sai
bwa sampe exon23_high_resolution_multi_ref.fa 1MT000037WB1_1_exon23_high_resolution_multi_ref.sai 1MT000037WB1_2_exon23_high_resolution_multi_ref.sai 1MT000037WB1_R1.fastq.gz 1MT000037WB1_R2.fastq.gz > 1MT000037WB1_exon23_high_resolution_multi_ref.sam
samtools view -bS 1MT000037WB1_exon23_high_resolution_multi_ref.sam > 1MT000037WB1_exon23_high_resolution_multi_ref.bam
samtools view -b -F 4 1MT000037WB1_exon23_high_resolution_multi_ref.bam > 1MT000037WB1_exon23_high_resolution_multi_ref_mappedreads.bam

samtools sort 1MT000037WB1_exon23_high_resolution_multi_ref_mappedreads.bam -o 1MT000037WB1_exon23_high_resolution_multi_ref_mappedreads_sorted.bam
samtools index 1MT000037WB1_exon23_high_resolution_multi_ref_mappedreads_sorted.bam
# rm 1MT000037WB1_exon23_high_resolution_multi_ref.sam
# rm 1MT000037WB1_exon23_high_resolution_multi_ref.bam
# rm 1MT000037WB1_exon23_high_resolution_multi_ref_mappedreads.bam
