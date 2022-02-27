/cygene/software/biosoftware/OptiType/OptiType-1.3.5/data


razers3 -i 95 -m 1 -dr 0 -o fished_1.bam /cygene/software/biosoftware/OptiType/OptiType-1.3.5/data/hla_reference_dna.fasta GB13_TU_CL1_P7-gDNA_2022112_L1_1.fq.gz
razers3 -i 95 -m 1 -dr 0 -o fished_2.bam /cygene/software/biosoftware/OptiType/OptiType-1.3.5/data/hla_reference_dna.fasta GB13_TU_CL1_P7-gDNA_2022112_L1_2.fq.gz

samtools bam2fq fished_1.bam > sample_1_fished.fastq
samtools bam2fq fished_2.bam > sample_2_fished.fastq
rm fished_1.bam fished_2.bam

python3 /cygene/software/biosoftware/OptiType/OptiType-1.3.5/OptiTypePipeline.py -i sample_1_fished.fastq sample_2_fished.fastq --dna -v -o ./

