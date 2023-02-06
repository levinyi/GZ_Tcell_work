# call HaplotypeCaller
/cygene/software/biosoftware/gatk/gatk-4.2.6.1/gatk --java-options "-Xmx30G" \
      HaplotypeCaller \
      -R /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta \
      -I Rawsample.hg38.sortedByCoord.bam \
      -L /cygene3/dushiyi/work/tmp-rbc/gene.bed \
      -O test.out.vcf \
      -contamination 0 \
      -G StandardAnnotation -G StandardHCAnnotation  \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      -ERC GVCF

