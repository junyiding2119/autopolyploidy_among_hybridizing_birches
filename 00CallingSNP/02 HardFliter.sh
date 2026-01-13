#!/bin/bash
reference=/home/nwlab/experiment_data/RADseq/reference/35080Betula_pendula.fa

ls GenotypeGVCFs/*gz > chr.list
gatk --java-options "-Xmx6G" MergeVcfs -I chr.list -O 02all.vcf.gz

#include SNP sites
gatk --java-options "-Xmx10g -Xms10g" SelectVariants -R $reference --select-type-to-include SNP --restrict-alleles-to ALL --max-nocall-fraction 0.1 -V 02all.vcf.gz -O 03all.allsnp.vcf.gz
#include non varient sites
gatk --java-options "-Xmx10g -Xms10g" SelectVariants -R $reference --select-type-to-include NO_VARIATION --restrict-alleles-to ALL --max-nocall-fraction 0.1 -V 02all.vcf.gz -O 03all.NO_VARIATION.vcf.gz

#filter SNP #DP 8-400
gatk --java-options "-Xmx6g -Xms6g" VariantFiltration --filter-expression "QUAL < 30" --filter-name "lowQUAL" --filter-expression "ExcessHet > 54.69" --filter-name "highExcessHet" --filter-expression "QD < 2.0" --filter-name "lowQD" --filter-expression "FS > 60.0" --filter-name "highFS" --filter-expression "MQ < 40.0" --filter-name "lowMQ" --filter-expression "MQRankSum < -12.5" --filter-name "lowMQRankSum" --filter-expression "ReadPosRankSum < -8.0" --filter-name "lowReadPosRankSum" --filter-expression "SOR > 3.0" --filter-name "highSOR" --genotype-filter-expression "DP < 8" --genotype-filter-name "lowDP" --genotype-filter-expression "DP > 400" --genotype-filter-name "highDP" -R $reference -V 03all.allsnp.vcf.gz -O 04all.allsnp.filter.vcf.gz
#filter non varient sites DP 4-400
gatk --java-options "-Xmx6g -Xms6g" VariantFiltration --genotype-filter-expression "DP < 4" --genotype-filter-name "lowDP" --genotype-filter-expression "DP > 400" --genotype-filter-name "highDP" -R $reference -V 03all.NO_VARIATION.vcf.gz -O 04all.NO_VARIATION.filter.vcf.gz

gatk --java-options "-Xmx8g -Xms8g" SelectVariants -R $reference --exclude-filtered --set-filtered-gt-to-nocall -V 04all.NO_VARIATION.filter.vcf.gz -O 05all.NO_VARIATION.QC.vcf.gz
gatk --java-options "-Xmx8g -Xms8g" SelectVariants -R $reference --exclude-filtered --set-filtered-gt-to-nocall -V 04all.allsnp.filter.vcf.gz -O 05all.allsnp.QC.vcf.gz

# fix NO_VARIATION vcf
bcftools view --exclude 'STRLEN(REF) >= 2' -Oz -o 05all.fix.NO_VARIATION.QC.vcf.gz 05all.NO_VARIATION.QC.vcf.gz
gatk IndexFeatureFile -I 05all.fix.NO_VARIATION.QC.vcf.gz

# merge NO_VARIATION and SNP vcf
gatk --java-options "-Xmx6G" MergeVcfs -I 05all.fix.NO_VARIATION.QC.vcf.gz -I 05all.allsnp.QC.vcf.gz -O 06all.snp.NO_VARIATION.QC.vcf.gz

# max nocall fraction 0.1
gatk --java-options "-Xmx8g -Xms8g" SelectVariants -R $reference --max-nocall-fraction 0.1 -V 06all.snp.NO_VARIATION.QC.vcf.gz -O 07all.snp.NO_VARIATION.QC.nocall0.1.vcf.gz

