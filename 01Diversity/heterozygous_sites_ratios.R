#!/bin/bash

vcf="/data_analysis/costataeEvolution/Result/00structure_pca/03.allInd.allChr.bisnp.QC.vcf.gz"
export vcf

header_line=$(zcat $vcf | head -61 | tail -1)
NF=$(echo "$header_line" | awk '{print NF}')
indices=$(seq 10 $NF)
ids=$(echo "$header_line" | cut -f 10-$NF)

process_sample() {
    local col_index=$1
    local sample_id=$2
    
    zgrep -v "##" $vcf | cut -f $col_index | grep -v "^\.\/\.\/\.\/\.:" | cut -f 2 -d ":" > ./heterozygous_sites_ratios/allele_frequency.${sample_id}.txt
}

export -f process_sample

paste <(echo "$indices") <(echo "$ids" | tr '\t' '\n') | parallel --colsep '\t' --progress -j 32 process_sample {1} {2}
