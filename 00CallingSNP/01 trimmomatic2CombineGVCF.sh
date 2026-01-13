#trimmomatic
trim=/home/nwlab/biosoft/software/Trimmomatic-0.39/trimmomatic-0.39.jar
Dir=/data_analysis/costataeEvolution/
cp ~/biosoft/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa .

for sample in {BJF04,BYS08,BYS20,CBS11,DTX14,DTX56,EDG23,FC02,FHSJ32,HR17,LJA06,LJS01,LJS17,LP01,LP08,LYL27,LYL43,LYL46,LYL58,PQG12,PQG25,QSLC10,CKX18-024,CKX18-006,SNJ18-010,TBS95,TBS98,TZZ60,TZZ78,WLP008,WLP022,XLA01,XLA32,XMX23}
do

fq_file1=~/ncbi_upload/2021AoB_buggsii/fastqFile/${sample}_R1.fq.gz
fq_file2=~/ncbi_upload/2021AoB_buggsii/fastqFile/${sample}_R2.fq.gz
outputR1=$Dir/clean_data/${sample}_R1.1.fq.gz
outputR2=$Dir/clean_data/${sample}_R2.2.fq.gz
upoutputR1=$Dir/clean_data/${sample}_R1unpaired.fq.gz
upoutputR2=$Dir/clean_data/${sample}_R2unpaired.fq.gz

java -classpath $trim org.usadellab.trimmomatic.TrimmomaticPE -phred33 $fq_file1 $fq_file2 $outputR1 $upoutputR1 $outputR2 $upoutputR2 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:40 ILLUMINACLIP:~/biosoft/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10
done

#bwa
Dir=/data_analysis/costataeEvolution/
reference=/home/nwlab/experiment_data/RADseq/reference/35080Betula_pendula.fa
for sample in {CKX18-009,CKX18-011,CKX18-019,CKX19-002,CKX19-006,CKX19-010,CKX19-011,CKX19-013,DBS001,DBS006,DBS010,DBS014,DBS017,DBS020,DCX017,DCX024,DCX030,DCX033,DCX043,DQ007,DQ013,DQ017,DQ021,DQ025,DTX009,DTX14,DTX018,DTX023,DTX026,DTX034,DTX035,DTX038,DTX044,DTX049,DTX56,DTX062,DTX066,EBY007,EBY010,EBY015,EBY022,EDG004,EDG008,EDG23,EDG026,EDG035,EDG039,EDG046,EDL010,EDL011}
do
bam_file=$Dir/alignment/${sample}.so.bam
fq_file1=$Dir/clean_data/${sample}_R1.1.fq.gz
fq_file2=$Dir/clean_data/${sample}_R2.2.fq.gz

bwa mem -t 8 -M -R "@RG\tID:${sample}\tSM:${sample}" $reference $fq_file1 $fq_file2 | samtools sort -@ 1 -o $bam_file 
done

#MarkDuplicates
Dir=/hpcfile/users/92016005/data_analysis/RADseq/costataeEvolution
for sample in {}
do
so_bam_file=$Dir/alignment/${sample}.so.bam
md_bam_file=$Dir/markDuplicate/${sample}.so.md.bam
metrics_file=$Dir/markDuplicate/${sample}.so.md.metrics

gatk --java-options "-Xmx2g -Xms2g" MarkDuplicates --REMOVE_DUPLICATES true -I $so_bam_file -O $md_bam_file -M $metrics_file --ASSUME_SORT_ORDER coordinate

#index bam
samtools index $md_bam_file
done

#calling SNP
#tetra
Dir=/data_analysis/costataeEvolution
reference=/home/nwlab/experiment_data/RADseq/reference/35080Betula_pendula.fa

for sample in {BX007,BX018,BX021,BX023,BX026,BX029,BX039,BYS003,BYS015,BYS20,CKX19-006,CKX19-010,CKX19-013,DBS001,DBS006,DBS010,DBS014,DBS017,DBS020,DTX035,DTX038,DTX044,DTX049,DTX062,DTX066,DTX14,DTX56,EDL011,EDL019,EDL025}

do
md_bam_file=$Dir/markDuplicate/${sample}.so.md.bam
haplotype_file=$Dir/callSNP/ermanii/${sample}.g.vcf.gz
gatk HaplotypeCaller --native-pair-hmm-threads 2 -R $reference -I $md_bam_file -ploidy 4 -O $haplotype_file --output-mode EMIT_ALL_CONFIDENT_SITES --minimum-mapping-quality 30 --min-base-quality-score 20 -ERC GVCF
done
#diploid
Dir=/data_analysis/costataeEvolution
reference=/home/nwlab/experiment_data/RADseq/reference/35080Betula_pendula.fa

for sample in {BJF011,BJF016,BJF019,BJF031,BJF037,BJF046,BJF053,BJF060,BJF065,EDG004,EDG008,EDG026,EDG035,EDG039,EDG046,EDG23,FHS021,FHS026,FHS033,FHS039,FHSJ012,HR007,HR013,HR031,HR034,HR047,HR17,HTY005,HTY014,HTY018,HTY020}

do
md_bam_file=$Dir/markDuplicate/${sample}.so.md.bam
haplotype_file=$Dir/callSNP/costata/${sample}.g.vcf.gz
gatk HaplotypeCaller --native-pair-hmm-threads 2 -R $reference -I $md_bam_file -ploidy 2 -O $haplotype_file --output-mode EMIT_ALL_CONFIDENT_SITES --minimum-mapping-quality 30 --min-base-quality-score 20 -ERC GVCF
done

#GenomicsDBImport
file=`cat /data_analysis/costataeEvolution/Result/01PolyAncestor_variant_calling/sample.list`

InputVCF=""
for V in $file
do
#echo ${V}
InputVCF="$InputVCF -V ${V}"
done

chrname="lcl\|Bpe_"
reference=/home/nwlab/experiment_data/RADseq/reference/35080Betula_pendula.fa

parallel -j 7 "gatk GenomicsDBImport -R $reference ${InputVCF} -L ${chrname}Chr{} --genomicsdb-workspace-path ./GenomicsDBImport/GenomicsDBImport_Chr{}" ::: `seq 1 14`

parallel -j 7 "gatk GenotypeGVCFs -R $reference -L ${chrname}Chr{} -V gendb://GenomicsDBImport/GenomicsDBImport_Chr{} -all-sites --sample-ploidy 4 -O ./GenotypeGVCFs/Chr{}.vcf.gz" ::: `seq 1 14`