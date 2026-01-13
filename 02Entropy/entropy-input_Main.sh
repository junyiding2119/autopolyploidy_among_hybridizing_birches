Rscript ./inputdataformat.R 05.allInd.non_coding.bisnp.QC.MAF0.01.nocall0.1.vcf

for r in {1..4}; do for K in {2..8};
do
nohup ./entropy-query.sh 05.allInd.non_coding.bisnp.QC.MAF0.01.nocall0.1.mpgl $K $r > log.K$K.r$r.txt &

done;done


for r in {3..5}; do for K in {2..8};
do
#Mem 458.50G CPU 41%
sbatch ./entropy-query.sh 05.allInd.non_coding.bisnp.QC.MAF0.01.nocall0.1.mpgl $K $r

done;done

#potential scale reduction factor,assess convergence
for i in {1..8};do estpost.entropy -p q -s 4 mcmcoutchain?.k${i}.hdf5 -o qmcmcdiag.${k}.txt;done

estpost.entropy -p gprob -s 0 mcmcoutchain?.k6.hdf5 -o genoest.txt