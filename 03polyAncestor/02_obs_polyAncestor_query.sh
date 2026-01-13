#!/bin/bash

#usage polyAncestor.sh species_of_the_tetraploid allotetraploid_id

######## What species is the polyploid
species=$1

######## Polyploidy sample ID
poly=$2

######## Initial path(working directory)
WorkDir=/data_analysis/costataeEvolution/Result/02PolyAncestor_main_analysis
SampleDir=$WorkDir/${species}/${poly}

mkdir -p $SampleDir

echo
echo "Species:" $species
echo "Sample ID:" $poly
echo "Output folder:" $SampleDir
echo


cd $SampleDir

mkdir MSA Iteration_cos Iteration_ash Iteration_bugg Iteration_len

bcftools view -s Alnus_glutinosa_ERR2103731,lenta_ERR2026272,nana_ERR2026268,occidentalis_ERR2026271,pendula_ERR2026097,populifolia_ERR2026269,EDL21-009-platy,XYB21-001-mcalli,TZZ78,CKX18-006,QSLC10,${poly} -O z -o ${poly}.snp.NO_VARIATION.QC.vcf.gz /data_analysis/costataeEvolution/Result/01PolyAncestor_variant_calling/07all.snp.NO_VARIATION.QC.nocall0.1.vcf.gz 

echo "subset ${poly} finished"

#convert vcf file to MSA
source ~/miniconda3/bin/activate pyvcf
python /data_analysis/costataeEvolution/Result/02PolyAncestor_main_analysis/vcf2fasta.py ${poly}.snp.NO_VARIATION.QC.vcf.gz ./MSA

#using trimAl to filter
for i in ./MSA/*.fa;do
   sed -i '2~2 s/N/-/g' ${i}
   bn=`basename ${i} .fa`
   trimal -in ${i} \
          -out ./MSA/${bn}.clean.1.fasta \
          -automated1 

   #convert to single line fasta
   seqkit seq -w 0 ./MSA/${bn}.clean.1.fasta > ./MSA/${bn}.clean.fasta
   rm ./MSA/${bn}.clean.1.fasta
   sed -i '2~2 s/-/N/g' ./MSA/${bn}.clean.fasta
done

#Check the length of MSA
total_length=0
file_count=0

for i in ./MSA/*.clean.fasta; do
    length=$(awk '{print length}' "$i" | head -2 | tail -1)
    total_length=$((total_length + length))
    file_count=$((file_count + 1))
done

if [ "$file_count" -gt 0 ]; then
    average_length=$(echo "scale=2; $total_length / $file_count" | bc)
    echo "number of RAD loci: $file_count" >> ./Log.MSA
    echo "total length: $total_length" >> ./Log.MSA
    echo "average length: $average_length" >> ./Log.MSA
else
    echo "Error cannot find the file"
fi

###################
## iteration 1 bugg

ref=bugg
for k in ./MSA/*.clean.fasta; do
   python3 /data_analysis/costataeEvolution/Result/02PolyAncestor_main_analysis/01s_polarizeTETRA.py $k 12 "CKX18-006" ${poly}
   echo MSAstep $k finished
done

mv ./MSA/*ALT.fasta ./Iteration_${ref}/
mv ./MSA/*stats ./Iteration_${ref}/

#Phylogenetic inference for each individual locus
ls ./Iteration_${ref}/*ALT.fasta > tree.lst
cat tree.lst | parallel -j 14 "iqtree --ufboot 1000 -alrt 1000 -s {} \
                                --seqtype DNA -m MFP -T 1 \
                                --seed 1234 --redo \
                                --safe --allnni -bnni -o Alnus_glutinosa_ERR2103731"
for i in {bionj,ALT.fasta.model.gz,contree,mldist,splits.nex,treefile,ufboot};do rm ./Iteration_${ref}/*${i} ;done

#Filter the parsimony-informative ect
/data_analysis/costataeEvolution/Result/02PolyAncestor_main_analysis/02_obs_polyAncestor_query.get_number_informative_MSA.sh ${ref} ${poly}

grep -v "FAIL" _AUX_${ref} | cut -f 1 -d " " > PASS.loci

#Summary all passed trees generated from each polarized MSA
for id in `cat PASS.loci`;do 
    cat ./Iteration_${ref}/${id}.clean.fa-ALT.fasta.iqtree | grep -A 2 'Consensus tree in newick format:' | tail -1
done > ${poly}Iteration_${ref}_astral_infile.txt


#polarized allotetraploid pairs with the diploid species
source ~/miniconda3/bin/activate base
Rscript --no-save /data_analysis/costataeEvolution/Result/02PolyAncestor_main_analysis/05_Polyploid_pairing_analysis.R \
                            ${poly}Iteration_${ref}_astral_infile.txt \
                            ${poly} \
                            "Alnus_glutinosa_ERR2103731" \
                             ./ \
                            "sister_ID_${poly}Iteration_${ref}.txt"
							

##################
## iteration 2 costata

ref=cos
for k in ./MSA/*.clean.fasta; do
   python3 /data_analysis/costataeEvolution/Result/02PolyAncestor_main_analysis/01s_polarizeTETRA.py $k 12 "QSLC10" ${poly}
   echo MSAstep $k finished
done

mv ./MSA/*ALT.fasta ./Iteration_${ref}/
mv ./MSA/*stats ./Iteration_${ref}/

#Phylogenetic inference for each individual locus
ls ./Iteration_${ref}/*ALT.fasta > tree.lst
cat tree.lst | parallel -j 14 "iqtree --ufboot 1000 -alrt 1000 -s {} \
                                --seqtype DNA -m MFP -T 1 \
                                --seed 1234 --redo \
                                --safe --allnni -bnni -o Alnus_glutinosa_ERR2103731"
for i in {bionj,ALT.fasta.model.gz,contree,mldist,splits.nex,treefile,ufboot};do rm ./Iteration_${ref}/*${i} ;done

#Summary all trees generated from each polarized MSA
for id in `cat PASS.loci`;do 
    cat ./Iteration_${ref}/${id}.clean.fa-ALT.fasta.iqtree | grep -A 2 'Consensus tree in newick format:' | tail -1
done > ${poly}Iteration_${ref}_astral_infile.txt

#polarized allotetraploid pairs with the diploid species
source ~/miniconda3/bin/activate base
Rscript --no-save /data_analysis/costataeEvolution/Result/02PolyAncestor_main_analysis/05_Polyploid_pairing_analysis.R \
                            ${poly}Iteration_${ref}_astral_infile.txt \
                            ${poly} \
                            "Alnus_glutinosa_ERR2103731" \
                             ./ \
                            "sister_ID_${poly}Iteration_${ref}.txt"
							
###################
## iteration 3 ash

ref=ash
for k in ./MSA/*.clean.fasta; do
   python3 /data_analysis/costataeEvolution/Result/02PolyAncestor_main_analysis/01s_polarizeTETRA.py $k 12 "TZZ78" ${poly}
   echo MSAstep $k finished
done

mv ./MSA/*ALT.fasta ./Iteration_${ref}/
mv ./MSA/*stats ./Iteration_${ref}/

#Phylogenetic inference for each individual locus
ls ./Iteration_${ref}/*ALT.fasta > tree.lst
cat tree.lst | parallel -j 14 "iqtree --ufboot 1000 -alrt 1000 -s {} \
                                --seqtype DNA -m MFP -T 1 \
                                --seed 1234 --redo \
                                --safe --allnni -bnni -o Alnus_glutinosa_ERR2103731"
for i in {bionj,ALT.fasta.model.gz,contree,mldist,splits.nex,treefile,ufboot};do rm ./Iteration_${ref}/*${i} ;done

for id in `cat PASS.loci`;do 
    cat ./Iteration_${ref}/${id}.clean.fa-ALT.fasta.iqtree | grep -A 2 'Consensus tree in newick format:' | tail -1
done > ${poly}Iteration_${ref}_astral_infile.txt

#polarized allotetraploid pairs with the diploid species
source ~/miniconda3/bin/activate base
Rscript --no-save /data_analysis/costataeEvolution/Result/02PolyAncestor_main_analysis/05_Polyploid_pairing_analysis.R \
                            ${poly}Iteration_${ref}_astral_infile.txt \
                            ${poly} \
                            "Alnus_glutinosa_ERR2103731" \
                             ./ \
                            "sister_ID_${poly}Iteration_${ref}.txt"
							

echo 
echo $poly "has finished."
date -u
echo

###################
## iteration 4 lenta

ref=len
for k in ./MSA/*.clean.fasta; do
   python3 /data_analysis/costataeEvolution/Result/02PolyAncestor_main_analysis/01s_polarizeTETRA.py $k 12 "lenta_ERR2026272" ${poly}
   echo MSAstep $k finished
done

if [ ! -d ./Iteration_len ]; then
mkdir ./Iteration_len
fi
done

mv ./MSA/*ALT.fasta ./Iteration_${ref}/
mv ./MSA/*stats ./Iteration_${ref}/

#Phylogenetic inference for each individual locus
ls ./Iteration_${ref}/*ALT.fasta > tree.lst
cat tree.lst | parallel -j 14 "iqtree --ufboot 1000 -alrt 1000 -s {} \
                                --seqtype DNA -m MFP -T 1 \
                                --seed 1234 --redo \
                                --safe --allnni -bnni -o Alnus_glutinosa_ERR2103731"
for i in {bionj,ALT.fasta.model.gz,contree,mldist,splits.nex,treefile,ufboot};do rm ./Iteration_${ref}/*${i} ;done

for id in `cat PASS.loci`;do 
    cat ./Iteration_${ref}/${id}.clean.fa-ALT.fasta.iqtree | grep -A 2 'Consensus tree in newick format:' | tail -1
done > ${poly}Iteration_${ref}_astral_infile.txt

#polarized allotetraploid pairs with the diploid species
source ~/miniconda3/bin/activate base
Rscript --no-save /data_analysis/costataeEvolution/Result/02PolyAncestor_main_analysis/05_Polyploid_pairing_analysis.R \
                            ${poly}Iteration_${ref}_astral_infile.txt \
                            ${poly} \
                            "Alnus_glutinosa_ERR2103731" \
                             ./ \
                            "sister_ID_${poly}Iteration_${ref}.txt"
							

echo 
echo $poly "has finished."
date -u
echo

