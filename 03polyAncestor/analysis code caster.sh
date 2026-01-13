bcftools view -s Alnus_glutinosa_ERR2103731,lenta_ERR2026272,nana_ERR2026268,occidentalis_ERR2026271,pendula_ERR2026097,populifolia_ERR2026269,EDL21-009-platy,XYB21-001-mcalli,TZZ78,CKX18-006,QSLC10,TBS100,CBS023,DCX024 -O z -o caster.snp.NO_VARIATION.QC.vcf.gz /data_analysis/costataeEvolution/Result/01PolyAncestor_variant_calling/07all.snp.NO_VARIATION.QC.nocall0.1.vcf.gz

gatk --java-options "-Xmx3g -Xms3g" IndexFeatureFile -I caster.snp.NO_VARIATION.QC.vcf.gz
gatk --java-options "-Xmx3g -Xms3g" VariantsToTable -V caster.snp.NO_VARIATION.QC.vcf.gz -GF GT -O nSNP_NO_VARIATION.miss0_1.phylo.txt

#check
grep "[ATGC*.]\{2,\}" nSNP_NO_VARIATION.miss0_1.phylo.txt

python3 table2caster.py nSNP_NO_VARIATION.miss0_1.phylo.txt caster.polyAncestor.fasta

caster-site -t 2 -i caster.polyAncestor.fasta -o caster.polyAncestor.tree --root Alnus_glutinosa_ERR2103731.GT

python pick_seq.py caster.polyAncestor.fasta caster.polyAncestor.oneSeq.fasta

sed -i "0~2s/\./N/g;0~2s/\*/N/g" caster.polyAncestor.oneSeq.fasta

raxml -T 8 -f d -g caster.polyAncestor.tree -m GTRGAMMA -s caster.polyAncestor.oneSeq.fasta -n caster.polyAncestor.branch.tree -p 12345 -o Alnus_glutinosa_ERR2103731.GT

#python
conda activate ete3
from ete3 import Tree

t = Tree('(alnus_SIM:80,(lenta_SIM:70,((((((ashburneri_SIM:10,ANCESTRAL_A:10):10,(costata_SIM:10,ANCESTRAL_C:10):10):10,ANCESTRAL_ANCac:30):10,(buggsii_SIM:20,ANCESTRAL_B:20):20):10,mcallisteri_SIM:50):10,((((pendula_SIM:20,platyphylla_SIM:20):10,populifolia_SIM:30):10,nana_SIM:40):10,occidentalis_SIM:50):10):10):10);')

t.set_outgroup('alnus_SIM')  

print('Non-ultrametric rooted tree:')
print(t.write())

t.convert_to_ultrametric()
print('Ultrametric rooted tree:')
print(t.write())