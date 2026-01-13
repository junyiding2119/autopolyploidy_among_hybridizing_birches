#!/bin/bash
model[1]="AAAA"
model[2]="BBBB"
model[3]="CCCC"
model[4]="AacAacAacAac_A"
model[5]="AacAacAacAac_C"
model[6]="AABB"
model[7]="AACC"
model[8]="CCBB"

pop[1]="albosinensis"
pop[2]="ermanii"
pop[3]="utilis"

OUTfile=summary_ABC_500_simulations.txt

for k in `seq 1 1 8`; do 
    for j in `seq 1 1 3`; do 

        AA_MAIN=/mnt/sda/dingjunyi/data_analysis/costataeEvolution/Result/05s_ABC/${model[k]}/${pop[j]}
        echo
        echo ${model[k]} ${pop[j]}
        echo $AA_MAIN
        echo
		if [ -d $AA_MAIN/$OUTfile ]; then rm $AA_MAIN/$OUTfile; fi
        
        echo "ABC_No" "pendula" "platyphylla" "pend_platy" "populifolia" "pop_out" "nana" "nana_out" "occidentalis" "occ_out" "cos" "ash" "cos_ash" "bug" "bug_out" "mca" "mca_out" "len" "basal" "pendula" "platyphylla" "pend_platy" "populifolia" "pop_out" "nana" "nana_out" "occidentalis" "occ_out" "cos" "ash" "cos_ash" "bug" "bug_out" "mca" "mca_out" "len" "basal" "pendula" "platyphylla" "pend_platy" "populifolia" "pop_out" "nana" "nana_out" "occidentalis" "occ_out" "cos" "ash" "cos_ash" "bug" "bug_out" "mca" "mca_out" "len" "basal" >> $AA_MAIN/$OUTfile
         
        for n in `seq 1 1 500`; do
        
              InFile=$AA_MAIN/$n/03_PAIRINGfrequencies/sister_ID_analysis_priors-JOINT_summary.txt
              sim_ash="$(cat $InFile | head -n1)"
              sim_bug="$(cat $InFile | tail -n2| head -n1)"
              sim_cos="$(cat $InFile | tail -n1)"
        
              echo $n $sim_ash $sim_bug $sim_cos >> $AA_MAIN/$OUTfile
        done
    done
done