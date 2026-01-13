#!/bin/bash
################## Remember initial path
SRCDIR_INI=$(pwd)

################## Folder containing results from ABC analysis
AA=/mnt/sda/dingjunyi/data_analysis/costataeEvolution/Result/05s_ABC


################## Simulation conditions (polyploidization models)
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
pop[4]="contact_zone_albo_erman"
pop[5]="contact_zone_uti_albo"

################# Input file (results from ABC simulations)


################################## Compute L2 distance

for j in `seq 1 1 8`; do
    for k in `seq 4 1 5`; do               # populations
    ################## priors folder
	    if [ $k -eq 4 ]; then
		    PP=/mnt/sda/dingjunyi/data_analysis/costataeEvolution/Result/05s_ABC/${model[$j]}/ermanii
            Infile=summary_ABC_500_simulations.txt
            AA_MAIN=${AA}/${model[$j]}/ermanii
            PP_GO=${PP}/ABC_PRIORS_${model[$j]}_ermanii
            RR_MAIN=/mnt/sda/dingjunyi/data_analysis/costataeEvolution/Result/05s_ABC_geographic_distribution/${model[$j]}/${pop[$k]}
	    elif [ $k -eq 5 ];then
		    PP=/mnt/sda/dingjunyi/data_analysis/costataeEvolution/Result/05s_ABC/${model[$j]}/albosinensis
            Infile=summary_ABC_500_simulations.txt
            AA_MAIN=${AA}/${model[$j]}/albosinensis
            PP_GO=${PP}/ABC_PRIORS_${model[$j]}_albosinensis
            RR_MAIN=/mnt/sda/dingjunyi/data_analysis/costataeEvolution/Result/05s_ABC_geographic_distribution/${model[$j]}/${pop[$k]}
		else
            PP=/mnt/sda/dingjunyi/data_analysis/costataeEvolution/Result/05s_ABC/${model[$j]}/${pop[$k]}
            Infile=summary_ABC_500_simulations.txt
            AA_MAIN=${AA}/${model[$j]}/${pop[$k]}
            PP_GO=${PP}/ABC_PRIORS_${model[$j]}_${pop[$k]}
            RR_MAIN=/mnt/sda/dingjunyi/data_analysis/costataeEvolution/Result/05s_ABC_geographic_distribution/${model[$j]}/${pop[$k]}
        fi
       ################## Observed pairing profiles
       PP_Obs=/mnt/sda/dingjunyi/data_analysis/costataeEvolution/Result/05s_ABC_geographic_distribution/00_PAIRING_PROFILES_${pop[$k]}.txt

       if [ ! -d $RR_MAIN ]; then
       mkdir -p $RR_MAIN
       fi

       echo
       echo ${model[$j]} $REP $NLOCUS
       echo $AA_MAIN
       echo $RR_MAIN
       echo $PP_GO
	   echo $PP_Obs
       echo

       python3 33_L2_ABC_results.py $AA_MAIN \
                                    $RR_MAIN \
                                    $PP_Obs \
                                    $PP_GO \
                                    $Infile


       # clean output file (remove square brackets; replace multiple spaces by single space)
       sed 's/[[]//g' $RR_MAIN/L2-norm_ABC_1000_simulations.txt > $RR_MAIN/_aux1.txt
       sed 's/[]]//g' $RR_MAIN/_aux1.txt | tr -s ' ' > $RR_MAIN/L2-norm_ABC_1000_simulations.txt
       rm -f $RR_MAIN/_aux1.txt
    done
done


echo
echo "Done!"
date -u
