#!/bin/bash

SRCDIR_INI=$(pwd) 

repstart=$1

repend=$2


Input_bestrun_file=$SRCDIR_INI/bootstrap/bestrun

echo "Start and End Pseudo data sets: $repstart $repend"

for n in `seq $repstart 1 $repend`; do 			

   echo
   echo "Pseudo data sets:" $n
  
   boot_file=$SRCDIR_INI/bootstrap/model${EVOLM}_maxL.bootstrap/model${EVOLM}_maxL.bootstrap_${n}

   $SRCDIR_INI/03fastsimcola26_confidence_intervals_query.sh $boot_file $Input_bestrun_file ${EVOLM}
   
   cd $boot_file
   
   ~/data_abalysis/DingJunYi/costataeEvolution/Result/07Demographic_modeling_Main/albosinensis/02fsc-selectbestrun.sh model${EVOLM}_maxL.bootstrap
   
   cd $SRCDIR_INI
done