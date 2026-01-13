######
#STEP1 generate simulated gene phylogenies with SimPhy
######
simphy -rs 25 \
        -rl f:50 \
        -rg 1 \
        -si f:2 \
        -sp f:400000 \
        -su ln:-17.27461,0.6931472 \
        -hh f:1 \
        -hs ln:1.5,1 \
        -hl ln:1.551533,0.6931472 \
        -hg ln:1.5,1 \
        -cs 55277 \
        -v 3 \
        -o /data_analysis/costataeEvolution/Result/04simphy_ILS/A_modILS_birch \
        -ot 0 \
        -op 1 \
        -od 1 \
        -s '(alnus_SIM:80,(lenta_SIM:70,((((((ashburneri_SIM:10,ANCESTRAL_A:10):10,(costata_SIM:10,ANCESTRAL_C:10):10):10,ANCESTRAL_ANCac:30):10,(buggsii_SIM:20,ANCESTRAL_B:20):20):10,mcallisteri_SIM:50):10,((((pendula_SIM:12,platyphylla_SIM:12):12,populifolia_SIM:24):12,nana_SIM:36):12,occidentalis_SIM:48):12):10):10);'
		
#####		
#STEP2 fix tree tips  
#####

#The 02_fix_tree_tips_lowILS_query.sh script can be found in https://github.com/LLN273/Complex-Polyploids/blob/main/09_Model_Testing_1_Polyploidization_and_hybridization_simulations/02_fix_tree_tips_lowILS_query.sh

#remember initial path
SRCDIR_INI=$(pwd)

# Input folder (root) [results from simphy]
AA=/data_analysis/costataeEvolution/Result/04simphy_ILS

# Input/output subfolder (different simulation conditions)
INfolder[1]=A_modILS_birch

## Sample list
SPL=/data_analysis/costataeEvolution/Result/04simphy_ILS/sample_list

# Number of species
NSPEC[1]=15

## number of replicates
REP=25

# number of loci/genes		
NLOCUS=50

for k in `seq 1 1 1`; do 			# different simulation conditions

      echo
      GO_NSPEC=${NSPEC[$k]}
      GO_FOLDER=${INfolder[$k]}
      GO_AA=${AA}/${GO_FOLDER}
      GO_RR=${AA}/${GO_FOLDER}_CLEAN
      mkdir -p $GO_RR

      cd $GO_RR
      rsync -ah $GO_AA/*.command .
      rsync -ah $GO_AA/*.db .
      rsync -ah $GO_AA/*.params .


      for r_aux in $( eval echo {01..${REP}} ); do		

         mkdir -p $GO_RR/$r_aux
         cd $GO_RR/$r_aux

         cd $SRCDIR_INI
         $SRCDIR_INI/02_fix_tree_tips_lowILS_query.sh $GO_AA/$r_aux \
                                                      $GO_RR/$r_aux \
                                                      $r_aux \
                                                      $GO_NSPEC \
                                                      $NLOCUS \
                                                      $SPL
                                          
         echo

      done
done 

#####
#decrease ILS
#####

#!/bin/bash


## Script used to decrease ILS by replacing random genes with genes reflecting true gene phylogeny



#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)

# Input folder (root) [results ASTRAL analysis]
AA=/data_analysis/costataeEvolution/Result/04simphy_ILS

# Input subfolder (different simulation conditions)
INfolder[1]=A_modILS_birch_CLEAN
INfolder[2]=${INfolder[1]}
INfolder[3]=${INfolder[1]}


# Output subfolder (different simulation conditions)
OUTfolder[1]=A_modILS_birch_CLEAN_reducedILS_30perct     # 30% genes replaced 
OUTfolder[2]=A_modILS_birch_CLEAN_reducedILS_40perct     # 40% genes replaced
OUTfolder[3]=A_modILS_birch_CLEAN_reducedILS_50perct     # 50% genes replaced

## new gene phylogeny (must replicate species tree if we aim at reducing ILS)
GP_GOLD=${AA}/00_reference_GENE_MAN.trees


## proportion of genes to be replaced
GP_prop[1]="0.3"
GP_prop[2]="0.4"
GP_prop[3]="0.5"

## number of replicates
REP=25

# number of loci/genes
NLOCUS=50


echo

for k in `seq 1 1 3`; do                      # different simulation conditions

      echo
      GO_FOLDER=${INfolder[$k]}
      GO_OUT=${OUTfolder[$k]}
      GO_AA=${AA}/${GO_FOLDER}
      GO_RR=${AA}/${GO_OUT}
      mkdir -p $GO_RR


      cd $GO_RR
      rsync -ah $GO_AA/*.command .
      rsync -ah $GO_AA/*.db .
      rsync -ah $GO_AA/*.params .


      for r_aux in $( eval echo {01..${REP}} ); do

         mkdir -p $GO_RR/$r_aux

         echo
                     echo
         echo $r_aux $NLOCUS
         echo ${GP_prop[$k]} $GP_GOLD
         echo $GO_AA/$r_aux
         echo $GO_RR/$r_aux


         cd $SRCDIR_INI
         $SRCDIR_INI/03_adjust_ILS_NEW_query.sh $GO_AA/$r_aux \
                                                      $GO_RR/$r_aux \
                                                      $r_aux \
                                                      $NLOCUS \
                                                      $GP_GOLD \
                                                      ${GP_prop[$k]}






         echo

      done
done


#####
#STEP4 INDELIble
#####

# path to INDELIble
INDELIBLE="/home/nwlab/biosoft/software/INDELibleV1.03/src/"
export PATH="$INDELIBLE:$PATH"

# simphy's INDELIble wrapper
IDW=/home/nwlab/biosoft/software/SimPhy_1.0.2/scripts/INDELIble_wrapper.pl


#### paths and folder names

# remember initial path
SRCDIR_INI=$(pwd)      

# reduced ILS

reduced_ILS=40

# Input folder
AA=/data_analysis/costataeEvolution/Result/04simphy_ILS/A_modILS_birch_CLEAN_reducedILS_${reduced_ILS}perct	# reduced ILS (40%); 1500bp read length


# Configuration file
CONFIG=${SRCDIR_INI}/INDELible_complex_1500.txt			# read length reduced to (mean=1500,st=100); used in P11 sims
		

# Run perl wrapper
# Run as: INDELIble_wrapper.pl directory input_config seed number-of-cores(0 means not use GNU Parallel)

perl $IDW $AA $CONFIG 22 0