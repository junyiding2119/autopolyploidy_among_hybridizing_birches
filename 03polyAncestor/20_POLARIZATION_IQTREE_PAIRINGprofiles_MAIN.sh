#!/bin/bash
#

## Read unpolarized MSAs and priors; generate polarized MSAs and apply polyploidization and hybridization model; run IQTREE2 and generate pairing profiles.

###################################### paths and folder names

###### remember initial path
SRCDIR_INI=$(pwd)                                           	 

###### input folder (MSAs, before polarization)
AA=${1:?msg}

###### output folder (root)
RR_MAIN=${2:?msg}

###### Simulation conditions (polyploidization models)
myMODEL=${3:?msg}

######  Number of replicates	
REP=${4:?msg}

######  Number of loci/genes		
NLOCUS=${5:?msg}

###### PRIORS
PP=${6:?msg}

SNIC_TMP=${7:?msg}


###### Gene list
GL_large=$SRCDIR_INI/00_GENE_LIST.txt




######################################  IQTree parameters

# Number of bootstraps
BS=1000

#### Evolutionary model
MODEL_IQTREE2=MFP		# MFP >>> ModelFinder used to determine best substitution model (but no correction for heterotachy)

####### MODEL RATE: List of rate heterogeneity among sites (when using Modelfinder and GHOST)
MODEL_IQTREE2_RATE="F*H"				

####### model name to be included in output file(s)
MODEL_outfilename=MFP_ModelFinder



######################################  Polyploidization and hybridization model

if [[ ${myMODEL} = "AAAA" ]] ; then
	
   ############## AAAA (B. ashburneri autopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="ANCESTRAL_A"	# ashburneri
   GO_PAR[2]="ANCESTRAL_A"	# ashburneri
   GO_PAR[3]="ANCESTRAL_A"	# ashburneri
   GO_PAR[4]="ANCESTRAL_A"	# ashburneri
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="pendula_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_C ANCESTRAL_C"	# block 0 
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_B ANCESTRAL_B"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_C ANCESTRAL_B"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_A ANCESTRAL_A"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 5	
   #read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_A ANCESTRAL_A"	# block 6
   #read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_A ANCESTRAL_A"	# block 7	
   #read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_A ANCESTRAL_A"	# block 8
   
   ### Extra introgression 
   Extra_introgr_species="PUBESCENS_PEN"

fi


if [[ ${myMODEL} = "CCCC" ]] ; then
	
   ############## CCCC (B. costata autopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="ANCESTRAL_C"	# costata
   GO_PAR[2]="ANCESTRAL_C"	# costata
   GO_PAR[3]="ANCESTRAL_C"	# costata
   GO_PAR[4]="ANCESTRAL_C"	# costata
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="pendula_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "ANCESTRAL_C ANCESTRAL_C ANCESTRAL_A ANCESTRAL_A"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "ANCESTRAL_C ANCESTRAL_C ANCESTRAL_B ANCESTRAL_B"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "ANCESTRAL_C ANCESTRAL_C ANCESTRAL_A ANCESTRAL_B"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "ANCESTRAL_C ANCESTRAL_C ANCESTRAL_C ANCESTRAL_C"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "ANCESTRAL_C ANCESTRAL_C ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 5	
   #read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "ANCESTRAL_C ANCESTRAL_C ANCESTRAL_C ANCESTRAL_C"	# block 6
   #read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "ANCESTRAL_C ANCESTRAL_C ANCESTRAL_C ANCESTRAL_C"	# block 7	
   #read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "ANCESTRAL_C ANCESTRAL_C ANCESTRAL_C ANCESTRAL_C"	# block 8
   
   ### Extra introgression
   Extra_introgr_species="PUBESCENS_PEN"

fi


if [[ ${myMODEL} = "BBBB" ]] ; then
	
   ############## BBBB (B. buggsii autopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="ANCESTRAL_B"	# buggsii
   GO_PAR[2]="ANCESTRAL_B"	# buggsii
   GO_PAR[3]="ANCESTRAL_B"	# buggsii
   GO_PAR[4]="ANCESTRAL_B"	# buggsii
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="pendula_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "ANCESTRAL_B ANCESTRAL_B ANCESTRAL_A ANCESTRAL_A"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "ANCESTRAL_B ANCESTRAL_B ANCESTRAL_C ANCESTRAL_C"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "ANCESTRAL_B ANCESTRAL_B ANCESTRAL_A ANCESTRAL_C"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "ANCESTRAL_B ANCESTRAL_B ANCESTRAL_B ANCESTRAL_B"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "ANCESTRAL_B ANCESTRAL_B ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 5	
   #read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "ANCESTRAL_B ANCESTRAL_B ANCESTRAL_B ANCESTRAL_B"	# block 6
   #read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "ANCESTRAL_B ANCESTRAL_B ANCESTRAL_B ANCESTRAL_B"	# block 7	
   #read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "ANCESTRAL_B ANCESTRAL_B ANCESTRAL_B ANCESTRAL_B"	# block 8
   
   ### Extra introgression
   Extra_introgr_species="PUBESCENS_PEN"

fi




if [ ${myMODEL} = "AACC" ] || [ ${myMODEL} = "AACC_HE" ] ; then
	
   ############## AACC (B. ashburneri/B. costata allopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="ANCESTRAL_A"	# ashburneri
   GO_PAR[2]="ANCESTRAL_A"	# ashburneri
   GO_PAR[3]="ANCESTRAL_C"	# costata
   GO_PAR[4]="ANCESTRAL_C"	# costata
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="pendula_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_B ANCESTRAL_B"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_C ANCESTRAL_C"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_C ANCESTRAL_B"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "ANCESTRAL_B ANCESTRAL_B ANCESTRAL_C ANCESTRAL_C"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 5	
   #read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_C ANCESTRAL_C"	# block 6
   #read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_C ANCESTRAL_C"	# block 7	
   #read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "POLY_ASH POLY_ASH POLY_COS POLY_COS"	# block 8
   
   ### Extra introgression 
   Extra_introgr_species="PUBESCENS_PEN"

fi


if [ ${myMODEL} = "AABB" ] || [ ${myMODEL} = "AABB_HE" ] ; then
	
   ############## AABB (B. ashburneri/B. buggsii allopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="ANCESTRAL_A"	# ashburneri
   GO_PAR[2]="ANCESTRAL_A"	# ashburneri
   GO_PAR[3]="ANCESTRAL_B"	# buggsii
   GO_PAR[4]="ANCESTRAL_B"	# buggsii
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="pendula_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_C ANCESTRAL_C"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_B ANCESTRAL_B"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_B ANCESTRAL_C"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "ANCESTRAL_C ANCESTRAL_C ANCESTRAL_B ANCESTRAL_B"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 5	
   #read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_B ANCESTRAL_B"	# block 6
   #read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_B ANCESTRAL_B"	# block 7	
   #read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "POLY_ASH POLY_ASH POLY_BUG POLY_BUG"	# block 8
   
   ### Extra introgression 
   Extra_introgr_species="PUBESCENS_PEN"

fi


if [ ${myMODEL} = "CCBB" ] || [ ${myMODEL} = "CCBB_HE" ] ; then
	
   ############## CCBB (B. costata/B. buggsii allopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="ANCESTRAL_C"	# costata
   GO_PAR[2]="ANCESTRAL_C"	# costata
   GO_PAR[3]="ANCESTRAL_B"	# buggsii
   GO_PAR[4]="ANCESTRAL_B"	# buggsii
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="pendula_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "ANCESTRAL_C ANCESTRAL_C ANCESTRAL_A ANCESTRAL_A"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "ANCESTRAL_C ANCESTRAL_C ANCESTRAL_B ANCESTRAL_B"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "ANCESTRAL_C ANCESTRAL_C ANCESTRAL_B ANCESTRAL_A"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "ANCESTRAL_A ANCESTRAL_A ANCESTRAL_B ANCESTRAL_B"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "ANCESTRAL_C ANCESTRAL_C ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 5	
   #read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "ANCESTRAL_C ANCESTRAL_C ANCESTRAL_B ANCESTRAL_B"	# block 6
   #read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "ANCESTRAL_C ANCESTRAL_C ANCESTRAL_B ANCESTRAL_B"	# block 7	
   #read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "POLY_COS POLY_COS POLY_BUG POLY_BUG"	# block 8
   
   ### Extra introgression 
   Extra_introgr_species="PUBESCENS_PEN"

fi


if [[ ${myMODEL} = "AacAacAacAac_A" ]] ; then
	
   ############## AacAacAacAac (ashburneri-costata ancestor, autopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="ANCESTRAL_ANCac"	# ashburneri-costata ancestor
   GO_PAR[2]="ANCESTRAL_ANCac"	# ashburneri-costata ancestor
   GO_PAR[3]="ANCESTRAL_ANCac"	# ashburneri-costata ancestor  
   GO_PAR[4]="ANCESTRAL_ANCac"	# ashburneri-costata ancestor
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="pendula_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_C ANCESTRAL_C"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_B ANCESTRAL_B"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_C ANCESTRAL_B"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 5	
   #read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 6
   #read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 7	
   #read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 8
   
   ### Extra introgression 

   Extra_introgr_species="ANCESTRAL_A"

fi

if [[ ${myMODEL} = "AacAacAacAac_C" ]] ; then
	
   ############## AacAacAacAac (ashburneri-costata ancestor, autopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="ANCESTRAL_ANCac"	# ashburneri-costata ancestor
   GO_PAR[2]="ANCESTRAL_ANCac"	# ashburneri-costata ancestor
   GO_PAR[3]="ANCESTRAL_ANCac"	# ashburneri-costata ancestor  
   GO_PAR[4]="ANCESTRAL_ANCac"	# ashburneri-costata ancestor
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="pendula_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_A ANCESTRAL_A"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_B ANCESTRAL_B"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_A ANCESTRAL_B"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 5	
   #read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 6
   #read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 7	
   #read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac ANCESTRAL_ANCac"	# block 8
   
   ### Extra introgression 

   Extra_introgr_species="ANCESTRAL_C"

fi
############################################################################  START
   
################### Read priors

unset PRIOR_LIST
x=1

while read -r LINE                                                                      
do
   PRIOR_LIST[$x]=$LINE
   let "x+=1"
done < ${PP}

N_PRIORS=${#PRIOR_LIST[@]}



################### Loci of interest
cat $GL_large | head -${NLOCUS} > $SNIC_TMP/00_gene_list_TRIM_MAX-00000.txt
GO_GL=$SNIC_TMP/00_gene_list_TRIM_MAX-00000.txt



###################
################### POLARIZATION
###################
echo
echo "1. POLARIZATION..." 
date -u

   # PRIOR_LIST[1]  ILS prior
   # PRIOR_LIST[2]  platy-pendula gene flow
   # PRIOR_LIST[3]  H0
   # PRIOR_LIST[4]  H1
   # PRIOR_LIST[5]  H2
   # PRIOR_LIST[6]  H3
   # PRIOR_LIST[7]  H4
   # PRIOR_LIST[8]  H5
   # PRIOR_LIST[9]  H6 extra species block

# ILS 
GO_ILS_reduction=${PRIOR_LIST[1]}
GO_FOLDER=A_modILS_birch_CLEAN_reducedILS_${GO_ILS_reduction}perct
	  
# input folder
GO_AA=${AA}/${GO_FOLDER}
	  
# Number of species (including root species)
GO_NSPEC=15
	  
# platyphylla to pendula gene flow [0-1]
GO_HYB_FRACTION_2pen=${PRIOR_LIST[2]}
	 
# Hybridization levels for each block
GO_HYBRY_FRACTION0=${PRIOR_LIST[3]}
GO_HYBRY_FRACTION1=${PRIOR_LIST[4]}
GO_HYBRY_FRACTION2=${PRIOR_LIST[5]}
GO_HYBRY_FRACTION3=${PRIOR_LIST[6]}
GO_HYBRY_FRACTION4=${PRIOR_LIST[7]}
GO_HYBRY_FRACTION5=${PRIOR_LIST[8]}

	  
# Hybridization levels for extra species
GO_HYBRY_EXTRA=${PRIOR_LIST[9]}


######## POLAR 1: Reference sequence: B. ashburneri
echo "Polarization 1... B. ashburneri"
	  
# Polarizing reference sequence ID
GO_REFspecies_1="ashburneri_SIM"
		 
# polarization label
NALT=1

# output folder
GO_RR=$SNIC_TMP/01_POLAR/ALT${NALT}
rm -rf $GO_RR
mkdir -p $GO_RR

#run in parallel for different replicates 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)
            $SRCDIR_INI/21_polarize_MSA_query.sh $GO_AA $GO_RR $GO_NSPEC $GO_FOLDER $r_aux $NLOCUS "${GO_PAR[@]}" $GO_REFspecies_1 $NALT $GF_2pend_species $GO_HYB_FRACTION_2pen "${GO_HYBR0[@]}" ${GO_HYBRY_FRACTION0} \
                                  "${GO_HYBR1[@]}" "${GO_HYBR2[@]}" "${GO_HYBR3[@]}" "${GO_HYBR4[@]}" "${GO_HYBR5[@]}"   \
								${GO_HYBRY_FRACTION1} ${GO_HYBRY_FRACTION2} ${GO_HYBRY_FRACTION3} ${GO_HYBRY_FRACTION4} ${GO_HYBRY_FRACTION5} \
								${Extra_introgr_species} ${GO_HYBRY_EXTRA} ${myMODEL} ${SNIC_TMP} &
done
		 
wait


######## POLAR 2: Reference sequence: B. buggsii
echo "Polarization 2... B. buggsii"
	  
# Polarizing reference sequence ID
GO_REFspecies_1="buggsii_SIM"
	 
# polarization label
NALT=2

# output folder
GO_RR=$SNIC_TMP/01_POLAR/ALT${NALT}
rm -rf $GO_RR
mkdir -p $GO_RR

#run in parallel for different replicates 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)
            $SRCDIR_INI/21_polarize_MSA_query.sh $GO_AA $GO_RR $GO_NSPEC $GO_FOLDER $r_aux $NLOCUS "${GO_PAR[@]}" $GO_REFspecies_1 $NALT $GF_2pend_species $GO_HYB_FRACTION_2pen "${GO_HYBR0[@]}" ${GO_HYBRY_FRACTION0} \
                               "${GO_HYBR1[@]}" "${GO_HYBR2[@]}" "${GO_HYBR3[@]}" "${GO_HYBR4[@]}" "${GO_HYBR5[@]}" \
							${GO_HYBRY_FRACTION1} ${GO_HYBRY_FRACTION2} ${GO_HYBRY_FRACTION3} ${GO_HYBRY_FRACTION4} ${GO_HYBRY_FRACTION5}  \
							${Extra_introgr_species} ${GO_HYBRY_EXTRA} ${myMODEL} ${SNIC_TMP} &
done
		 
wait
	  
	  
######## POLAR 3: Reference sequence: B. costata
echo "Polarization 3... B. costata"
	  
# Polarizing reference sequence ID
GO_REFspecies_1="costata_SIM"
 
# polarization label
NALT=3

# output folder
GO_RR=$SNIC_TMP/01_POLAR/ALT${NALT}
rm -rf $GO_RR
mkdir -p $GO_RR

#run in parallel for different replicates 
		for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)
            $SRCDIR_INI/21_polarize_MSA_query.sh $GO_AA $GO_RR $GO_NSPEC $GO_FOLDER $r_aux $NLOCUS "${GO_PAR[@]}" $GO_REFspecies_1 $NALT $GF_2pend_species $GO_HYB_FRACTION_2pen "${GO_HYBR0[@]}" ${GO_HYBRY_FRACTION0} \
                               "${GO_HYBR1[@]}" "${GO_HYBR2[@]}" "${GO_HYBR3[@]}" "${GO_HYBR4[@]}" "${GO_HYBR5[@]}" \
							${GO_HYBRY_FRACTION1} ${GO_HYBRY_FRACTION2} ${GO_HYBRY_FRACTION3} ${GO_HYBRY_FRACTION4} ${GO_HYBRY_FRACTION5}  \
							${Extra_introgr_species} ${GO_HYBRY_EXTRA} ${myMODEL} ${SNIC_TMP} &
done
		 
wait
  


###################
################### IQ-TREE
###################
echo
echo "2. IQ-TREE..." 
date -u
		  
AA2=$SNIC_TMP/01_POLAR
RR2=$SNIC_TMP/02_IQTREE
	  
#GO_FOLDER=${INfolder[$k]}
#GO_POLAR=ALT
#GO_NLOCUS=${NLOCUS}
	  
#input file suffix
InFileSUF=TRUE_CLEAN-ALT.fasta
		  
# Outgroup sequence
OUTG="alnus_SIM"
	  
######## POLAR 1: Reference sequence: B. ashburneri
#echo "Polarization 1..."
		  
# polarization label
NALT=1

GO_AA2=${AA2}/ALT${NALT}
GO_RR2=${RR2}/ALT${NALT}
rm -rf $GO_RR2
mkdir -p $GO_RR2
			   
#run in parallel for different replicates 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)
	  		      $SRCDIR_INI/22_IQ-TREE2_SINGLEGenes_ultrafastBS_SimData_query.sh $GO_AA2 $InFileSUF $GO_RR2 $BS ${MODEL_outfilename} $MODEL_IQTREE2 $MODEL_IQTREE2_RATE $NLOCUS $OUTG $r_aux $NALT ${SNIC_TMP} &
done
			   
wait
			
			
			   
######## POLAR 2: Reference sequence: B. buggsii
#echo "Polarization 2..."
	  
# polarization label
NALT=2

GO_AA2=${AA2}/ALT${NALT}
GO_RR2=${RR2}/ALT${NALT}
rm -rf $GO_RR2
mkdir -p $GO_RR2
   
#run in parallel for different replicates 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)
		    $SRCDIR_INI/22_IQ-TREE2_SINGLEGenes_ultrafastBS_SimData_query.sh $GO_AA2 $InFileSUF $GO_RR2 $BS ${MODEL_outfilename} $MODEL_IQTREE2 $MODEL_IQTREE2_RATE $NLOCUS $OUTG $r_aux $NALT ${SNIC_TMP} &
done
		   
wait  
			   
			   
			   
######## POLAR 3: Reference sequence: B. costata
#echo "Polarization 3..."
	  
# polarization label
NALT=3

GO_AA2=${AA2}/ALT${NALT}
GO_RR2=${RR2}/ALT${NALT}
rm -rf $GO_RR2
mkdir -p $GO_RR2

#run in parallel for different replicates 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)
		    $SRCDIR_INI/22_IQ-TREE2_SINGLEGenes_ultrafastBS_SimData_query.sh $GO_AA2 $InFileSUF $GO_RR2 $BS ${MODEL_outfilename} $MODEL_IQTREE2 $MODEL_IQTREE2_RATE $NLOCUS $OUTG $r_aux $NALT ${SNIC_TMP} &
done
		 
wait    
			   

			   
###################
################### ## Create file with IQ-TEST2 consensus tree obtained for each gene | Use either maximum number of genes available (MAX) or random gene lists previously created
###################			   
echo
echo "3. Glean IQ-TREE gene trees..." 
date -u
	  
AA3=$SNIC_TMP/02_IQTREE
RR3=${RR_MAIN}/03_PAIRINGfrequencies
		 
# input file suffix
#file_aux=""
InFileSUF=TRUE_CLEAN_UFBoot_${MODEL_outfilename}.iqtree
			   
### read gene list

unset GENE_LIST
x=1

while read -r LINE                                                                      
do
   GENE_LIST[$x]=$LINE
   let "x+=1"
 done < ${GO_GL}

N_GENES=${#GENE_LIST[@]}
		 
		 
######## POLAR 1: Reference sequence: B. pendula
#echo "Polarization 1..."
	  
# polarization label
NALT=1
		 
GO_AA3=${AA3}/ALT${NALT}
GO_RR3=${RR3}/ALT${NALT}
rm -rf $GO_RR3
		 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)

   mkdir -p $GO_RR3/$r_aux
  
   ########output file
   OUTfile=GENE_TREES_${r_aux}_UFBoot_MAX-${N_GENES}.newick
   rm -f ${GO_RR3}/$r_aux/$OUTfile
			   
   for i in `seq 1 1 $N_GENES`; do                 # cycle through genes
			   
                  #input file
                  GO_GENE=${GENE_LIST[$i]}
                  InFile=data_${GO_GENE}_${InFileSUF}
                  #InFile=data_${GO_GENE}_TRUE_CLEAN_UFBoot_${MODEL_outfilename}.iqtree
			
                  # get IQ-TEST2 consensus tree for each gene           
                  cat $GO_AA3/$r_aux/${InFile} | grep -A 2 'Consensus tree in newick format:' | tail -1 >> ${GO_RR3}/$r_aux/$OUTfile
   done
done
			
# get genes from all replicates
cat ${GO_RR3}/*/GENE_TREES_*_UFBoot_MAX-*.newick > ${GO_RR3}/GENE_TREES_ALL_REPLICATES.newick
		 


######## POLAR 2: Reference sequence: B. nana
#echo "Polarization 2..."
	  
# polarization label
NALT=2
	 
GO_AA3=${AA3}/ALT${NALT}
GO_RR3=${RR3}/ALT${NALT}
rm -rf $GO_RR3
	 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)

   mkdir -p $GO_RR3/$r_aux

   ########output file
   OUTfile=GENE_TREES_${r_aux}_UFBoot_MAX-${N_GENES}.newick
   rm -f ${GO_RR3}/$r_aux/$OUTfile
		    
   for i in `seq 1 1 $N_GENES`; do                 # cycle through genes
		   
               #input file
               GO_GENE=${GENE_LIST[$i]}
               InFile=data_${GO_GENE}_${InFileSUF}
               #InFile=data_${GO_GENE}_TRUE_CLEAN_UFBoot_${MODEL_outfilename}.iqtree
		
               # get IQ-TEST2 consensus tree for each gene           
               cat $GO_AA3/$r_aux/${InFile} | grep -A 2 'Consensus tree in newick format:' | tail -1 >> ${GO_RR3}/$r_aux/$OUTfile
   done
done
		
# get genes from all replicates
cat ${GO_RR3}/*/GENE_TREES_*_UFBoot_MAX-*.newick > ${GO_RR3}/GENE_TREES_ALL_REPLICATES.newick
	 
		   
		   
######## POLAR 3: Reference sequence: B. humilis
#echo "Polarization 3..."
	  
# polarization label
NALT=3
	 
GO_AA3=${AA3}/ALT${NALT}
GO_RR3=${RR3}/ALT${NALT}
rm -rf $GO_RR3
	 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)

   mkdir -p $GO_RR3/$r_aux

   ########output file
   OUTfile=GENE_TREES_${r_aux}_UFBoot_MAX-${N_GENES}.newick
   rm -f ${GO_RR3}/$r_aux/$OUTfile
		      
   for i in `seq 1 1 $N_GENES`; do                 # cycle through genes
		   
               #input file
               GO_GENE=${GENE_LIST[$i]}
               InFile=data_${GO_GENE}_${InFileSUF}
               #InFile=data_${GO_GENE}_TRUE_CLEAN_UFBoot_${MODEL_outfilename}.iqtree
		
               # get IQ-TEST2 consensus tree for each gene           
               cat $GO_AA3/$r_aux/${InFile} | grep -A 2 'Consensus tree in newick format:' | tail -1 >> ${GO_RR3}/$r_aux/$OUTfile
   done
done
		
# get genes from all replicates
cat ${GO_RR3}/*/GENE_TREES_*_UFBoot_MAX-*.newick > ${GO_RR3}/GENE_TREES_ALL_REPLICATES.newick
	   
###################
################### For each gene family, ID species closer to polyploid
###################	
echo
echo "4. Pairing frequencies..." 
date -u

AA4=${RR_MAIN}/03_PAIRINGfrequencies
	  
## tetraploid IDs
TETRAPLOID="TETRAPLOID_SIM"
	  
## outgroup
outgroup_ID="alnus_SIM"
		 

######## POLAR 1: Reference sequence: B. pendula
#echo "Polarization 1..."
	  
# polarization label
NALT=1
		  
GO_AA4=${AA4}/ALT${NALT}
GO_RR4=$GO_AA4
		  
## output file
OUTFile="sister_ID_analysis_priors-ALT${NALT}.txt"
		 
Rscript --no-save $SRCDIR_INI/24_IQTREE_gene_tree_DISTANCE_ALT.R \
				                                     ${GO_AA4}/GENE_TREES_ALL_REPLICATES.newick \
                                                                     $TETRAPLOID \
                                                                     $outgroup_ID \
                                                                     $GO_RR4 \
                                                                     $OUTFile
		  
$SRCDIR_INI/23_pairing_profile_summary.sh $GO_RR4 $OUTFile
		  

######## POLAR 2: Reference sequence: B. nana
#echo "Polarization 2..."
		  
# polarization label
NALT=2
		  
GO_AA4=${AA4}/ALT${NALT}
GO_RR4=$GO_AA4
		  
## output file
OUTFile="sister_ID_analysis_priors-ALT${NALT}.txt"
		 
Rscript --no-save $SRCDIR_INI/24_IQTREE_gene_tree_DISTANCE_ALT.R \
				                                     ${GO_AA4}/GENE_TREES_ALL_REPLICATES.newick \
                                                                     $TETRAPLOID \
                                                                     $outgroup_ID \
                                                                     $GO_RR4 \
                                                                     $OUTFile
			 
$SRCDIR_INI/23_pairing_profile_summary.sh $GO_RR4 $OUTFile
		  
		  
######## POLAR 3: Reference sequence: B. humilis
#echo "Polarization 3..."
		  
# polarization label
NALT=3
		  
GO_AA4=${AA4}/ALT${NALT}
GO_RR4=$GO_AA4
		  
## output file
OUTFile="sister_ID_analysis_priors-ALT${NALT}.txt"
		 
Rscript --no-save $SRCDIR_INI/24_IQTREE_gene_tree_DISTANCE_ALT.R \
				                                     ${GO_AA4}/GENE_TREES_ALL_REPLICATES.newick \
                                                                     $TETRAPLOID \
                                                                     $outgroup_ID \
                                                                     $GO_RR4 \
                                                                     $OUTFile
			 
$SRCDIR_INI/23_pairing_profile_summary.sh $GO_RR4 $OUTFile
		  

####### Joint pairing profile (all polarization geometries)

cat $AA4/ALT1/sister_ID_analysis_priors-ALT1_summary.txt $AA4/ALT2/sister_ID_analysis_priors-ALT2_summary.txt $AA4/ALT3/sister_ID_analysis_priors-ALT3_summary.txt | grep -v '^p'> $AA4/sister_ID_analysis_priors-JOINT_summary.txt		  



###### Clean scratch disk
rm -f $SNIC_TMP/_*
rm -f $SNIC_TMP/data*
rm -f $SNIC_TMP/*.py

rm -f $SNIC_TMP/*_UFBoot_*
rm -f $SNIC_TMP/*fa
rm -f $SNIC_TMP/*fasta
rm -f $SNIC_TMP/AUX*
rm -rf $SNIC_TMP/ALT*
rm -rf $SNIC_TMP/0*

######################################################### END
#echo
#echo "-------------20_POLARIZATION_IQTREE_PAIRINGprofiles_MAIN.sh done-------------"
#date -u






