#!/bin/bash
#

### Script used to perform initial model optimization using simulated annealing.
### Script must be run independently for each polyploidization model, Species, and chain 
### There is no cross-talk between chains

###################################### paths and folder names

################## Remember initial path
SRCDIR_INI=$(pwd)   

################## Folder containing simulated (unpolarized) MSAs (root)
AA=/data_analysis/costataeEvolution/Result/04simphy_ILS

################## Simulation conditions (polyploidization models)
myMODEL[1]="AAAA"    # B. ashburneri, autopolyploid
myMODEL[2]="CCCC"    # B. costata, autopolyploid
myMODEL[3]="BBBB"    # B. buggsii, autopolyploid
myMODEL[4]="AacAacAacAac_A"    # ancestor of B. ashburneri-B. costata ,autopolyploid ,ashburneri gene flow
myMODEL[5]="AacAacAacAac_C"    # ancestor of B. ashburneri-B. costata ,costata,ashburneri gene flow

myMODEL[6]="AACC"    # B. ashburneri/B. costata, allopolyploid
myMODEL[7]="AABB"    # B. ashburneri/B. buggsii, allopolyploid
myMODEL[8]="CCBB"    # B. costata/B. buggsii, allopolyploid

################## Species
spe[1]="ermanii"
spe[2]="albosinensis"
spe[3]="utilis"

################## Homoeologous Exchange across subgenomes  (1: allowed;  0: not allowed) >> valid for allopolyploid models only
HomoeologousExchange=0

################## Number of chains
N_chains=4

################## Prior flag
#Priorflag="user"		# priors provided by user (PP file, below)
Priorflag="random"		# random generated priors

################## Number of replicates
REP=25		

################## Number of loci/genes		
NLOCUS=50

################## ILS prior
ILS_PRIOR=40

################## pendula-platyphylla gene flow (1: allowed;  0: not allowed)
geneFlowFlag=0







for k in `seq 2 1 2`; do               # models

   for n in `seq 1 1 1`; do               # species

      for j in `seq 1 1 $N_chains`; do               # chain
	  
         GO_myMODEL=${myMODEL[$k]}
         GO_spe=${spe[1]}
         chain_number=${j}
         SNIC_TMP=/data_analysis/costataeEvolution/Result/TMP/TMP_${GO_spe}_${GO_myMODEL}
		 if [ ! -d $SNIC_TMP ]; then
		     mkdir -p $SNIC_TMP
	     fi
		 
         ################## Observed pairing profiles (each line in input file correspond to a different polarization geometry, in the following order: B. pendula; B. nana; B. humilis; B. platyphylla)
         PP_Obs=/data_analysis/costataeEvolution/Result/05ABC_Simulated_annealing/00_PAIRING_PROFILES_${GO_spe}.txt	

         ##################  Output folder (root); chain number is the last digit in the path

         if [[ ${HomoeologousExchange} = 0 ]] ; then
            RR_MAIN=/data_analysis/costataeEvolution/Result/05ABC_Simulated_annealing/${GO_myMODEL}/${GO_spe}/${chain_number}
            echo ""
         fi

         if [[ ${HomoeologousExchange} = 1 ]] ; then
            RR_MAIN=/data_analysis/costataeEvolution/Result/05ABC_Simulated_annealing/${GO_myMODEL}_HE/${GO_spe}/${chain_number}
            echo ""
         fi

         rm -rf $RR_MAIN
         mkdir -p $RR_MAIN
		 

         ################## Initial priors (when provided by the user)
         #PP=/crex1/proj/snic2017-7-149/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/PRIORS_simAnnealing_userDefined/00_PRIORS_ILS_gf_H1_H7_SIMANN_${GO_myMODEL}_${GO_spe}.txt    		
         PP=""

         if [[ ${PP} = "" ]] ; then
            PP=$SRCDIR_INI/PPempty.txt
            if [[ -e "$PP" ]]; then
               echo $PP "file exists."
            else
               touch $PP
            fi
         fi


         echo "MODEL:" $GO_myMODEL
         echo "Species:" $GO_spe
         echo "Chain:" $chain_number
         echo "Input folder:" $AA
         echo "Output folder:" $RR_MAIN
         echo "Observed pairing profiles:" $PP_Obs
         echo "Number of replicates:" $REP
         echo "Number of loci:" $NLOCUS
         echo "ILS prior:" $ILS_PRIOR
         if [[ ${Priorflag} = "random" ]] ; then
            echo "Initial priors: random generated"
         fi
         if [[ ${Priorflag} = "user" ]] ; then
            echo "Initial priors:" $PP
         fi
         echo "HomoeologousExchange flag:" $HomoeologousExchange
         echo "pendula-platyphylla gene flow flag:" $geneFlowFlag
         echo "TMP dir:" $SNIC_TMP


         $SRCDIR_INI/05_ABC_SimAnnealing_query.sh $AA \
                                                  $RR_MAIN \
                                                  $PP_Obs \
                                                  $REP \
                                                  $NLOCUS \
                                                  $ILS_PRIOR \
                                                  $GO_myMODEL \
                                                  $HomoeologousExchange \
                                                  $geneFlowFlag \
                                                  $Priorflag \
                                                  $PP \
												  $SNIC_TMP
                                            
											
      done
   done
done
