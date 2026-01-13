#!/bin/bash
#


## Glean gene tree counts from outputs produced by 24_IQTREE_gene_tree_DISTANCE_ALT.R


#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

# input folder 
GO_AA=${1:?msg}

# input file
In_file=${2:?msg}




#echo
GO_RR=${GO_AA}

## output file
outFile=${In_file%.txt}_summary.txt
rm -f $GO_RR/$outFile
echo "pendula" "platyphylla" "pend_platy" "populifolia" "pop_out" "nana" "nana_out" "occidentalis" "occ_out" "cos" "ash" "cos_ash" "bug" "bug_out" "mca" "mca_out" "len" "basal" > $GO_RR/$outFile


# get counts (%)
cd $GO_AA
TC_PEND="$(grep -w 'pendula_SIM' ${In_file} | cut -f3)"
TC_PLATY="$(grep -w 'platyphylla_SIM' ${In_file} | cut -f3)"
TC_PEND_PLATY="$(grep -w 'pendula_SIM_platyphylla_SIM' ${In_file} | cut -f3)"
TC_POP="$(grep -w 'populifolia_SIM' ${In_file} | cut -f3)"
TC_POP_OUT="$(grep -w 'pendula_SIM_platyphylla_SIM_populifolia_SIM' ${In_file} | cut -f3)"
TC_NANA="$(grep -w 'nana_SIM' ${In_file} | cut -f3)"
TC_NANA_OUT="$(grep -w 'nana_SIM_pendula_SIM_platyphylla_SIM_populifolia_SIM' ${In_file} | cut -f3)"
TC_OCC="$(grep -w 'occidentalis_SIM' ${In_file} | cut -f3)"
TC_OCC_OUT="$(grep -w 'nana_SIM_occidentalis_SIM_pendula_SIM_platyphylla_SIM_populifolia_SIM' ${In_file} | cut -f3)"
TC_COS="$(grep -w 'costata_SIM' ${In_file} | cut -f3)"
TC_ASH="$(grep -w 'ashburneri_SIM' ${In_file} | cut -f3)"
TC_COS_ASH="$(grep -w 'ashburneri_SIM_costata_SIM' ${In_file} | cut -f3)"
TC_BUG="$(grep -w 'buggsii_SIM' ${In_file} | cut -f3)"
TC_BUG_OUT="$(grep -w 'ashburneri_SIM_buggsii_SIM_costata_SIM' ${In_file} | cut -f3)"
TC_MCA="$(grep -w 'mcallisteri_SIM' ${In_file} | cut -f3)"
TC_MCA_OUT="$(grep -w 'ashburneri_SIM_buggsii_SIM_costata_SIM_mcallisteri_SIM' ${In_file} | cut -f3)"
TC_LEN="$(grep -w 'lenta_SIM' ${In_file} | cut -f3)"
TC_BASAL="$(grep -w 'alnus_SIM_ashburneri_SIM_buggsii_SIM_costata_SIM_lenta_SIM_mcallisteri_SIM_nana_SIM_occidentalis_SIM_pendula_SIM_platyphylla_SIM_populifolia_SIM' ${In_file} | cut -f3)"

if [[ $TC_PEND = "" ]]; then TC_PEND=0; fi
if [[ $TC_PLATY = "" ]]; then TC_PLATY=0; fi
if [[ $TC_PEND_PLATY = "" ]]; then TC_PEND_PLATY=0; fi
if [[ $TC_POP = "" ]]; then TC_POP=0; fi
if [[ $TC_POP_OUT = "" ]]; then TC_POP_OUT=0; fi
if [[ $TC_NANA = "" ]]; then TC_NANA=0; fi
if [[ $TC_NANA_OUT = "" ]]; then TC_NANA_OUT=0; fi
if [[ $TC_OCC = "" ]]; then TC_OCC=0; fi
if [[ $TC_OCC_OUT = "" ]]; then TC_OCC_OUT=0; fi
if [[ $TC_COS = "" ]]; then TC_COS=0; fi
if [[ $TC_ASH = "" ]]; then TC_ASH=0; fi
if [[ $TC_COS_ASH = "" ]]; then TC_COS_ASH=0; fi
if [[ $TC_BUG = "" ]]; then TC_BUG=0; fi
if [[ $TC_BUG_OUT = "" ]]; then TC_BUG_OUT=0; fi
if [[ $TC_MCA = "" ]]; then TC_MCA=0; fi
if [[ $TC_MCA_OUT = "" ]]; then TC_MCA_OUT=0; fi
if [[ $TC_LEN = "" ]]; then TC_LEN=0; fi
if [[ $TC_BASAL = "" ]]; then TC_BASAL=0; fi

echo $TC_PEND $TC_PLATY $TC_PEND_PLATY $TC_POP $TC_POP_OUT $TC_NANA $TC_NANA_OUT $TC_OCC $TC_OCC_OUT $TC_COS $TC_ASH $TC_COS_ASH $TC_BUG $TC_BUG_OUT $TC_MCA $TC_MCA_OUT $TC_LEN $TC_BASAL>> $GO_RR/$outFile
