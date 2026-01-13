#!/bin/bash
#generate Summary pairing matrix
#### paths and folder names


SRCDIR_INI=$(pwd)

# input file
In_file=${1:?msg}

Ind=${2:?msg}

Species=${3:?msg}

Ref=${4:?msg}



GO_RR=${SRCDIR_INI}/${Species}/${Ind}

## output file
outFile=pairing_matrix_${Ind}_${Ref}.txt
rm -f $GO_RR/$outFile
echo "pendula" "platyphylla" "pend_platy" "populifolia" "pop_out" "nana" "nana_out" "occidentalis" "occ_out" "cos" "ash_EDL" "ash_TZZ" "ash_EDL_TZZ" "cos_ash" "bug" "bug_out" "mca" "mca_out" "len" "basal" > $GO_RR/$outFile


# get counts (%)
cd $GO_RR
TC_PEND="$(grep -w 'pendula_ERR2026097' ${In_file} | cut -f3)"
TC_PLATY="$(grep -w 'EDL21-009-platy' ${In_file} | cut -f3)"
TC_PEND_PLATY="$(grep -w 'EDL21-009-platy_pendula_ERR2026097' ${In_file} | cut -f3)"
TC_POP="$(grep -w 'populifolia_ERR2026269' ${In_file} | cut -f3)"
TC_POP_OUT="$(grep -w 'EDL21-009-platy_pendula_ERR2026097_populifolia_ERR2026269' ${In_file} | cut -f3)"
TC_NANA="$(grep -w 'nana_ERR2026268' ${In_file} | cut -f3)"
TC_NANA_OUT="$(grep -w 'EDL21-009-platy_nana_ERR2026268_pendula_ERR2026097_populifolia_ERR2026269' ${In_file} | cut -f3)"
TC_OCC="$(grep -w 'occidentalis_ERR2026271' ${In_file} | cut -f3)"
TC_OCC_OUT="$(grep -w 'EDL21-009-platy_nana_ERR2026268_occidentalis_ERR2026271_pendula_ERR2026097_populifolia_ERR2026269' ${In_file} | cut -f3)"
TC_COS="$(grep -w 'QSLC10' ${In_file} | cut -f3)"
TC_ASH_EDL="$(grep -w 'EDL028' ${In_file} | cut -f3)"
TC_ASH_TZZ="$(grep -w 'TZZ78' ${In_file} | cut -f3)"
TC_ASH_EDL_TZZ="$(grep -w 'EDL028_TZZ78' ${In_file} | cut -f3)"
TC_COS_ASH="$(grep -w 'EDL028_QSLC10_TZZ78' ${In_file} | cut -f3)"
TC_BUG="$(grep -w 'CKX18-006' ${In_file} | cut -f3)"
TC_BUG_OUT="$(grep -w 'CKX18-006_EDL028_QSLC10_TZZ78' ${In_file} | cut -f3)"
TC_MCA="$(grep -w 'XYB21-001-mcalli' ${In_file} | cut -f3)"
TC_MCA_OUT="$(grep -w 'CKX18-006_EDL028_QSLC10_TZZ78_XYB21-001-mcalli' ${In_file} | cut -f3)"
TC_LEN="$(grep -w 'lenta_ERR2026272' ${In_file} | cut -f3)"
TC_BASAL="$(grep -w 'Alnus_glutinosa_ERR2103731_CKX18-006_EDL028_EDL21-009-platy_lenta_ERR2026272_nana_ERR2026268_occidentalis_ERR2026271_pendula_ERR2026097_populifolia_ERR2026269_QSLC10_TZZ78_XYB21-001-mcalli' ${In_file} | cut -f3)"

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

echo $TC_PEND $TC_PLATY $TC_PEND_PLATY $TC_POP $TC_POP_OUT $TC_NANA $TC_NANA_OUT $TC_OCC $TC_OCC_OUT $TC_COS $TC_ASH_EDL $TC_ASH_TZZ $TC_ASH_EDL_TZZ  $TC_COS_ASH $TC_BUG $TC_BUG_OUT $TC_MCA $TC_MCA_OUT $TC_LEN $TC_BASAL>> $GO_RR/$outFile