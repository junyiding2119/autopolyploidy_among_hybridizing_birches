#!/bin/bash
######## Reference species
ref=$1

######## Polyploidy sample ID
poly=$2

NSPEC=12

declare -a GENE_LIST

cd ./Iteration_${ref}

for LogFile in *log;do basename $LogFile .clean.fa-ALT.fasta.log;done > 00_gene_list.txt

i=1
while IFS= read -r LINE
do
   GENE_LIST[${i}]="${LINE}"
   ((i++))
done < 00_gene_list.txt

N_GENES=${#GENE_LIST[@]}


for i in `seq 1 1 $N_GENES`; do                 # cycle through genes
    
	GO_GENE=${GENE_LIST[$i]}
    #input file
    InFileSUF=${GO_GENE}.clean.fa-ALT.fasta.log

    # get number of informative sites
    NPS_AUX="$(cat ${InFileSUF} | grep 'parsimony-informative' | cut -f1 -d' ')"
    NPS[i]=$NPS_AUX

    if [[ "${NPS[$i]}" -ge 5 ]]
    then                  
       NPSF[i]=OK
    else
       NPSF[i]=FAIL
    fi

    #echo ${NPS[$i]} ${NPSF[$i]}

    # get number of sequences in alignment
    NSEQ[i]="$(cat ${InFileSUF} | grep 'Alignment has' | cut -f3 -d' ')"

    if [[ "${NSEQ[$i]}" = ${NSPEC} ]]
    then                  
       NSEQF[i]=OK
    else
       NSEQF[i]=FAIL
    fi

    # check whether there are sequences with identical nucleotide composition
    INC_AUX="$(grep 'is identical to' ${InFileSUF})"

    if [[ "${INC_AUX}" = "" ]]
    then                  
       INC[i]=OK
    else
       INC[i]=FAIL
    fi

    # check whether bs analysis converged
    BS_AUX="$(grep 'WARNING: bootstrap analysis did not converge. You should rerun with higher number of iterations (-nm option)' ${InFileSUF})"

    if [[ "${BS_AUX}" = "" ]]
    then                  
       BS[i]=OK
    else
       BS[i]=FAIL
    fi

    # check whether contains only gaps
    BS_AUX="$(grep 'contains only gaps or missing data' ${InFileSUF})"

    if [[ "${BS_AUX}" = "" ]]
    then                  
       OG[i]=OK
    else
       OG[i]=FAIL
    fi

done

# save to file
rm -f ../_AUX_${ref}
for i in `seq 1 1 $N_GENES`; do
  echo ${GENE_LIST[$i]} ${NPS[$i]} ${NPSF[$i]} ${INC[$i]} ${BS[$i]} ${NSEQ[$i]} ${NSEQF[$i]} ${OG[$i]} >> ../_AUX_${ref}
done

cd ..
