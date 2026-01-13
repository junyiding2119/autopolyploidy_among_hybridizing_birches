#!/bin/bash

Species=ermanii
Ref=ash #ash,bugg,cos,platy
for Ind in {}
do
In_file=sister_ID_${Ind}Iteration_${Ref}.txt
./02s_Summary_pairing_matrix_query.sh $In_file $Ind $Species $Ref

done