#!/bin/bash
#
#SBATCH -J fsc_loopy
#SBATCH -n 6
#SBATCH -t 24:00:00
#SBATCH -A uppmax2025-2-230
#SBATCH -M pelle

##### load modules
#module load bioinfo-tools

# local installation of fastSimcoal2 (ver 2.7.0.9)
fsc=/proj/snic2020-6-184/private/DingJY/software/fsc27_linux64/fsc27093

PREFIX=$1 #model name
repstart=$2
repend=$3

DIR_INI=$(pwd)
DIR_INPUT=$DIR_INI
FILE_INPUT1_0=${DIR_INPUT}/6Pop_jointDAFpop1_0.obs
FILE_INPUT2_0=${DIR_INPUT}/6Pop_jointDAFpop2_0.obs
FILE_INPUT2_1=${DIR_INPUT}/6Pop_jointDAFpop2_1.obs
FILE_INPUT3_0=${DIR_INPUT}/6Pop_jointDAFpop3_0.obs
FILE_INPUT3_1=${DIR_INPUT}/6Pop_jointDAFpop3_1.obs
FILE_INPUT3_2=${DIR_INPUT}/6Pop_jointDAFpop3_2.obs
FILE_INPUT4_0=${DIR_INPUT}/6Pop_jointDAFpop4_0.obs
FILE_INPUT4_1=${DIR_INPUT}/6Pop_jointDAFpop4_1.obs
FILE_INPUT4_2=${DIR_INPUT}/6Pop_jointDAFpop4_2.obs
FILE_INPUT4_3=${DIR_INPUT}/6Pop_jointDAFpop4_3.obs
FILE_INPUT5_0=${DIR_INPUT}/6Pop_jointDAFpop5_0.obs
FILE_INPUT5_1=${DIR_INPUT}/6Pop_jointDAFpop5_1.obs
FILE_INPUT5_2=${DIR_INPUT}/6Pop_jointDAFpop5_2.obs
FILE_INPUT5_3=${DIR_INPUT}/6Pop_jointDAFpop5_3.obs
FILE_INPUT5_4=${DIR_INPUT}/6Pop_jointDAFpop5_4.obs


tpl_Input=${DIR_INPUT}/${PREFIX}/${PREFIX}.tpl
est_Input=${DIR_INPUT}/${PREFIX}/${PREFIX}.est

dos2unix $tpl_Input $est_Input
echo
echo "Model: " $PREFIX
echo "Start rep: " $repstart
echo "End rep: " $repend
echo "Polyploidy species: " $species
echo $tpl_Input
echo $est_Input
echo
echo
echo "Starting fsc jobs ..."
date -u
echo

###### Main
if [ ! -d "$PREFIX" ]; then
mkdir -p $PREFIX
fi

cd $PREFIX

if [ ! -f "${PREFIX}_jointDAFpop1_0.obs" ]; then
cp $FILE_INPUT1_0 ${PREFIX}_jointDAFpop1_0.obs
cp $FILE_INPUT2_0 ${PREFIX}_jointDAFpop2_0.obs
cp $FILE_INPUT2_1 ${PREFIX}_jointDAFpop2_1.obs
cp $FILE_INPUT3_0 ${PREFIX}_jointDAFpop3_0.obs
cp $FILE_INPUT3_1 ${PREFIX}_jointDAFpop3_1.obs
cp $FILE_INPUT3_2 ${PREFIX}_jointDAFpop3_2.obs
cp $FILE_INPUT4_0 ${PREFIX}_jointDAFpop4_0.obs
cp $FILE_INPUT4_1 ${PREFIX}_jointDAFpop4_1.obs
cp $FILE_INPUT4_2 ${PREFIX}_jointDAFpop4_2.obs
cp $FILE_INPUT4_3 ${PREFIX}_jointDAFpop4_3.obs
cp $FILE_INPUT5_0 ${PREFIX}_jointDAFpop5_0.obs
cp $FILE_INPUT5_1 ${PREFIX}_jointDAFpop5_1.obs
cp $FILE_INPUT5_2 ${PREFIX}_jointDAFpop5_2.obs
cp $FILE_INPUT5_3 ${PREFIX}_jointDAFpop5_3.obs
cp $FILE_INPUT5_4 ${PREFIX}_jointDAFpop5_4.obs
echo no obs
else
	echo obs exit
fi

 for i in `seq $repstart $repend`
 do
   mkdir run${i}
   cp $tpl_Input $est_Input ${PREFIX}_jointDAFpop1_0.obs ${PREFIX}_jointDAFpop2_0.obs ${PREFIX}_jointDAFpop2_1.obs ${PREFIX}_jointDAFpop3_0.obs ${PREFIX}_jointDAFpop3_1.obs ${PREFIX}_jointDAFpop3_2.obs ${PREFIX}_jointDAFpop4_0.obs ${PREFIX}_jointDAFpop4_1.obs ${PREFIX}_jointDAFpop4_2.obs ${PREFIX}_jointDAFpop4_3.obs ${PREFIX}_jointDAFpop5_0.obs ${PREFIX}_jointDAFpop5_1.obs ${PREFIX}_jointDAFpop5_2.obs ${PREFIX}_jointDAFpop5_3.obs ${PREFIX}_jointDAFpop5_4.obs ./run$i"/"
   cd run${i}
   $fsc -t ${PREFIX}.tpl -e ${PREFIX}.est -n100000 -M -d -0 -l 10 -L 40 -q -c 6 -B 6
   cd ..

   echo "$PREFIX run${i} finished"
 done

echo
echo "-------------fsc has finished-------------"
date -u
echo
#${DIR_INPUT}/fsc-selectbestrun.sh ${PREFIX}
