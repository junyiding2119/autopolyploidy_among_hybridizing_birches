#!/bin/bash
#SBATCH -J fsc_loopy
#SBATCH -n 12
#SBATCH -t 5:00:00
#SBATCH -A uppmax2025-2-230
#SBATCH -M pelle

#pwd:/proj/snic2020-6-184/private/DingJY/data_analysis/costataeEvolution/Result/08fastsimcoal/goodness-of-fit/model02_maxL.bootstrap

fsc=/proj/snic2020-6-184/private/DingJY/software/fsc27_linux64/fsc27093

$fsc -i model02_maxL.bootstrap.par -n 5000 -j -q -c 12 -B 12

# -j  --jobs              : output one simulated or bootstrapped SFS per file
#                           in a separate directory for easier analysis
#                           which means simulate one file each time

# -q Quiet Mode

# -n  --numsims 1000      : number of simulations to perform
#                           Also applies for parameter estimation


