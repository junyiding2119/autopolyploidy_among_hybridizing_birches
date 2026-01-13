#!/bin/bash
#SBATCH --job-name=arlsum
#SBATCH --array=2-10000%50
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=logs/arlsum_%a.out


#pwd:/proj/snic2020-6-184/private/DingJY/data_analysis/costataeEvolution/Result/08fastsimcoal/goodness-of-fit/model02_maxL.bootstrap

arlsumstat=/proj/snic2020-6-184/private/DingJY/data_analysis/costataeEvolution/Result/08fastsimcoal/goodness-of-fit/arlsumstat_linux/arlsumstat3522_64bit


SIM_DIR="model02_maxL.bootstrap"
OUT_DIR="stats_output"

# get .arp file name
# file name will be model02_maxL_1_1.arp, model02_maxL_1_2.arp, ...
ARP_FILE="${SIM_DIR}/model02_maxL.bootstrap_1_${SLURM_ARRAY_TASK_ID}.arp"

if [ ! -d "$OUT_DIR" ]; then
mkdir -p ${OUT_DIR}
fi

if [[ ! -f "$ARP_FILE" ]]; then
    echo "ERROR: File not found: $ARP_FILE"
    exit 1
fi

# no head
$arlsumstat "$ARP_FILE" "${OUT_DIR}/stats_${SLURM_ARRAY_TASK_ID}.txt" 0 0

echo "Done: $ARP_FILE"