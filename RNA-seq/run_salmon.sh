#!/bin/bash
# 
#SBATCH -J run_salmon
#SBATCH -o run_salmon-%j.out
#SBATCH -n 8

IDX="/home/lamj2/data/ref/gencode_vM1_index/"

module load Salmon

mkdir quants

for fn in ./fastq/{3567..3572};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i $IDX -l A \
         -1 ${fn}*R1_001.fastq.gz \
         -2 ${fn}*R2_001.fastq.gz \
         -p 8 --validateMappings -o ./quants/${samp}_quant
done 