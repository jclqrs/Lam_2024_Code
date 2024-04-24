#!/bin/bash
#
#SBATCH -J pol2_gb
#SBATCH -n 16
#SBATCH -o log_pol2_gb_%j.txt
#SBATCH -t 48:00:00

module load petagene/current

multiBamSummary BED-file \
     --bamfiles $(ls $PWD | grep "_filtered.bam$")\
     --BED ../longest_genebody.bed\
     --smartLabels \
     -out ../longest_genebody.npz --outRawCounts ../longest_genebody.tab \
     -v \
     -p 16

