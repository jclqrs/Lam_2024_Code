#!/bin/bash
#
#SBATCH -J chip_pipeline
#SBATCH -n 1
#SBATCH --mem-per-cpu=16G  
#SBATCH -o log_htseqcount_%j.txt
#SBATCH -t 48:00:00


module load petagene


SAMPLE_REF="../mm9.refGene.gtf"
OUT='../YY1aid_ttseq_sample_counts.txt'


htseq-count -f bam \
			-r pos \
			-t transcript \
			-i gene_id \
			-s reverse \
			$(ls | grep "\.bam") \
			${SAMPLE_REF} > ${OUT}