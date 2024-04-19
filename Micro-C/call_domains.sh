#!/bin/bash
#
#SBATCH -J calldomains
#SBATCH -n 8
#SBATCH -o calldomains_%j.txt
#SBATCH -t 48:00:00
#SBATCH --mem-per-cpu=64G


date

conda init bash 
source ~/miniconda3/etc/profile.d/conda.sh

conda activate py38

while getopts m: option
do 
	case "${option}"
		in
		m) mcool=${OPTARG};;
	esac
done

# get abs path of mcool 
mcool=$(readlink -f ${mcool})

outdir=$(basename ${mcool} | cut -f1 -d .)_domains
mkdir ${outdir}
cd ${outdir}


cool10k=${mcool}"::/resolutions/10000"

num_jobs="\j"  # The prompt escape for number of jobs currently running


# loop through chromosomes 
chrs=( $(seq 1 19))
chrs=("${chrs[@]}" "X")

# dump contact info into table for rgmap for each chromosome
for chrnum in ${chrs[@]}
do
	while (( ${num_jobs@P} >= 8 )); do
	wait -n
	done

	chrom=chr${chrnum}

	python fetch_contacts_for_rgmap.py -i ${cool10k} -c ${chrom} &
done
wait

conda activate Renv

for chrnum in ${chrs[@]}
do
	while (( ${num_jobs@P} >= 8 )); do
	wait -n
	done
	chrom=chr${chrnum}

	rgmap_domains.R ${chrom}_pixels.txt ${chrom} &
done
wait

