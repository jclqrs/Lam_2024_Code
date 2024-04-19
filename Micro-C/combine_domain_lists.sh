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

while getopts m:e:t:n:f:u:o: option
do 
	case "${option}"
		in
		m) cool1=${OPTARG};;
		e) exp1=${OPTARG};;
		t) tad1=${OPTARG};;
		n) cool2=${OPTARG};;
		f) exp2=${OPTARG};;
		u) tad2=${OPTARG};;
		o) out=${OPTARG};;
	esac
done

domain_mergefilter.py -m ${cool1} -e ${exp1} -t ${tad1} -n ${cool2} -f ${exp2} -u ${tad2} -o ${out}