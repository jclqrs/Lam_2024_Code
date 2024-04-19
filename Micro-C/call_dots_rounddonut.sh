#!/bin/bash
#
#SBATCH -J call_rounddots
#SBATCH -o call_rounddots-%j.out
#SBATCH -n 16
#SBATCH --mem-per-cpu=32G

date
conda init bash
source ~/miniconda3/etc/profile.d/conda.sh

conda activate py38

#MCLR = mcool file
set -x

while getopts m: option
do 
	case "${option}"
		in
		m) MCLR=${OPTARG};;
	esac
done

CHROMFILE=mm9_chr1toX.bed
OUTPREFIX=$(basename ${MCLR} | cut -d . -f 1)


RES="::/resolutions/5000"
CLR=${MCLR}${RES}
OUT=${OUTPREFIX}_5k

# NAME OUTPUT FILES
EXP=${OUT}_expected-cis.tsv
OUT_LOOPS=${OUT}_loops

# Only make expected if it does not exist in same folder
if [ ! -f "$EXP" ]; then
	# MAKE EXPECTED 
	cooltools expected-cis -p 16 -o $EXP \
		--view ${CHROMFILE}  \
		--clr-weight-name weight \
		$CLR  
fi

python round_dot_caller.py -i ${MCLR} -e ${EXP} -r 5000

# 10K
RES="::/resolutions/10000"
CLR=${MCLR}${RES}
OUT=${OUTPREFIX}_10k

# NAME OUTPUT FILES
EXP=${OUT}_expected-cis.tsv
OUT_LOOPS=${OUT}_loops

# Only make expected if it does not exist in same folder
if [ ! -f "$EXP" ]; then
	# MAKE EXPECTED 
	cooltools expected-cis -p 16 -o $EXP \
		--view ${CHROMFILE}  \
		--clr-weight-name weight \
		$CLR  
fi

python round_dot_caller.py -i ${MCLR} -e ${EXP} -r 10000
