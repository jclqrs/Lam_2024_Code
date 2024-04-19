#!/bin/bash
#
#SBATCH -J chip_pipeline
#SBATCH -n 64
#SBATCH -o log_chip_pipeline_%j.txt
#SBATCH -t 48:00:00

module load SAMtools

date

conda init bash 
source ~/miniconda3/etc/profile.d/conda.sh

conda activate py38

while getopts c: option
do 
	case "${option}"
		in
		c) config=${OPTARG};;
	esac
done

mkdir ./bam
mkdir ./bigwig
mkdir ./macs2results

# blacklist file
blist=/home/lamj2/data/ref/blacklist.bed

echo "PROCESS INPUTS"
##################################################################
###########                           					##########
###########         PROCESS INPUT            			##########                   
###########                           					##########                   
###########                           					##########                   
##################################################################
                                  

# parallelize alignment to speed things up

num_jobs="\j"  # The prompt escape for number of jobs currently running

for fq in $(ls ./input/*R1*)
do
	while (( ${num_jobs@P} >= 8 )); do
	wait -n
	done
	echo $fq
	# get sample name
	ID=$(basename ${fq} | cut -f1 -d_)
	sample=$(basename ${fq} | cut -f2 -d_)
	fq1=$(ls  ./input/ | grep ${ID} | grep _R1_)
	fq2=$(ls ./input/ | grep ${ID} | grep _R2_)
	echo ${ID}
	echo $fq1
	echo $fq2


	bowtie2 -p 8 --local \
		-x ~/data/ref/mm9_bt2/mm9 \
		-1  ./input/$fq1 -2 ./input/$fq2 \
		-S ~/scr/chip_input/bam/${sample}.sam &
done
wait

for fq in $(ls ./input/*R1*)
do
# get sample name
sample=$(basename ${fq} | cut -f2 -d_)


# remove duplicates
samtools collate -o ~/scr/chip_input/bam/${sample}_collate.bam ~/scr/chip_input/bam//${sample}.sam -@16
samtools fixmate -m ~/scr/chip_input/bam/${sample}_collate.bam ~/scr/chip_input/bam/${sample}_fixmate.bam -@16
samtools sort -o ~/scr/chip_input/bam/${sample}_positionsort.bam ~/scr/chip_input/bam/${sample}_fixmate.bam -@16
samtools markdup -r -s -f ~/scr/chip_input/bam/${sample}_stats.txt \
	~/scr/chip_input/bam/${sample}_positionsort.bam ~/scr/chip_input/bam/${sample}_nodups.bam -@16

# filter MAPQ < 20
samtools view -bhSq 20 ~/scr/chip_input/bam/${sample}_nodups.bam -@16 > ~/scr/chip_input/bam/${sample}_filtered.bam

# index
samtools index ~/scr/chip_input/bam/${sample}_filtered.bam -@16

# make rpkm normalized bigwig
bamCoverage --bam ~/scr/chip_input/bam/${sample}_filtered.bam \
	-o ~/scr/chip_input/bigwig/${sample}_normed_rpm.bw \
	--binSize 20 \
	--smoothLength 60 \
	--normalizeUsing BPM \
	--ignoreDuplicates \
	--centerReads \
	--minMappingQuality 20 \
	--extendReads 300 \
	-bl ${blist} \
	--ignoreForNormalization chrX \
	-p 16 \
	--verbose

done
wait

# move stats files to their own folder
mkdir ~/scr/chip_input/stats
mv ~/scr/chip_input/bam/*stats.txt ~/scr/chip_input/stats




##################################################################
###########                           					##########
###########         PROCESS SAMPLES            			##########                   
###########                           					##########                   
###########                           					##########                   
##################################################################


for fq in $(ls ./fastq/*R1*)
do
	while (( ${num_jobs@P} >= 8 )); do
	wait -n
	done
	echo $fq
	# get sample name
	ID=$(basename ${fq} | cut -f1 -d_)
	sample=$(grep ${ID} ${config} | cut -f2 -d,)
	input=$(grep ${ID} ${config} | cut -f3 -d,)
	fq1=$(ls  ./fastq/ | grep ${ID} | grep _R1_)
	fq2=$(ls ./fastq/ | grep ${ID} | grep _R2_)
	echo ${ID}
	echo $fq1
	echo $fq2


	bowtie2 -p 8 --local \
		-x ~/data/ref/mm9_bt2/mm9 \
		-1  ./fastq/$fq1 -2 ./fastq/$fq2 \
		-S ./bam/${sample}.sam &
done
wait

for fq in $(ls ./fastq/*R1*)
do

echo $fq
# get sample name
ID=$(basename ${fq} | cut -f1 -d_)
sample=$(grep ${ID} ${config} | cut -f2 -d,)
input=$(grep ${ID} ${config} | cut -f3 -d,)
fq1=$(ls  ./fastq/ | grep ${ID} | grep _R1_)
fq2=$(ls ./fastq/ | grep ${ID} | grep _R2_)
echo ${ID}
echo $fq1
echo $fq2

#remove duplciates
samtools collate -o ./bam/${sample}_collate.bam ./bam/${sample}.sam -@16
samtools fixmate -m ./bam/${sample}_collate.bam ./bam/${sample}_fixmate.bam -@16
samtools sort -o ./bam/${sample}_sort.bam ./bam/${sample}_fixmate.bam -@16
samtools markdup -r -s -f ./bam/${sample}_stats.txt \
	./bam/${sample}_sort.bam ./bam/${sample}_nodups.bam -@16

# filter MAPQ < 20
samtools view -bhSq 20 ./bam/${sample}_nodups.bam -@16 > ./bam/${sample}_filtered.bam

# index
samtools index ./bam/${sample}_filtered.bam -@16


# link sample with control bam using list in config csv
#input_bam=${input}_sorted.bam

# Macs2
macs2 callpeak -t ./bam/${sample}_filtered.bam \
	-c ~/scr/chip_input/bam/${input}_filtered.bam \
	-g 1.87e9 \
	-n ${sample}_q05 \
	-f BAMPE \
	--outdir macs2results

# make rpkm normalized bigwig
bamCoverage --bam ./bam/${sample}_filtered.bam \
	-o ./bigwig/${sample}_normed_rpm.bw \
	--binSize 20 \
	--smoothLength 60 \
	--normalizeUsing BPM \
	--ignoreDuplicates \
	--centerReads \
	--minMappingQuality 20 \
	--extendReads 300 \
	-bl ${blist} \
	--ignoreForNormalization chrX \
	-p 16 \
	--verbose

done

rm ./bam/*_nodups.bam
rm ./bam/*_sort
rm ./bam/*_fixmate.bam
rm ./bam/*_collate.bam
rm ./bam/*.sam