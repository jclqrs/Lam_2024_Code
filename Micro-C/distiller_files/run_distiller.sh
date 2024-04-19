#!/bin/bash
#
#SBATCH -J pipeline_distiller
#SBATCH -o log_distiller-%j.out
#SBATCH -t 2-0:00 

# Clone distiller pipeline with command below
# before running this script. 
nextflow clone open2c/distiller-nf ./


module load singularity
nextflow run distiller.nf -params-file project.yml -profile cluster

