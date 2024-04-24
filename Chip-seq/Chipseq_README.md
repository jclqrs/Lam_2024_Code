# Walkthrough for ChIP-seq analysis

# Set up python environment
The `spec-file.txt` file can be used with conda to install the software required to run the loop quantification script on a Linux platform.

`conda create --name <env> --file <this file>`

For other platforms, you can manually install the packages listed in the file (and in the Methods section). Importantly, the loop quantification script uses cooltools version 0.5.1. 

# ChIP-seq pipeline
`sbatch chip_pipeline_2.sh -c CONFIG.csv`

The `chip_pipeline_2.sh` script includes steps for aligning and de-duplicating samples, calling peaks with macs2, and creating bigwigs with deeptools. 

Before running the pipeline, you have to download paired-end sample fastq files to a directory named `fastq` and input fastq files to a directory named `input`

You also need to create a CONFIG file where column 1 contains unique sample ID (which must be present in the fastq filename), column2 contains the sample names, and column 3 contains the name of the input file that should be used for peak calling. An example is provide. 

# Pol II ChIP-seq raw counts
For Pol II ChIP-seq, the script `fetch_gb_counts.sh` was used to quantify reads at gene bodies. A bed file containing gene body regions is provided (`longest_genebody.bed`).
