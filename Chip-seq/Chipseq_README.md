# Walkthrough for ChIP-seq analysis

# ChIP-seq pipeline for alignment, peak calling, and bigwig generation
`sbatch chip_pipeline_2.sh -c CONFIG.csv`

The `chip_pipeline_2.sh` script includes steps for aligning and de-duplicating samples, calling peaks with macs2, and creating bigwigs with deeptools. 

Before running the pipeline, you have to download paired-end sample fastq files to a directory named `fastq` and input fastq files to a directory named `input`

You also need to create a CONFIG file where column 1 contains unique sample ID (which must be present in the fastq filename), column2 contains the sample names, and column 3 contains the name of the input file that should be used for peak calling. An example config file is provided: `example_config_h3k27ac_chip.csv`. 

# Combining peak calls across replicates
We used Diffbind to merge peak calls and standardize peak sizes across replicates. We also used Diffbind to quantify counts at merged, standardized peaks. Example code to do the above steps is found in the Jupyter notebook `chipseq_diffbind_example.ipynb`. The file `config_diffbind_q05.csv` is an example of a configuration file that can be used with Diffbind.

# Pol II ChIP-seq raw counts
For Pol II ChIP-seq, the script `fetch_gb_counts.sh` was used to quantify reads at gene bodies. A bed file containing gene body regions is provided (`longest_genebody.bed`).

# Differential binding analysis
The raw counts generated with Diffbind (or counts generated with `fetch_gb_counts.sh` in the case of Pol II ChIP-seq) were used as input into DESeq2. See Jupyter notebook `../RNA-seq/RNAseq_DESeq2_example_YY1AID_timecourse.ipynb` for the settings and commands used. 



