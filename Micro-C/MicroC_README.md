# Walkthrough of Micro-C analysis 
# Set up python environment
The `env/py38.yml` file can be used with conda to install the software required to run the loop quantification script on a Linux platform.

For other platforms, you can manually install the packages listed in the file (and in the Methods section). Importantly, the loop quantification script uses cooltools version 0.5.1. 

# Input files
The .mcool files and loop lists deposited on GEO can also be used as input for the scripts below.
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247254

# Generate contact maps with distiller
1) Download distiller pipeline files

`nextflow clone mirnylab/distiller-nf .`

2) Add sample names and corresponding datafiles to project.yml. An example with the settings we used is provided.
3) Make appropriate changes to the config files. Our config files are located in the distiller_files directory.
4) Run the pipeline

`sh run_distiller.sh`

# Generate contact versus distance curves
Also known as P(s) curve. This script plots the log-binned decay of interactions with genomic distance.

`python logbin_contact_v_distance.py -i MCOOL_DIR -o OUTFILE.PDF`

MCOOL_DIR is the path to a directory that contains multiple .mcool files to be analyzed. 

# Generating loop calls
Generate expected values and call loops at 5k and 10k resolution. Instead of classic local neighborhoods, we used a 'rounded donut' and round lower left (allows for more loop calls near the diagonal). 

`call_dot_rounddonut.sh -m MCOOL`

# Quantify loop strength across multiple heatmaps
First, create a csv file detailing the heatmaps that should be quantified. It should have a header row with three columns: name, map, expected.

* name = contact map name for the output table columns; 
* map = path to corresponding mcool file; 
* expected = path to corresponding expected file

The example `config_loops10k.csv` is provided.

`python loops_fetch_ooe.py -c ${config} -l ${looplist} -r ${resolution} -p ${cores}`

* config = the config file as described above
* looplist = a csv file that includes loop coordinates and a header row containing bedpe columns ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
* resolution = resolution of contact maps
* cores = number of processors

e.g.,

`python loops_fetch_ooe.py -c config.csv -l loop_list.csv -r 10000 -p 16`

### Output

Loop strengths will be saved in a new csv file named "[config]_ooe_table.csv". This will have the original loop list columns plus the following additional columns containing loop quantifications:

* peak_obs_[name]: observed contacts for the loop
* peak_obs_[name]: observed/expected contacts for the loop
* peak_ooe_donut_[name]: "loop strength" (observed/expected_local) for the loop

The metric we call "loop strength" (observed / local-adjusted expected) was originally defined in Rao et al., 2014:

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5635824/figure/F3/

For our quantification, we calculated the local-adjusted expected based on a "rounded donut" as defined and depicted in the cooltools documentation:

https://cooltools.readthedocs.io/en/latest/notebooks/dots.html#Calling-dots-with-a-%22rounded-donut%22-kernel

# Call domains with rGMAP
`sh call_domains.sh -m CONTACTMAP.mcool`

The first step is creating a new file that is compatible with rGMAP. The second step runs rGMAP in two iterations, one optimized for TADs and one for smaller TADs / subTADs. This script runs on 10kb resolution heatmaps. The resolution is hard-coded in the script, so change it manually if you want to try a different resolution. Note: this uses 8 cores by default. 

# Create merged/filtered domain and boundary lists
`sh combine_domain_lists.sh -m CONTACTMAP1.mcool -e EXPECTED1.tsv -t CONTACTMAP1_domains -n CONTACTMAP2.mcool -f EXPECTED2.tsv -u CONTACTMAP2_domains`

The previous step should have made a directory containing the TAD calls (should be labeled "CONTACTMAP_domains") 
Given 2 conditions (e.g. untreated and auxin), use this script to make a list of boundaries, domains, and insulation scores for each condition. It will also create a merged list of boundaries from both conditions.



