# Code for Lam et al., 2024)

# Set up python environment
The `spec-file.txt` file can be used with conda to install the software required to run the loop quantification script on a Linux platform.

`conda create --name <env> --file <this file>`

For other platforms, you can manually install the packages listed in the file (and in the Methods section). Importantly, the loop quantification script uses cooltools version 0.5.1. 

# Create a config file
First, create a csv file detailing the heatmaps that should be quantified. It should have a header row with three columns: name, map, expected.

* name = contact map name for the output table columns; 
* map = path to corresponding mcool file; 
* expected = path to corresponding expected file

The example `config_loops10k.csv` is provided. 

# Quantify loop strength across multiple heatmaps

`python loops_fetch_ooe.py -c ${config} -l ${looplist} -r ${resolution} -p ${cores}`

* config = the config file as described above
* looplist = a csv file that includes loop coordinates and a header row containing bedpe columns ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
* resolution = resolution of contact maps
* cores = number of processors

e.g.,

`python loops_fetch_ooe.py -c config.csv -l loop_list.csv -r 10000 -p 16`

# Input files


.mcool files deposited on GEO can be used as input for this script.
