#!/usr/bin/env Rscript
# Run rGMAP on single chromosme 
# arg = name of output from juicer dump 
# must be for a single chromosome 10kbp res kR norm

args = commandArgs(trailingOnly=TRUE)

library(rGMAP);

# read output text file from juicer dump (10kbp KR normalized)
hic_dumpfile <- args[1]
outname <- args[2]

# run rGMAP
res <- rGMAP(hic_dumpfile, dom_order=3)

# write to csv
outname <- paste(args[2], "sweep1_hierTADS.csv", sep="_")
write.csv(res$hierTads, outname, row.names=FALSE) 



# run second sweep of rGMAP to capture smaller structures
res <- rGMAP(hic_dumpfile, dom_order=3, maxDistInBin=50)
outname <- paste(args[2], "sweep2_hierTADS.csv", sep="_")
write.csv(res$hierTads, outname, row.names=FALSE) 
