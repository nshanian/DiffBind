#!/usr/bin/env Rscript

# The `diffbind_atc.R` script is written to run in an `R` environment as a part of SLURM submission on a high performance computing (HPC) cluster. 
# It is NOT written to run locally using `R` or `RStudio` since the necessary `.bam` files called in the `sampleSheet` are typically too large. 
# However, you may try running it the 'Terminal' portion of `RStudio` once you connect to a cluster that can run `R`. 

# `diffbind_atac.R` `R` sciprt will run as part of `diffbind.sh` shell script, available in [ATAC-seq](https://github.com/nshanian/ATAC-seq) repository.

# Download and modify the following paths in `diffbind.sh`:
  
# dir=/<path_to_working_directory>/diffbind
# sampleSheet=<path_to_working_directory>/diffbind/prop_atac_diffbind.csv
# name=atac

# Save this `R` script in your working directory, together with the modified `diffbind.sh` shell script and `sampleSheet` metafile listing covariates and experimental factors. 
# Examples of both are provided in this repository. 

# Run the `diffind.sh` shell script from terminal using `sbatch` `diffbind.sh` and make sure the directory also contains either the raw data in `.bam` and `.narrowPeak` file formats (or paths to directories of raw data in the metafile), the `sampleSheet` metafile in `.csv` format, and the `diffbind_chip.R` script. 
# Make sure that `R` v3.6 or higher is available on the HPC cluster. Please see instructions in `diffind.sh` provided in this repository for more details.

# The script will first perform differential binding analysis based on a `DESeq2` method (`edgeR` is also available as an option). 
# A minimum of two replicates per condition or factor are necessary for differential analysis. 
# Once differential binding analysis is complete, a correlation matrix and scatter XY plot in `.pdf` format will be generated, along with a final report in `.csv` format. 

#  As a final output, an `.RData` data object will also be generated that will contain all the necessary information and can be loaded into RStudio without the `sampleSheet` or final report.

# `DiffBind` package is required. 

# If install.packages("DiffBind") command fails in newer versions of R, uncomment and run the commands below:
# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("DiffBind")

args = commandArgs(trailingOnly=TRUE)

# set working directory
setwd(args[1])

library(DiffBind)

# read in sampleSheet
diff <- dba(sampleSheet=args[2])
diff

# count reads in peaks
diff <- dba.count(diff, summits=250)
diff

# set up contrast
diff <- dba.contrast(diff, categories=c(DBA_FACTOR,DBA_CONDITION), minMembers=2)
diff

# perform Deseq2-based differential analysis
diff <- dba.analyze(diff)
diff

# save your work
save.image(paste(args[3], ".RData", sep=""))

# make cross correlation matrix
pdf(paste("heatmap3_", args[3], ".pdf", sep=""))
plot(diff, contrast=1)
dev.off()

# make scatter XY plot, flip axis
pdf(paste("plotXY_", args[3], ".pdf", sep=""))
dba.plotMA(diff, bXY=TRUE, bFlip=TRUE)
dev.off()

# write csv file of differentially bound regions
report <- dba.report(diff)
write.csv(report, paste(args[3], "_report.csv", sep=""))

# save your work again
save.image(paste(args[3], ".RData", sep=""))


