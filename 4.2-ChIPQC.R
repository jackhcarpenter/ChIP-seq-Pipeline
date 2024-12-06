################################################################################
# Loading dependencies
################################################################################

library(BiocManager)
BiocManager::install(c("GenomicFeatures", "ChIPQC", "txdbmaker" , "biomaRt"))

## Load libraries

library(BiocParallel)
library(GenomicFeatures)
library(ChIPQC)
library(GenomeInfoDb)
library(GenomicRanges)

################################################################################
# Setting up environment
################################################################################

# Setting workingdir

setwd("your/working/dir")

# Loading in sample data (a metadata file needs to be generated in the 
# appropriate format)

samples <- read.csv("metadata/samplesheet.csv")

# Check the format is as expected

View(samples)

################################################################################
# Main CMDs
################################################################################

## Create ChIPQC object

# make sure the ChIPQC runs in serial instead of parallel otherwise it crashes
register(SerialParam())

chipObj <- ChIPQC(samples, BPPARAM = SerialParam())

# Generate a ChIPQC Report

ChIPQCreport(chipObj, reportName="ChIP QC report: Col0 and TCP4", reportFolder="ChIPQCreport")

#################################################################################
# End
#################################################################################