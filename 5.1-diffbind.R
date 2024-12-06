################################################################################
# Loading dependencies
################################################################################

library(BiocManager)

## Install libraries
BiocManager::install(c("GenomeInfoDbData", "DiffBind"))
install.packages(c("tidyverse", "BiocParallel", "ggplot2"))

## Load libraries
library(DiffBind)
library(tidyverse)
library(BiocParallel)
library(ggplot2)

################################################################################
# Setting up environment
################################################################################

# Setting workingdir

setwd("your/working/dir")

## Reading in peaksets and associated metadata
# Note: control first and sample at the end

samples <- read.csv("metadata/samplesheetTCP4s.csv")
samples
dbObj <- dba(sampleSheet=samples)  #build a consensus set
dbObj #have a look at consensus set


################################################################################
# Main CMDs
################################################################################

## Compute count information for each of the peaks in the consensus set
#Note: increase memory to 500 000 with memory.limit(size = 500000) before executing

dbObj <- dba.count(dbObj, 
                    bParallel = FALSE, 
                    bUseSummarizeOverlaps = TRUE)
dbObj

dba.plotPCA(dbObj, 
            attributes=DBA_CONDITION, 
            label = DBA_ID)

plot(dbObj)

## Establishing a contrast (which samples we want to compare to one another)
## In this case, the contrast is the IP samples against the NOAB (-ve) control

dbObj.contrast <- dba.contrast(dbObj, 
                                categories = DBA_CONDITION, 
                                minMembers = 2)

## In this case, the analysis will only be run with DESEQ2 but method can be set to
## method = DBA_ALL_METHODS to use EDGER as well

dbObj.contrast <- dba.analyze(dbObj.contrast, 
                                method = DBA_DESEQ2)

## Examining which method is better and what peak distribution is between samples
## ONLY RUN THIS WHEN METHODS SET TO DBA_ALL_METHODS

dba.show(dbObj.contrast, 
            bContrasts = T)

dba.plotPCA(dbObj.contrast, 
            contrast = 1, 
            method = DBA_DESEQ2, 
            attributes = DBA_CONDITION, 
            label = DBA_ID)

dba.plotPCA(dbObj.contrast, 
            contrast = 1, 
            method = DBA_EDGER, 
            attributes = DBA_CONDITION, 
            label = DBA_ID)

dba.plotVenn(dbObj.contrast, 
                contrast=1, 
                method = DBA_ALL_METHODS)

dba.plotMA(dbObj.contrast, 
            method = DBA_DESEQ2)

dba.plotMA(dbObj.contrast, 
            bXY = TRUE)

pvals <- dba.plotBox(dbObj.contrast)


## Extracting results

res_deseq <- dba.report(dbObj.contrast, 
                        method = DBA_DESEQ2, 
                        contrast = 1, th = 1)
# note: result is a GRanges object

## Write out file

out <- as.data.frame(res_deseq) #convert GRanges object to dataframe
write.table(out, 
            file = "results/samplesTCP4s_DESeq2.txt", 
            sep = "\t", 
            quote = F, 
            row.names = F)


## Create bed files for ip-enriched keeping only significant peaks (p < 0.05)

ip_enrich <- out %>%
  filter(FDR < 0.05 & Fold < 0) %>%
  select(seqnames, start, end)

## Write out bed file

write.table(ip_enrich, 
            file = "results/TCP4enriched.bed", 
            sep = "\t", 
            quote = F, 
            row.names = F, 
            col.names = F)


#################################################################################
# End
#################################################################################