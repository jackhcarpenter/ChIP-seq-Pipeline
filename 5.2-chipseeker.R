################################################################################
# Loading dependencies
################################################################################

library(BiocManager)

## Install libraries

BiocManager::install(c("ChIPpeakAnno", 
                        "ChIPseeker", 
                        "GenomicFeatures", 
                        "biomaRt", 
                        "biomaRt"))

## Load libraries
library(ChIPpeakAnno)
library(ChIPseeker)
library(GenomicFeatures)
library(biomaRt)
library(gprofiler2)
library(tidyverse)

################################################################################
# Setting up environment
################################################################################

## Setting workingdir

setwd("your/working/dir")

## Making dir for plots

dir.create("plots")

dir.create("beds")

## Create tx database for Araport annotation (Note: example for A. thaliana. It will be a different mart if it's not a plant)

txdb <- makeTxDbFromBiomart(biomart = "plants_mart", 
                            host = "https://plants.ensembl.org", 
                            dataset = "athaliana_eg_gene")
## Load data

samplefiles <- list.files(pattern = ".bed", 
                            full.names = T)
                            
samplefiles <- as.list(samplefiles)

names(samplefiles) <- "TCP4"


################################################################################
# Main CMDs
################################################################################

## Annotation

peakAnnoList <- lapply(samplefiles, 
                        annotatePeak, 
                        TxDb = txdb, 
                        tssRegion = c(-5000,5000), 
                        verbose = FALSE)


## Coverage of genome for each sample

peak_sample <- readPeakFile(samplefiles[[1]])

## Making the plot
coverage_sample <- covplot(peak_sample, weightCol = peak_sample$V5)

## Having a look
coverage_sample

## Saving plots as .pdf and .png
ggsave("coverage_sample.pdf", plot = coverage_sample, path = "plots/", dpi = 800)
ggsave("coverage_sample.png", plot = coverage_sample, path = "plots/", dpi = 800)

## Same for control or second sample in samplefiles
peak_control <- readPeakFile(samplefiles[[2]])
coverage_control <- covplot(peak_control, weightCol = peak_control$V5)
coverage_control
ggsave("coverage_control.pdf", plot = coverage_control, path = "plots", dpi = 800)
ggsave("coverage_control.png", plot = coverage_control, path = "plots", dpi = 800)

## Retrieving all typical promoters from txdb annotation

promoter <- getPromoters(TxDb = txdb, upstream = 5000, downstream = 5000)

## Getting peak information from around the promoter regions

tagMatrix_sample <- getTagMatrix(peak_sample, windows = promoter)

## Plotting

avgPlot_sample <- plotAvgProf(tagMatrix_sample, 
                                xlim = c(-5000, 5000), 
                                conf = 0.95, 
                                resample = 1000, 
                                xlab = "Genomic Region (5' -> 3')", 
                                ylab = "Read Count Frequency")

## Saving plots as .pdf and .png
ggsave("avgPlot_sample.pdf", plot = avgPlot_sample, path = "plots", dpi = 800)
ggsave("avgPlot_sample.png", plot = avgPlot_sample, path = "plots", dpi = 800)

## Same for control or second sample in samplefiles
tagMatrix_control <- getTagMatrix(peak_control, windows = promoter)
avgPlot_control <- plotAvgProf(tagMatrix_control, xlim = c(-5000, 5000), conf = 0.95, resample = 1000, xlab = "Genomic Region (5' -> 3')", ylab = "Read Count Frequency")
ggsave("avgPlot_control.pdf", plot = avgPlot_control, path = "plots", dpi = 800)
ggsave("avgPlot_control.png", plot = avgPlot_control, path = "plots", dpi = 800)


## Visualization of genomic feature representation

comparison_bar <- plotAnnoBar(peakAnnoList)

comparison_to_TSS <- plotDistToTSS(peakAnnoList, title = "Distribution of TF-binding loci relative to TSS")

bar2_1 <- plotDistToTSS(peakAnnoList[[1]], title = "Distribution of TF-binding loci relative to TSS")
ggsave("bar2_1.pdf", plot = bar2_1, path = "plots", dpi = 800)
ggsave("bar2_1.png", plot = bar2_1, path = "plots", dpi = 800)

bar2_2 <- plotDistToTSS(peakAnnoList[[2]], title = "Distribution of TF-binding loci relative to TSS")
ggsave("bar2_2.pdf", plot = bar2_2, path = "plots", dpi = 800)
ggsave("bar2_2.png", plot = bar2_2, path = "plots", dpi = 800)

peak1 <- peakAnnoList[[1]]
plotAnnoPie(peak1)

peak2 <- peakAnnoList[[2]]
plotAnnoPie(peak2)


## Retrieve annotation
sample_annot <- data.frame(peakAnnoList[["sample"]]@anno)
control_annot <- data.frame(peakAnnoList[["control"]]@anno)

## Write down gene names in a table
sample.df <- data.frame(peakAnnoList[["sample"]])
write.table(sample.df, "beds/sample_peaks.txt", sep = "\t")
control.df <- data.frame(peakAnnoList[["control"]])
write.table(control.df, "plots/control_peaks.txt", sep = "\t")

            
## Functional enrichment (GO)
gconv <- gconvert(query = sample.df$geneId, organism = "athaliana", #Note: change to your organism of interest 
                  target = "ENSG", mthreshold = Inf, filter_na = TRUE)
colnames(gconv)[2] <- "geneId"
sample.annot <- merge(sample.df, gconv, by='geneId')
write.table(sample.annot, "beds/annotated_sample_peaks.txt", sep = "\t")


gconv <- gconvert(query = control.df$geneId, organism = "athaliana", target = "ENSG",
                  mthreshold = Inf, filter_na = TRUE)
colnames(gconv)[2] <- "geneId"
control.annot <- merge(control.df, gconv, by='geneId')
write.table(control.annot, "beds/annotated_control_peaks.txt", sep = "\t")
