################################################################################
################################ Session Info. #################################
## R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
## Copyright (C) 2020 The R Foundation for Statistical Computing
## Platform: x86_64-pc-linux-gnu (64-bit)

## Created by Ruwaa I. Mohamed
## Date: Saturday, May 30, 2020.
################################################################################
## Installing Bioconductor and the required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("DESeq2")              # packageVersion: 1.28.1
################################################################################
## Loading the required libraries/scripts
source('scripts/get_df_from_gdc_files.R')

library(DESeq2)
################################################################################
## Reading the sample sheet
#mirna.sample_sheet <- read.csv("saved_objects/mirna.sample_sheet.csv", header=TRUE)
mirna.sample_sheet <- read.csv("raw_data/gdc_sample_sheet.2020-07-08.tsv", sep="\t", header=TRUE)
mirna.sample_sheet <- mirna.sample_sheet[! mirna.sample_sheet$Sample.Type == "Metastatic",]
mirna.sample_sheet$Sample.Type <- as.factor(mirna.sample_sheet$Sample.Type)
mirna.sample_sheet$Sample.Type <- relevel(mirna.sample_sheet$Sample.Type, "Solid Tissue Normal")
levels(mirna.sample_sheet$Sample.Type) <- c("Normal", "Tumor")
table(mirna.sample_sheet$Sample.Type)
################################################################################
## Reading the files
mirna.exp.df <- get_df_from_gdc_mirna(mirna.path="raw_data/miRNA-seq", file_ext=".mirnas.quantification.txt$")
mirna.exp.df.sub <- mirna.exp.df[, colnames(mirna.exp.df) %in% mirna.sample_sheet$File.ID]

# comparing with the RNA-seq data
summary(mirna.sample_sheet$Case.ID %in% sample_sheet$Case.ID)
file.ids <- mirna.sample_sheet[mirna.sample_sheet$Case.ID %in% sample_sheet$Case.ID,]$File.ID
mirna.exp.df.sub <- mirna.exp.df.sub[, colnames(mirna.exp.df.sub) %in% file.ids]
mirna.sample_sheet <- mirna.sample_sheet[mirna.sample_sheet$File.ID %in% file.ids,]
################################################################################
## Preparing DESeq2 Data
count.data.mirna <- mirna.exp.df.sub

col.data.mirna <- mirna.sample_sheet[, colnames(mirna.sample_sheet) %in% c("File.ID", "Sample.Type", "Sample.ID")]

# reording 
new.order <- match(colnames(count.data.mirna), col.data.mirna$File.ID)
count.data.mirna <- count.data.mirna[order(new.order)]
rm(new.order)

# convert to matrix (numbers only)
count.data.mirna <- apply (count.data.mirna, 2, as.integer)
rownames(count.data.mirna) <- rownames(mirna.exp.df.sub)

# check the final order match
all(colnames(count.data.mirna) == col.data.mirna$File.ID)
################################################################################
## Building DESeqDataSet
dds.mirna <- DESeqDataSetFromMatrix(countData=count.data.mirna , colData=col.data.mirna , design=~Sample.Type)
dds.mirna
saveRDS(dds.mirna, "saved_objects/dds.mirna.rds")
dds.mirna <- readRDS("saved_objects/dds.mirna.rds")

## Collapse Technical Replicates
ddsCollapsed.mirna <- collapseReplicates(dds.mirna, groupby=dds.mirna$Sample.ID)
ddsCollapsed.mirna
saveRDS(ddsCollapsed.mirna, "saved_objects/ddsCollapsed.mirna.rds")
ddsCollapsed.mirna <- readRDS("saved_objects/ddsCollapsed.mirna.rds")

## to check that collapsing worked!
original <- rowSums(counts(dds.mirna)[, dds.mirna$Sample.ID == "TCGA-A7-A0DB-01A"])
all(original == counts(ddsCollapsed.mirna)[,"TCGA-A7-A0DB-01A"])
rm(original)
################################################################################
## Run DESEQ2
dds.mirna.run <- DESeq(ddsCollapsed.mirna)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 109 genes
# -- DESeq argument 'minReplicatesForReplace' = 7
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
dds.mirna.run

saveRDS(dds.mirna.run, file="saved_objects/dds.mirna.run.rds")
dds.mirna.run <- readRDS("saved_objects/dds.mirna.run.rds")
################################################################################
## Creating the results (RES) object
res.mirna <- results(dds.mirna.run, contrast=c("Sample.Type", "Tumor", "Normal"), alpha=0.05,lfcThreshold=1)
summary(res.mirna)
# out of 1600 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 1.00 (up)    : 87, 5.4%
# LFC < -1.00 (down) : 64, 4%
# outliers [1]       : 0, 0%
# low counts [2]     : 676, 42%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

## Removing incomplete genes (No need for this step!)
# res.mirna.2 <- res.mirna[complete.cases(res.mirna),]
# summary(res.mirna.2)

## saving to and reading from RDS object
saveRDS(res.mirna, file="saved_objects/res.mirna.rds")
res.mirna <- readRDS("saved_objects/res.mirna.rds")

## Chkeckin the distribution of the LFC and p-adjusted value
summary(res.mirna$log2FoldChange)
hist(res.mirna$log2FoldChange, main="Distribution of the Log2 fold change in the miRNA results", xlab="Log2 Fold Change")
summary(res.mirna$padj)
hist(res.mirna$padj, main="Distribution of the adjusted p-value in the miRNA results", xlab="Adjusted p-value")
