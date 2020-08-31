################################################################################
################################################################################
## R version 4.0.1 (2020-06-06) -- "See Things Now"
## Copyright (C) 2020 The R Foundation for Statistical Computing
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04 LTS

## BiocManager version: 1.30.10
## Bioconductor version: 3.11
## Date: Wednesday, June 17, 2020.
## Created by Ruwaa I. Mohamed
################################################################################
################################################################################
## Installing Bioconductor and the required libraries
# use packageVersion('') to check the versio of all installed packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::version()
install.packages("BiocManager")

BiocManager::install("DESeq2")              # packageVersion: 1.28.1
BiocManager::install("org.Hs.eg.db")        # packageVersion: 3.11.4
BiocManager::install('EnhancedVolcano')     # packageVersion: 1.6.0
BiocManager::install("biomaRt")             # packageVersion: 2.44.1

# BiocManager::install("genefilter")        # packageVersion: 1.70.0    # already installed with DESeq2
# BiocManager::install("apeglm")              # packageVersion: 1.10.0

install.packages("stringr")       # packageVersion: 1.4.0
install.packages("gplots")        # packageVersion: 3.0.3
install.packages("ggplot2")       # packageVersion: 3.3.2
install.packages("RColorBrewer")  # packageVersion: 1.1.2

# install.packages("cluster")     # packageVersion: 2.1.0
# install.packages("biclust")     # packageVersion: 2.0.2
# install.packages("ggdendro")    # packageVersion: 0.1.20
################################################################################
## Loading the required libraries/scripts
source('scripts/get_df_from_gdc.R')
source('scripts/ENSG_to_symbol.R')
source('scripts/aggregate_rows.R')
source('scripts/proteins.R')
source('scripts/plotting_gene.R')
source('scripts/to_tiff.R')

library(gplots)
library(ggplot2)
library(RColorBrewer)
# library(cluster)
# library(biclust)
# library(ggdendro)

library(DESeq2)
library(EnhancedVolcano)
library(biomaRt)
# library(apeglm)
# library(genefilter)
################################################################################
## Reading the RNA-seq files (1,222 files)
rna.data.path <- "raw_data/RNA-seq/"

## reading the expression values from the raw data and saving them to RDS object (Don't use this is unless you don't have the RDS object)!
rna.exp.df <- get_df_from_gdc(rna.data.path, 'htseq.counts.gz$')
rna.exp.df <- rna.exp.df[1:(nrow(rna.exp.df)-5),]
saveRDS(rna.exp.df, file = "saved_objects/exp.rna.rds")

## reading the expression values from the previously saved RDS object (use this instead)!
rna.exp.df <- readRDS(file = "saved_objects/exp.rna.rds")
View(rna.exp.df) # DON'T!
################################################################################
## Reading the sample sheet
sample.sheet <- read.csv("raw_data/females_only.ductal-and-lobular-only.gdc_sample_sheet.2020-06-17.tsv", sep="\t", header=TRUE)
table(sample.sheet$Sample.Type)

## dropping metastatic samples from the sample sheet
sample.sheet <- sample.sheet[! sample.sheet$Sample.Type == "Metastatic",]
table(sample.sheet$Sample.Type)

## dropping metastatic samples from the RNA-seq expression dataframe
rna.exp.df <- rna.exp.df[, colnames(rna.exp.df) %in% sample.sheet$File.ID]

saveRDS(rna.exp.df, file = "saved_objects/exp.rna.subsetted.rds")
rna.exp.df <- readRDS("saved_objects/exp.rna.subsetted.rds")
################################################################################
## Reading the clinical data 
clinical <- read.csv("raw_data/clinical.project-TCGA-BRCA.2020-03-05/clinical.tsv", header=TRUE, sep="\t", na="--")
clinical <- clinical[, colSums(is.na(clinical)) != nrow(clinical)]
clinical <- clinical[clinical$submitter_id %in% sample.sheet$Case.ID ,]
write.csv(clinical, 'saved_objects/clinical.all.csv', row.names=FALSE)
write.csv(clinical[clinical$submitter_id %in% sample.sheet$Case.ID,], 'saved_objects/clinical_filtered.csv', row.names=FALSE)
################################################################################
## subsetting from clinical and mixing with sample sheet
clinical <- clinical[clinical$primary_diagnosis %in% c("Infiltrating duct carcinoma, NOS", "Lobular carcinoma, NOS"),]
clinical <- unique(clinical[,c("submitter_id", "primary_diagnosis")])
mapperIDs <- match(sample.sheet$Case.ID, clinical$submitter_id)
sample.sheet.clinical <- cbind(sample.sheet, primary_diagnosis=clinical$primary_diagnosis[mapperIDs])
sample.sheet.clinical <- sample.sheet.clinical[! is.na(sample.sheet.clinical$primary_diagnosis),]
rm(mapperIDs)
sample.sheet <- sample.sheet.clinical

write.csv(sample.sheet, "saved_objects/sample.sheet.csv", row.names=FALSE)
sample.sheet <- read.csv("saved_objects/sample.sheet.csv", header=TRUE)
write.csv(clinical, "saved_objects/clinical.csv", row.names=FALSE)
clinical <- read.csv("saved_objects/clinical.csv", header=TRUE)
## 992 individual, 2 columns (with mixed)
## 964 individual (without mixed)

# drop from the RNA-exp DF too
rna.exp.df <- rna.exp.df[, colnames(rna.exp.df) %in% sample.sheet$File.ID]
saveRDS(rna.exp.df, file = "saved_objects/rna.exp.df.final.rds")
rna.exp.df <- readRDS("saved_objects/rna.exp.df.final.rds")
################################################################################
## Changing ENSG to SYMBOL
rna.exp.df.sym <- ENSG_to_symbol(rna.exp.df)
saveRDS(rna.exp.df.sym, file = "saved_objects/rna.exp.df.sym.rds")
rna.exp.df.sym <- readRDS(file = "saved_objects/rna.exp.df.sym.rds")
## 60,483 ENSG where reduced to 25,596 gene symbol (org.Hs.eg.db version 3.11.4)

rna.exp.df.sym.agg <- aggregate_rows(rna.exp.df.sym, agg.var=rna.exp.df.sym$symbol) # This step takes too long to run!
saveRDS(rna.exp.df.sym.agg, file = "saved_objects/rna.exp.df.sym.agg.rds")
rna.exp.df.sym.agg <- readRDS(file="saved_objects/rna.exp.df.sym.agg.rds")
## They're aggregated to 25,531 genes
rm(rna.exp.df.sym)
rm(rna.exp.df)
################################################################################
## Preparing the inputs of DESEQ2
count.data <- rna.exp.df.sym.agg

col.data <- sample.sheet[, colnames(sample.sheet) %in% c("File.ID", "Sample.ID", "Sample.Type", "primary_diagnosis")]

col.data$Sample.Type <- as.factor(col.data$Sample.Type)
col.data$Sample.Type <- relevel(col.data$Sample.Type, "Solid Tissue Normal")
levels(col.data$Sample.Type) <- c("Normal", "Tumor")
table(col.data$Sample.Type)
## 96 Normal and 974 Tumor

col.data$primary_diagnosis <- as.factor(col.data$primary_diagnosis)
levels(col.data$primary_diagnosis) <- c("IDC", "LC")
table(col.data$primary_diagnosis)
## 860 IDC, and 210 LC

# re-ording 
new.order <- match(colnames(count.data), col.data$File.ID)
count.data <- count.data[order(new.order)]
rm(new.order)

# convert to matrix (numbers only)
count.data <- apply (count.data, 2, as.integer)
rownames(count.data) <- rownames(rna.exp.df.sym.agg)

# check the final order match
all(colnames(count.data) == col.data$File.ID)
################################################################################
## Building DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=count.data , colData=col.data , design=~primary_diagnosis+Sample.Type+primary_diagnosis:Sample.Type)
dds

dds$group <- factor(paste(dds$primary_diagnosis, dds$Sample.Type, sep="_"))
design(dds) <- ~ group
dds

table(dds$group)
## IDC_ Normal: 89, IDC_Tumor: 771
## LC_Normal: 7, LC_Tumor:203

saveRDS(dds, "saved_objects/dds.rds")
dds <- readRDS("saved_objects/dds.rds")

## Collapse Technical Replicates
ddsCollapsed <- collapseReplicates(dds, groupby=dds$Sample.ID)
ddsCollapsed
saveRDS(ddsCollapsed, "saved_objects/ddsCollapsed.rds")
ddsCollapsed <- readRDS("saved_objects/ddsCollapsed.rds")

## to check that collapsing worked!
original <- rowSums(counts(dds)[, dds$Sample.ID == "TCGA-A7-A26J-01A"])
all(original == counts(ddsCollapsed)[,"TCGA-A7-A26J-01A"])
rm(original)

## testing
table(ddsCollapsed$group)
## IDC_ Normal: 89, IDC_Tumor: 767 (-4)
## LC_Normal: 7, LC_Tumor:202 (-1)
################################################################################
################################################################################
## Run DESEQ2: on the 1070 samples (on the new 'group' column)
dds.run <- DESeq(ddsCollapsed)
# The following steps were run:
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 3670 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
dds.run
# resultsNames(dds.run)
saveRDS(dds.run, file="saved_objects/dds.run.group.rds")
dds.run <- readRDS("saved_objects/dds.run.group.rds")
################################################################################
################################################################################
## To subset based on BC subtype, run the following section

## Data Normalization for plotting: 1. VST Normalization
dds.vsd <- varianceStabilizingTransformation(dds.run, blind=FALSE)
dds.vsd

saveRDS(dds.vsd, "saved_objects/dds.vsd.rds")
dds.vsd <- readRDS("saved_objects/dds.vsd.rds")

source('scripts/idc.R')
source('scripts/lc.R')

## For P-val 0.05 and LFC 1, out of 25,175 genes
## IDC: 2,261 up-regulated, 1,619 down-regulated
## LC: 162 up-regulated, 119 down-regulated
## Mixed: 387 up-regulated, 620 down-regulated
################################################################################
################################################################################
### Saving all figures from the plots tab at once 
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="saved_objects/figures/")
################################################################################