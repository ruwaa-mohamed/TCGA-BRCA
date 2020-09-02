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
BiocManager::install('EnhancedVolcano')     # packageVersion: 1.6.0

install.packages("RColorBrewer")  # packageVersion: 1.1.2
install.packages("gplots")        # packageVersion: 3.0.3
install.packages("openxlsx")      # packageVersion: 4.1.5
install.packages("VennDiagram")   # packageVersion: 1.6.20
################################################################################
## Loading the required libraries/scripts
source('scripts/get_df_from_gdc_files.R')
source("scripts/proteins.R")

# library(RColorBrewer)
# library(gplots)
# library(ggplot2)
library(openxlsx)

library(DESeq2)
# library(EnhancedVolcano)
################################################################################
## Reading the sample sheet 
## (after adding clinical data and sync with RNA-seq sample.sheet)
mirna.sample.sheet <- read.csv('saved_objects/mirna.sample.sheet.csv', header=TRUE)

## As factor and Releveling
mirna.sample.sheet$Sample.Type <- as.factor(mirna.sample.sheet$Sample.Type)
levels(mirna.sample.sheet$Sample.Type)
table(mirna.sample.sheet$Sample.Type)
# 88 Normal and 962 Tumor

mirna.sample.sheet$primary_diagnosis <- as.factor(mirna.sample.sheet$primary_diagnosis)
levels(mirna.sample.sheet$primary_diagnosis)
table(mirna.sample.sheet$primary_diagnosis)
# 842 IDC and 208 LC
################################################################################
## Reading the miRNA-seq files (1,207 files)
mirna.exp.df <- get_df_from_gdc_mirna(mirna.path="raw_data/miRNA-seq", file_ext=".mirnas.quantification.txt$")
mirna.exp.df <- mirna.exp.df[, colnames(mirna.exp.df) %in% mirna.sample.sheet$File.ID]

saveRDS(mirna.exp.df, "saved_objects/mirna.exp.df.rds")
mirna.exp.df <- readRDS("saved_objects/mirna.exp.df.rds")
################################################################################
## Preparing DESeq2 Data
mirna.count.data <- mirna.exp.df

mirna.col.data <- mirna.sample.sheet[, colnames(mirna.sample.sheet) %in% c("File.ID", "Sample.ID", "Sample.Type", "primary_diagnosis")]

# reording 
new.order <- match(colnames(mirna.count.data), mirna.col.data$File.ID)
mirna.count.data <- mirna.count.data[order(new.order)]
rm(new.order)

# convert to matrix (numbers only)
mirna.count.data <- apply (mirna.count.data, 2, as.integer)
rownames(mirna.count.data) <- rownames(mirna.exp.df)

# check the final order match
all(colnames(mirna.count.data) == mirna.col.data$File.ID)

# removing original DF
rm(mirna.exp.df)
rm(mirna.sample.sheet)
################################################################################
## Building DESeqDataSet
mirna.dds <- DESeqDataSetFromMatrix(countData=mirna.count.data, colData=mirna.col.data , design=~Sample.Type)
mirna.dds
saveRDS(mirna.dds, "saved_objects/mirna.dds.rds")
mirna.dds <- readRDS("saved_objects/mirna.dds.rds")
rm(mirna.count.data, mirna.col.data)

mirna.dds$group <- factor(paste(mirna.dds$primary_diagnosis, mirna.dds$Sample.Type, sep="_"))
table(mirna.dds$group)
# IDC_Normal 82, IDC_Tumor 760
# LC_Normal 6, LC_Tumor 202
design(mirna.dds) <- ~ group
mirna.dds
# 1881 miRs and 1050 samples

## Collapse Technical Replicates
mirna.ddsCollapsed <- collapseReplicates(mirna.dds, groupby=mirna.dds$Sample.ID)
mirna.ddsCollapsed
# 1881 miRs and 1045 samples
# IDC_Normal 82, IDC_Tumor 756
# LC_Normal 6, LC_Tumor 201
table(mirna.ddsCollapsed$group)
saveRDS(mirna.ddsCollapsed, "saved_objects/mirna.ddsCollapsed.rds")
mirna.ddsCollapsed <- readRDS("saved_objects/mirna.ddsCollapsed.rds")

## to check that collapsing worked!
original <- rowSums(counts(mirna.dds)[, mirna.dds$Sample.ID == "TCGA-A7-A0DB-01A"])
all(original == counts(mirna.ddsCollapsed)[,"TCGA-A7-A0DB-01A"])
rm(original)

rm(mirna.dds)
################################################################################
################################################################################
## Run DESEQ2 (Subtypes)
mirna.dds.run <- DESeq(mirna.ddsCollapsed)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 101 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
mirna.dds.run

rm(mirna.ddsCollapsed)

saveRDS(mirna.dds.run, file="saved_objects/mirna.dds.run.rds")
mirna.dds.run <- readRDS("saved_objects/mirna.dds.run.rds")
################################################################################
## varianceStabilizingTransformation
mirna.dds.vsd <- varianceStabilizingTransformation(mirna.dds.run, blind=FALSE)
mirna.dds.vsd

saveRDS(mirna.dds.vsd, file="saved_objects/mirna.dds.vsd.rds")
mirna.dds.vsd <- readRDS("saved_objects/mirna.dds.vsd.rds")
################################################################################
################################################################################
## mirTarBase (experimentally validated only)
## Reading the database
hsa_MTI <- read.xlsx("raw_data/miRNA-DBs/hsa_MTI.xlsx")
substring(hsa_MTI$miRNA, 7) = 'r'

mir.gene.pairs <- unique(hsa_MTI[, c("miRNA", "Target.Gene")])
write.table(mir.gene.pairs, "saved_objects/mir.gene.pairs.csv", sep=",", quote=FALSE)

mir.gene.pairs.rep <- mir.gene.pairs[mir.gene.pairs$Target.Gene %in% repair.genes$symbol,]
write.table(mir.gene.pairs.rep, "saved_objects/mir.gene.pairs.rep.csv", sep=",", quote=FALSE)

rm(hsa_MTI, mir.gene.pairs, mir.gene.pairs.rep)
################################################################################