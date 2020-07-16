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
# BiocManager::install("genefilter")        # packageVersion: 1.70.0    # already installed with DESeq2
BiocManager::install("apeglm")              # packageVersion: 1.10.0

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

library(DESeq2)
library(EnhancedVolcano)
library(apeglm)
# library(genefilter)

library(gplots)
library(ggplot2)
library(RColorBrewer)
# library(cluster)
# library(biclust)
# library(ggdendro)
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
sample_sheet <- read.csv("raw_data/females_only.ductal-and-lobular-only.gdc_sample_sheet.2020-06-17.tsv", sep="\t", header=TRUE)
table(sample_sheet$Sample.Type)

## dropping metastatic samples from the sample sheet
sample_sheet <- sample_sheet[! sample_sheet$Sample.Type == "Metastatic",]
table(sample_sheet$Sample.Type)

## dropping metastatic samples from the RNA-seq expression dataframe
rna.exp.df.subsetted <- rna.exp.df[, colnames(rna.exp.df) %in% sample_sheet$File.ID]

saveRDS(rna.exp.df.subsetted, file = "saved_objects/exp.rna.subsetted.rds")
rna.exp.df.subsetted <- readRDS("saved_objects/exp.rna.subsetted.rds")
################################################################################
## Reading the clinical data 
clinical <- read.csv("raw_data/clinical.project-TCGA-BRCA.2020-03-05/clinical.tsv", header=TRUE, sep="\t", na="--")
clinical <- clinical[, colSums(is.na(clinical)) != nrow(clinical)]
clinical <- clinical[clinical$submitter_id %in% sample_sheet$Case.ID ,]
write.csv(clinical, 'saved_objects/clinical.all.csv', row.names=FALSE)
write.csv(clinical[clinical$submitter_id %in% sample_sheet$Case.ID,], 'saved_objects/clinical_filtered.csv', row.names=FALSE)

exposure <- read.csv("raw_data/clinical.project-TCGA-BRCA.2020-03-05/exposure.tsv", header=TRUE, sep="\t", na="--")
exposure <- exposure[, colSums(is.na(exposure)) != nrow(exposure)]
exposure <- exposure[exposure$submitter_id %in% sample_sheet$Case.ID ,]
# we don't need the alcohol too! (all not reported)

family_history <- read.csv("raw_data/clinical.project-TCGA-BRCA.2020-03-05/family_history.tsv", header=TRUE, sep="\t", na="--")
# Empty file!
################################################################################
## subsetting from clinical and mixing with sample sheet
clinical.types <- clinical[clinical$primary_diagnosis %in% c("Infiltrating duct and lobular carcinoma", "Infiltrating duct carcinoma, NOS", "Lobular carcinoma, NOS"),]
y <- unique(clinical.types[,c("submitter_id", "primary_diagnosis")])
mapperIDs <- match(sample_sheet$Case.ID, y$submitter_id)
sample.sheet.clinical.types <- cbind(sample_sheet, primary_diagnosis=y$primary_diagnosis[mapperIDs])
sample.sheet.clinical.types.uniq <- sample.sheet.clinical.types[! is.na(sample.sheet.clinical.types$primary_diagnosis),]
sample.sheet <- sample.sheet.clinical.types.uniq

# drop from the RNA-exp DF too
rna.exp.df.subsetted <- rna.exp.df.subsetted[, colnames(rna.exp.df.subsetted) %in% sample.sheet.clinical.types.uniq$File.ID]
################################################################################
## Changing ENSG to SYMBOL
rna.exp.df.sym <- ENSG_to_symbol(rna.exp.df.subsetted)
saveRDS(rna.exp.df.sym, file = "saved_objects/rna.exp.df.sym.rds")
rna.exp.df.sym <- readRDS(file = "saved_objects/rna.exp.df.sym.rds")
## 60,483 ENSG where reduced to 25,596 gene symbol (org.Hs.eg.db version 3.11.4)

rna.exp.df.sym.agg <- aggregate_rows(rna.exp.df.sym, agg.var=rna.exp.df.sym$symbol) # This step takes too long to run!
saveRDS(rna.exp.df.sym.agg, file = "saved_objects/rna.exp.df.sym.agg.rds")
rna.exp.df.sym.agg <- readRDS(file="saved_objects/rna.exp.df.sym.agg.rds")
rm(rna.exp.df.sym)
################################################################################
## Preparing the inputs of DESEQ2
count.data <- rna.exp.df.sym.agg

col.data <- sample.sheet[, colnames(sample.sheet) %in% c("File.ID", "Sample.ID", "Sample.Type", "primary_diagnosis")]

col.data$Sample.Type <- as.factor(col.data$Sample.Type)
col.data$Sample.Type <- relevel(col.data$Sample.Type, "Solid Tissue Normal")
levels(col.data$Sample.Type) <- c("Normal", "Tumor")
table(col.data$Sample.Type)

col.data$primary_diagnosis <- as.factor(col.data$primary_diagnosis)
levels(col.data$primary_diagnosis) <- c("Mixed", "IDC", "LC")

# reording 
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
################################################################################
## Run DESEQ2 
s <- Sys.time()
dds.run.type <- DESeq(ddsCollapsed)
e <- Sys.time()
# The following steps were run:
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 3267 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
e-s
dds.run.type
resultsNames(dds.run.type)
saveRDS(dds.run.type, file="saved_objects/dds.run.type.rds")
dds.run.type <- readRDS("saved_objects/dds.run.type.rds")
## testing
table(dds.run.type$group)
################################################################################
################################################################################
## To subset based on BC subtype, run the following scripts

## Data Normalization for plotting: 1. VST Normalization
dds.vsd <- varianceStabilizingTransformation(dds.run, blind=FALSE)
dds.vsd

source('scripts/idc.R')
source('scripts/lc.R')
source('scripts/mixed.R')

## For P-val 0.05 and LFC 1, out of 25,175 genes
## IDC: 2,261 up-regulated, 1,619 down-regulated
## LC: 162 up-regulated, 119 down-regulated
## Mixed: 387 up-regulated, 620 down-regulated

## continue the rest of this scipt only if you want to ignore BC subtyps and deal with tumor vs. normal only
################################################################################
################################################################################
## Creating the results (RES) object
design(ddsCollapsed) <- ~Sample.Type
s <- Sys.time()
dds.run.type <- DESeq(ddsCollapsed)
e <- Sys.time()
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 3906 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
e-s
rm(s)
rm(e)

res <- results(dds.run.type, contrast=c("Sample.Type", "Tumor", "Normal"), alpha=0.05,lfcThreshold=1)
summary(res)
# out of 25,175 genes: 2,220 up-regulated and 1,608 down-regulated

## Removing incomplete genes (No need for this step!)
res.complete <- res[complete.cases(res),]
summary(res.complete)

## saving to and reading from RDS object
saveRDS(res, file="saved_objects/res.rds")
res <- readRDS("saved_objects/res.rds")

## Chkeckin the distribution of the LFC and p-adjusted value
summary(res$log2FoldChange)
hist(res$log2FoldChange, main="Distribution of the Log2 fold change in the results", xlab="Log2 Fold Change")
summary(res$padj)
hist(res$padj, main="Distribution of the adjusted p-value in the results", xlab="Adjusted p-value")
################################################################################
## Converting the results to DF
res.df <- as.data.frame(res)
res.df <- res.df[order(res.df$padj),]
write.csv(res.df, "saved_objects/all.results.sorted.csv")
################################################################################
## Extracting DEGs
res.df.degs <- res.df[res.df$padj<0.05 & abs(res.df$log2FoldChange)>1,]
write.csv(res.df.degs, "saved_objects/degs-p.01-LFC2.csv")

res.degs <- res.complete[res.complete$padj<0.05 & abs(res.complete$log2FoldChange)>1,]
summary(res.degs)

summary(res.degs$padj)
hist(res.degs$padj, main="Distribution of the adjusted p-value in the DEGs", xlab="Adjusted p-value")
################################################################################
## Repair genes subsetting
res.rep <- res[rownames(res) %in% repair.genes$symbol,]
summary(res.rep)
## 285 genes, 34 up-regulated, 0 down-regulated

## Extracting DEGs
res.degs.rep <- res.degs[rownames(res.degs) %in% repair.genes$symbol,]
summary(res.degs.rep)

## Chkeckin the distribution of the LFC and p-adjusted value
summary(res.rep$log2FoldChange)
hist(res.rep$log2FoldChange, main="Distribution of the Log2 fold change in the results (Repair Genes Only)", xlab="Log2 Fold Change")
summary(res.rep$padj)
hist(res.rep$padj, main="Distribution of the adjusted p-value in the results (Repair Genes Only)", xlab="Adjusted p-value")

summary(res.degs.rep$log2FoldChange)
hist(res.degs.rep$log2FoldChange, main="Distribution of the Log2 fold change in the DEGs (Repair Genes Only)", xlab="Log2 Fold Change")
summary(res.degs.rep$padj)
hist(res.degs.rep$padj, main="Distribution of the adjusted p-value in the DEGs (Repair Genes Only)", xlab="Adjusted p-value")
################################################################################
## Normalization Methods Available:
## 1. rlog()  --> takes too long and not applicable for our dataset (not as sensetive as vst and takes too long and fails)
## 2. vst() or varianceStabilizingTransformation() --> the first is just a wrapper for the second. 
## 3. normTransform() --> log2 transformation or any other desired fuction. Will not be used!
## 4. lfcShrink
################################################################################
## Data Normalization for plotting: 1. VST Normalization
dds.type.vsd <- varianceStabilizingTransformation(dds.run.type, blind=FALSE)
dds.type.vsd
# the matrix of transformed values is stored in assay(dds.type.vsd)

# Extracting all DEGs
dds.type.vsd.degs <- dds.type.vsd[rownames(dds.type.vsd) %in% rownames(res.degs),]
dds.type.vsd.degs

# Extracting all repair and repair DEGs
dds.type.vsd.rep <- dds.type.vsd[rownames(dds.type.vsd) %in% repair.genes$symbol,]
dds.type.vsd.rep.degs <- dds.type.vsd.rep[row.names(dds.type.vsd.rep) %in% rownames(res.degs.rep)]
################################################################################
## Data Normalization for plotting: 2. lfcShrink Normalization
resultsNames(dds.run.type)
res.lfc <- lfcShrink(dds.run.type, coef=2, res=res, type = "apeglm")
summary(res.lfc)

## Chkeckin the distribution of the p-adjusted value
summary(res.lfc$padj)
hist(res.lfc$padj, main="Distribution of the adjusted p-value in the lfc-shrinked (apeglm) results", xlab="Adjusted p-value")
################################################################################
## Data Normalization for plotting: 3. normTransform
# dds.nrd <- normTransform(dds.run.type)
# dds.nrd
# saveRDS(dds.nrd, "saved_objects/dds.nrd.rds")
# dds.nrd <- readRDS("saved_objects/dds.nrd.rds")
################################################################################
## Plotting: 1. MAplot
## not-normalized data
plotMA(res, alpha=0.05, main="MA plot of the not-mormalized DESeq Results")

## normalized data
plotMA(res.lfc, alpha=0.05, main="MA plot of the lfcschrink-normalized DESeq Results")

# ## repair genes (normalized and not normalized)
# plotMA(res.rep, alpha=0.05, main="MA plot of the not-normalized DESeq Results for repair genes only")
# plotMA(res.lfc.rep, alpha=0.05, main="MA plot of the lfcschrink-normalized DESeq Results for repair genes only")
# 
# ## DEGs
# plotMA(res.degs, alpha=0.05, main="MA plot of the not-normalized all DEGs")
# plotMA(res.lfc.degs, alpha=0.05, main="MA plot of the lfc-shrunk all DEGs")
# 
# plotMA(res.degs.rep, alpha=0.05, main="MA plot of the not-normalized repair DEGs")
# plotMA(res.lfc.rep.degs, alpha=0.05, main="MA plot of the lfc-shrunk repair DEGs")
################################################################################
## Plotting: 2. EnhancedVolcano Plot
# not normalized data
EnhancedVolcano(res, 
                lab = rownames(res), 
                x = 'log2FoldChange', 
                y = 'pvalue',
                title = "Volcano plot of the not-normalized DESeq results",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                legendPosition = "right")

# # normalized data
# EnhancedVolcano(res.lfc, 
#                 lab = rownames(res.lfc), 
#                 x = 'log2FoldChange', 
#                 y = 'pvalue',
#                 title = "Volcano plot of the lfcschrink-normalized DESeq Results",
#                 pCutoff = 0.05,
#                 FCcutoff = 1.0,
#                 legendPosition = "right")
################################################################################
## Plotting: 3. PCA
# should we use normalized data instead? The manual says yes!
# what is the best ntop value? size of the DEGs?

## All genes
pca.plot.vsd <- plotPCA(dds.type.vsd, intgroup="Sample.Type", ntop = 25531, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd, "percentVar"))
ggplot(pca.plot.vsd, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vst-normalized DESeq DataSet (ntop = 25531)")
rm(percentVar)
rm(pca.plot.vsd)

# all repair
pca.plot.vsd.rep <- plotPCA(dds.type.vsd.rep, intgroup="Sample.Type", ntop = 285, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.rep, "percentVar"))
ggplot(pca.plot.vsd.rep, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vsd-normalized DESeq DataSet of the repair genes only (ntop = 285)")
rm(percentVar)
rm(pca.plot.vsd.rep)

## all DEGs
pca.plot.vsd.degs <- plotPCA(dds.type.vsd.degs, intgroup="Sample.Type", ntop = 3828, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.degs, "percentVar"))
ggplot(pca.plot.vsd.degs, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vst-normalized DEGs (ntop = 3828)")
rm(percentVar)
rm(pca.plot.vsd.degs)

## PCA for Repair DEGs
pca.plot.vsdrep.degs <- plotPCA(dds.type.vsd.rep.degs, intgroup="Sample.Type", ntop = 34, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsdrep.degs, "percentVar"))
ggplot(pca.plot.vsdrep.degs, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vst-normalized repair DEGs (ntop = 34)")
rm(percentVar)
rm(pca.plot.vsdrep.degs)
################################################################################
## Plotting: 4. Dispersion estimate
# estimateDispersions(dds)
plotDispEsts(dds.run.type)
################################################################################
## Plotting: 5. meanSdPlot
# meanSdPlot()
################################################################################
## Plotting: 6. Counts plot
# genes.list <- c("PCLAF", "ISG15", "EXO1")
# genes.list <- rownames(res.degs.rep)
# for (gen in genes.list){
#   plotCounts(dds.run.type, gene=gen, intgroup = "Sample.Type") # returnData = FALSE
# }
# rm(genes.list)
# plotCounts(dds.run.type, gene=gen, intgroup = "Sample.Type")
################################################################################
## Plotting: 7. Boxplot with dots
list(assays(dds.run.type)) # counts mu H cooks replaceCounts replaceCooks
genes.list <- rownames(res.degs.rep)

par(mfrow=c(3,3))
for (gen in genes.list){
  plotting_gene(dds.run.type, gen)
}
par(mfrow=c(1,1))
rm(genes.list)
# list(assays(dds.type.vsd)) # counts mu H cooks replaceCounts replaceCooks
# boxplot(t(assays(dds.type.vsd["EXO1"])[[1]])~dds.type.vsd$Sample.Type, 
#         range=0, las=1, log='y', 
#         xlab="Groups", ylab="Counts", main="EXO1 (vsd-normalized)",
#         col=c("darksalmon", "darkred"))
# stripchart(t(assays(dds.type.vsd["EXO1"])[[1]])~dds.type.vsd$Sample.Type, 
#            vertical=TRUE, method='jitter', add=TRUE, pch=23, col="black", cex=0.5)
################################################################################
## Plotting: 8. Heatmap & Dendogram
colors2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)

# for repair degs
heatmap.2(assay(dds.type.vsd.rep.degs), col=colors2,
          scale="row", trace="none", labCol=substring(colnames(dds.type.vsd.rep.degs), 6),
          # dendrogram = "row",
          Colv=order(dds.type.vsd.rep.degs$Sample.Type), Rowv=TRUE,
          ColSideColors = c(Primary.Tumor="darkgreen", Solid.Tissue.Normal="orange")[colData(dds.type.vsd.rep.degs)$Sample.Type],
          key =TRUE, key.title="Heatmap Key",
          main="Heat map of the top 34 differentially expressed repair genes")

# # for all degs
# heatmap.2(assay(dds.type.vsd.degs), col=colors2,
#           scale="row", trace="none", labCol=substring(colnames(dds.type.vsd.degs), 6),
#           # dendrogram = "row",
#           Colv=order(dds.type.vsd.degs$Sample.Type), Rowv=TRUE,
#           ColSideColors = c(Primary.Tumor="darkgreen", Solid.Tissue.Normal="orange")[colData(dds.type.vsd.degs)$Sample.Type],
#           key =TRUE, key.title="Heatmap Key",
#           main="Heat map of the 1,054 differentially expressed genes")
################################################################################
################################################################################
### Saving all figures from the plots tab at once 
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="saved_objects/figures/")
################################################################################