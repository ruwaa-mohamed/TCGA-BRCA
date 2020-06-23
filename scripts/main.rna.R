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

col.data <- sample_sheet[, colnames(sample_sheet) %in% c("File.ID", "Sample.Type", "Sample.ID")]
col.data$Sample.Type <- as.factor(col.data$Sample.Type)
col.data$Sample.Type <- relevel(col.data$Sample.Type, "Solid Tissue Normal")
levels(col.data$Sample.Type) <- c("Solid.Tissue.Normal", "Primary.Tumor")
table(col.data$Sample.Type)

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
dds <- DESeqDataSetFromMatrix(countData=count.data , colData=col.data , design=~Sample.Type)
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
################################################################################
## Run DESEQ2 
dds.run <- DESeq(ddsCollapsed)     ## This step takes almost 4 hours to run
dds.run
# The following steps were run:
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 3919 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
saveRDS(dds.run, file="saved_objects/dds.run.rds")
dds.run <- readRDS("saved_objects/dds.run.rds")
################################################################################
## Creating the results (RES) object
res <- results(dds.run, contrast=c("Sample.Type", "Primary.Tumor", "Solid.Tissue.Normal"), alpha=0.05,lfcThreshold=1)
summary(res)
## 25,177 genes, 2,233 up-regulated, 1,639 down-regulated

## Removing incomplete genes (No need for this step!)
# res <- res[complete.cases(res),]
# summary(res)

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
res.df.degs <- res.df[res.df$padj<0.05 & abs(res.df$log2FoldChange)>log2(2),]
write.csv(res.degs, "saved_objects/degs-p.05-LFC1.csv")

res.degs <- res[rownames(res) %in% rownames(res.df.degs),]
summary(res.degs)

summary(res.degs$log2FoldChange)
hist(res.degs$log2FoldChange, main="Distribution of the Log2 fold change in the DEGs", xlab="Log2 Fold Change")
summary(res.degs$padj)
hist(res.degs$padj, main="Distribution of the adjusted p-value in the DEGs", xlab="Adjusted p-value")
################################################################################
## Repair genes subsetting (only one of the two following sets!)

# dds.run.rep <- dds.run[rownames(dds.run) %in% repair.genes$symbol,]
# res.rep <- results(dds.run.rep, contrast=c("Sample.Type", "Primary.Tumor", "Solid.Tissue.Normal"), alpha=0.05,lfcThreshold=1)
# summary(res.rep)
# res.rep <- res.rep[complete.cases(res.rep),]
# summary(res.rep)

res.rep <- res[rownames(res) %in% repair.genes$symbol,]
summary(res.rep)
## 285 genes, 35 up-regulated, 0 down-regulated
res.rep.degs <- res.rep[res.rep$padj<0.05 & abs(res.rep$log2FoldChange)>log2(2),]
summary(res.rep.degs)

## Chkeckin the distribution of the LFC and p-adjusted value
summary(res.rep$log2FoldChange)
hist(res.rep$log2FoldChange, main="Distribution of the Log2 fold change in the results (Repair Genes Only)", xlab="Log2 Fold Change")
summary(res.rep$padj)
hist(res.rep$padj, main="Distribution of the adjusted p-value in the results (Repair Genes Only)", xlab="Adjusted p-value")

## Extracting DEGs
res.rep.degs <- res.degs[rownames(res.degs) %in% repair.genes$symbol,]
summary(res.rep.degs)

summary(res.rep.degs$log2FoldChange)
hist(res.rep.degs$log2FoldChange, main="Distribution of the Log2 fold change in the DEGs (Repair Genes Only)", xlab="Log2 Fold Change")
summary(res.rep.degs$padj)
hist(res.rep.degs$padj, main="Distribution of the adjusted p-value in the DEGs (Repair Genes Only)", xlab="Adjusted p-value")
################################################################################
## Normalization Methods Available:
## 1. rlog()  --> takes too long and not applicable for our dataset (not as sensetive as vst and takes too long and fails)
## 2. vst() or varianceStabilizingTransformation() --> the first is just a wrapper for the second. 
## 3. normTransform() --> log2 transformation or any other desired fuction. Will not be used!
## 4. lfcShrink
################################################################################
## Data Normalization for plotting: 1. VST Normalization
dds.vsd <- varianceStabilizingTransformation(dds.run, blind=FALSE)
dds.vsd
# the matrix of transformed values is stored in assay(dds.vsd)

# saving to RDS object
saveRDS(dds.vsd, "saved_objects/dds.vsd.rds")
dds.vsd <- readRDS("saved_objects/dds.vsd.rds")

# Extracting all DEGs
dds.vsd.degs <- dds.vsd[rownames(dds.vsd) %in% rownames(res.degs),]
dds.vsd.degs

# Extracting all repair and repair DEGs
dds.vsd.rep <- dds.vsd[rownames(dds.vsd) %in% repair.genes$symbol,]
dds.vsd.rep.degs <- dds.vsd.rep[row.names(dds.vsd.rep) %in% rownames(res.rep.degs)]
################################################################################
## Data Normalization for plotting: 2. lfcShrink Normalization
resultsNames(dds.run)
res.lfc <- lfcShrink(dds.run, coef=2, res=res, type = "apeglm") # change the coef with resultsNames(dds.run)
summary(res.lfc)
# should we specify lfcThreshold and get s-value instead of p-value or not?

saveRDS(res.lfc, "saved_objects/res.lfc.rds")
res.lfc <- readRDS("saved_objects/res.lfc.rds")

## Chkeckin the distribution of the LFC and p-adjusted value
summary(res.lfc$log2FoldChange)
hist(res.lfc$log2FoldChange, main="Distribution of the Log2 fold change in the lfc-shrinked (apeglm) results", xlab="Log2 Fold Change")
summary(res.lfc$padj)
hist(res.lfc$padj, main="Distribution of the adjusted p-value in the lfc-shrinked (apeglm) results", xlab="Adjusted p-value")

## all DEGS
res.lfc.degs <- res.lfc[rownames(res.lfc) %in% rownames(res.degs),]
summary(res.lfc.degs)
summary(res.lfc.degs$log2FoldChange)
hist(res.lfc.degs$log2FoldChange, main="Distribution of the Log2 fold change in the lfc-shrinked All DEGs", xlab="Log2 Fold Change")
summary(res.lfc.degs$padj)
hist(res.lfc.degs$padj, main="Distribution of the adjusted p-value in the lfc-shrinked All DEGs", xlab="Adjusted p-value")


## repair genes
res.lfc.rep <- res.lfc[rownames(res.lfc) %in% repair.genes$symbol,]
summary(res.lfc.rep)

summary(res.lfc.rep$log2FoldChange)
hist(res.lfc.rep$log2FoldChange, main="Distribution of the Log2 fold change in the lfc-shrinked results (Repair Genes Only)", xlab="Log2 Fold Change")
summary(res.lfc.rep$padj)
hist(res.lfc.rep$padj, main="Distribution of the adjusted p-value in the lfc-shrinked results (Repair Genes Only)", xlab="Adjusted p-value")

## Repair genes DEGs
res.lfc.rep.degs <- res.lfc.rep[rownames(res.lfc.rep) %in% rownames(res.rep.degs),]
summary(res.lfc.rep.degs)

summary(res.lfc.rep.degs$log2FoldChange)
hist(res.lfc.rep.degs$log2FoldChange, main="Distribution of the Log2 fold change in the lfc-shrinked repair DEGs", xlab="Log2 Fold Change")
summary(res.lfc.rep.degs$padj)
hist(res.lfc.rep.degs$padj, main="Distribution of the adjusted p-value in the lfc-shrinked repair DEGs", xlab="Adjusted p-value")

################################################################################
## Data Normalization for plotting: 3. normTransform
# dds.nrd <- normTransform(dds.run)
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
# plotMA(res.rep.degs, alpha=0.05, main="MA plot of the not-normalized repair DEGs")
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

# normalized data
EnhancedVolcano(res.lfc, 
                lab = rownames(res.lfc), 
                x = 'log2FoldChange', 
                y = 'pvalue',
                title = "Volcano plot of the lfcschrink-normalized DESeq Results",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                legendPosition = "right")

# # what about the repair genes?
# EnhancedVolcano(res.rep, 
#                 lab = rownames(res.rep), 
#                 x = 'log2FoldChange', 
#                 y = 'pvalue',
#                 title = "Volcano plot of the not-normalized repair genes results",
#                 pCutoff = 0.05,
#                 FCcutoff = 1.0,
#                 legendPosition = "right")
# 
# EnhancedVolcano(res.lfc.rep, 
#                 lab = rownames(res.lfc.rep), 
#                 x = 'log2FoldChange', 
#                 y = 'pvalue',
#                 title = "Volcano plot of the lfcschrink-normalized repair genes results",
#                 pCutoff = 0.05,
#                 FCcutoff = 1.0,
#                 legendPosition = "right")
################################################################################
## Plotting: 3. PCA
# should we use normalized data instead? The manual says yes!
# what is the best ntop value? size of the DEGs?

# a. using varianceStabilizingTransformation
pca.plot.vsd <- plotPCA(dds.vsd, intgroup="Sample.Type", ntop = 25531, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd, "percentVar"))
ggplot(pca.plot.vsd, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vst-normalized DESeq DataSet (ntop = 25531)")
rm(percentVar)

# b. using normTransform
# pca.plot.nrd <- plotPCA(dds.nrd, intgroup="Sample.Type", ntop = 25531, returnData=TRUE)
# percentVar <- round(100 * attr(pca.plot.nrd, "percentVar"))
# ggplot(pca.plot.nrd, aes(PC1, PC2, color=Sample.Type)) + 
#   geom_point(size=1) + stat_ellipse(type = "norm") + 
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   # xlab(percentage[1]) + ylab(percentage[2]) +
#   ggtitle("PCA plot of the log2-normalized DESeq DataSet")
# rm(percentVar)

# what about the repair genes?
pca.plot.vsd.rep <- plotPCA(dds.vsd.rep, intgroup="Sample.Type", ntop = 285, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.rep, "percentVar"))
ggplot(pca.plot.vsd.rep, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vsd-normalized DESeq DataSet of the repair genes only (ntop = 285)")
rm(percentVar)
################################################################################
## Plotting: 4. Dispersion estimate
# estimateDispersions(dds)
plotDispEsts(dds.run)
################################################################################
## Plotting: 5. meanSdPlot
# meanSdPlot()
################################################################################
## Plotting: 6. Counts plot
genes.list <- c("PCLAF", "ISG15", "EXO1")
for (gen in genes.list){
  plotCounts(dds.run, gene=gen, intgroup = "Sample.Type") # returnData = FALSE
}
rm(genes.list)
plotCounts(dds.run, gene=gen, intgroup = "Sample.Type")
################################################################################
## Plotting: 7. Boxplot with dots
list(assays(dds.run)) # counts mu H cooks replaceCounts replaceCooks
boxplot(t(assays(dds.run["EXO1"])[["counts"]])~dds.run$Sample.Type, 
        range=0, las=1, log='y',
        xlab="Groups", ylab="Counts",
        col=c("darksalmon", "darkred"))
stripchart(t(assays(dds.run["EXO1"])[["counts"]])~dds.run$Sample.Type, 
           vertical=TRUE, method='jitter', add=TRUE, pch=16, col=c("firebrick", "orange")) 

list(assays(dds.vsd)) # counts mu H cooks replaceCounts replaceCooks
boxplot(t(assays(dds.vsd["EXO1"])[[1]])~dds.vsd$Sample.Type, 
        range=0, las=1, log='y',
        xlab="Groups", ylab="Counts", main="EXO1 (vsd-normalized)",
        col=c("darksalmon", "darkred"))
stripchart(t(assays(dds.vsd["EXO1"])[[1]])~dds.vsd$Sample.Type, 
           vertical=TRUE, method='jitter', add=TRUE, pch=16, col=c("firebrick", "orange")) 
################################################################################
## Plotting: 8. Heatmap & Dendogram
colors2 <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(5)
# dds.vsd.rep.degs.2 <- dds.vsd.rep.degs
# colnames(dds.vsd.rep.degs.2) <- substring(colnames(dds.vsd.rep.degs), 6)
heatmap.2(assay(dds.vsd.rep.degs), col=colors2,
          scale="row", trace="none", labCol=substring(colnames(dds.vsd.rep.degs), 6),
          # dendrogram = "row",
          Colv=order(dds.vsd.rep.degs$Sample.Type), Rowv=TRUE,
          ColSideColors = c(Primary.Tumor="darkgreen", Solid.Tissue.Normal="orange")[colData(dds.vsd.rep.degs)$Sample.Type],
          key =TRUE, key.title="Heatmap Key",
          main="Heat map of the differentially expressed repair genes")

################################################################################
################################################################################
### Saving all figures from the plots tab at once 
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="saved_objects/figures/")
################################################################################
