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
source('scripts/to_tiff.R')

library(RColorBrewer)
library(gplots)
library(ggplot2)
library(openxlsx)

library(DESeq2)
library(EnhancedVolcano)
################################################################################
## Reading the sample sheet 
## (after adding clinical data and sync with RNA-seq sample.sheet)
mirna.sample.sheet <- read.csv('saved_objects/mirna.sample.sheet.csv', header=TRUE)

## As factor and Releveling
mirna.sample.sheet$Sample.Type <- as.factor(mirna.sample.sheet$Sample.Type)
levels(mirna.sample.sheet$Sample.Type)
table(mirna.sample.sheet$Sample.Type)
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
mirna.col.data$primary_diagnosis <- as.factor(mirna.col.data$primary_diagnosis)
levels(mirna.col.data$primary_diagnosis)

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
mirna.dds <- DESeqDataSetFromMatrix(countData=mirna.count.data , colData=mirna.col.data , design=~Sample.Type)
mirna.dds
saveRDS(mirna.dds, "saved_objects/mirna.dds.rds")
mirna.dds <- readRDS("saved_objects/mirna.dds.rds")

## Collapse Technical Replicates
mirna.ddsCollapsed <- collapseReplicates(mirna.dds, groupby=mirna.dds$Sample.ID)
mirna.ddsCollapsed
saveRDS(mirna.ddsCollapsed, "saved_objects/mirna.ddsCollapsed.rds")
mirna.ddsCollapsed <- readRDS("saved_objects/mirna.ddsCollapsed.rds")

## to check that collapsing worked!
original <- rowSums(counts(mirna.dds)[, mirna.dds$Sample.ID == "TCGA-A7-A0DB-01A"])
all(original == counts(mirna.ddsCollapsed)[,"TCGA-A7-A0DB-01A"])
rm(original)

## organizing two studies (All and subtypes)
mirna.ddsCollapsed.TN <- mirna.ddsCollapsed
design(mirna.ddsCollapsed.TN)

mirna.ddsCollapsed$group <- factor(paste(mirna.ddsCollapsed$primary_diagnosis, mirna.ddsCollapsed$Sample.Type, sep="_"))
design(mirna.ddsCollapsed) <- ~ group
design(mirna.ddsCollapsed)
table(mirna.ddsCollapsed$group)
# IDC_Normal    IDC_Tumor    LC_Normal     LC_Tumor Mixed_Normal  Mixed_Tumor 
# 82          756            6          201            8           27
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
# -- replacing outliers and refitting for 89 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
mirna.dds.run

saveRDS(mirna.dds.run, file="saved_objects/mirna.dds.run.rds")
mirna.dds.run <- readRDS("saved_objects/mirna.dds.run.rds")

mirna.dds.vsd <- varianceStabilizingTransformation(mirna.dds.run, blind=FALSE)
mirna.dds.vsd
################################################################################
################################################################################
## Run DESEQ2 (ALL)
mirna.dds.run.TN <- DESeq(mirna.ddsCollapsed.TN)
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
mirna.dds.run.TN

saveRDS(mirna.dds.run.TN, file="saved_objects/mirna.dds.run.TN.rds")
mirna.dds.run.TN <- readRDS("saved_objects/mirna.dds.run.TN.rds")
################################################################################
## Creating the results (RES) object (ALL)
mirna.res.TN <- results(mirna.dds.run.TN, contrast=c("Sample.Type", "Tumor", "Normal"), alpha=0.05,lfcThreshold=1)
summary(mirna.res.TN)
# out of 1599 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 1.00 (up)    : 86, 5.4%
# LFC < -1.00 (down) : 61, 3.8%
# outliers [1]       : 0, 0%
# low counts [2]     : 614, 38%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

## saving to and reading from RDS object
saveRDS(mirna.res.TN, file="saved_objects/mirna.res.TN.rds")
mirna.res.TN <- readRDS("saved_objects/mirna.res.TN.rds")

## Removing incomplete genes
mirna.res.TN.complete <- mirna.res.TN[complete.cases(mirna.res.TN),]
summary(mirna.res.TN.complete)

## Chkeckin the distribution of the LFC and p-adjusted value
summary(mirna.res.TN$log2FoldChange)
to_tiff(
  hist(mirna.res.TN$log2FoldChange, main="Distribution of the Log2 Fold Change in miRNA results", xlab="Log2 Fold Change"), 
  "miRNA-TN-Hist-L2FC.tiff")
summary(mirna.res.TN$padj)
to_tiff(
  hist(mirna.res.TN$padj, main="Distribution of the adjusted p-value in the miRNA results", xlab="Adjusted p-value"),
  "miRNA-TN-Hist-Adjusted-pval.tiff")
################################################################################
## Getting the DEMs
mirna.res.TN.dems <- mirna.res.TN.complete[abs(mirna.res.TN.complete$log2FoldChange)>1 & mirna.res.TN.complete$padj<0.05,]
summary(mirna.res.TN.dems)

summary(mirna.res.TN.dems$padj)
to_tiff(
  hist(mirna.res.TN.dems$padj, main="Distribution of the adjusted p-value in the miRNA DEGs", xlab="Adjusted p-value"),
  "miRNA-TN-Hist-dems-Adjusted-pval.tiff")

write.csv(mirna.res.TN.dems, file="saved_objects/mirna.dems.all.csv", quote=FALSE)
write.table(rownames(mirna.res.TN.dems), file="saved_objects/mirna.res.TN.dems.tsv", row.names = FALSE, col.names=FALSE, quote=FALSE)

## as data frame
# mirna.res.TN.dems.df <- as.data.frame(mirna.res.TN.dems)
# mirna.res.TN.dems.df <- mirna.res.TN.dems.df[order(mirna.res.TN.dems.df$padj),]
################################################################################
## Data Normalization for plotting: 1. VST Normalization
mirna.dds.TN.vsd <- varianceStabilizingTransformation(mirna.dds.run.TN, blind=FALSE)
mirna.dds.TN.vsd

## saving to and reading from RDS object
saveRDS(mirna.dds.TN.vsd, file="saved_objects/mirna.dds.TN.vsd.rds")
mirna.dds.TN.vsd <- readRDS("saved_objects/mirna.dds.TN.vsd.rds")

## VSD for DEMs
mirna.dds.TN.vsd.dems <- mirna.dds.TN.vsd[rownames(mirna.dds.TN.vsd) %in% rownames(mirna.res.TN.dems),]
mirna.dds.TN.vsd.dems
################################################################################
## Plotting: 1. MA plot
to_tiff(
  plotMA(mirna.res.TN, alpha=0.05, main="MA plot of miRNA Results"),
  "mirna-MAplot.tiff")
################################################################################
## Plotting: 2. EnhancedVolcano
tiff("saved_objects/figures/Jul22/mirna-volcano.tiff", width = 18, height = 21, units = 'cm', res = 300)
EnhancedVolcano(mirna.res.TN, 
                        lab = rownames(mirna.res.TN), 
                        x = 'log2FoldChange', 
                        y = 'pvalue',
                        title = "Volcano plot of all miRs",
                        pCutoff = 0.05,
                        FCcutoff = 1.0)
dev.off()
################################################################################
## Plotting: 3. PCA
## All genes
tiff("saved_objects/figures/Jul22/mirna-pca-all.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca <- plotPCA(mirna.dds.TN.vsd, intgroup="Sample.Type", ntop = 1881, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca, "percentVar"))
ggplot(mirna.pca, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vsd-normalized DESeq DataSet of the miRNA data (ntop = 1881)")
rm(mirna.pca, percentVar)
dev.off()

## All DEMs
tiff("saved_objects/figures/Jul22/mirna-pca-dems.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca.dems <- plotPCA(mirna.dds.TN.vsd.dems, intgroup="Sample.Type", ntop = 147, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca.dems, "percentVar"))
ggplot(mirna.pca.dems, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vsd-normalized miRNA DEMs (ntop = 147)")
rm(mirna.pca.dems, percentVar)
dev.off()
################################################################################
## Plotting: 4. Dispersion estimate
to_tiff(plotDispEsts(mirna.dds.run.TN), "mirna-plotDispEsts.tiff")
################################################################################
## Plotting: 7. Boxplot with dots
dem.list <- rownames(mirna.res.TN.dems)

par(mfrow=c(3,3))
for (dem in dem.list){
  plotting_gene(mirna.dds.run.TN, dem)
}
par(mfrow=c(1,1))
rm(dem, dem.list)
################################################################################
## Plotting: 8. Heatmap
colors2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)
heatmap.2(assay(mirna.dds.TN.vsd.dems), col=colors2,
          scale="row", trace="none", labCol=substring(colnames(mirna.dds.TN.vsd.dems), 6),
          Colv=order(mirna.dds.TN.vsd.dems$Sample.Type), Rowv=TRUE,
          ColSideColors = c(Tumor="darkgreen", Normal="orange")[colData(mirna.dds.TN.vsd.dems)$Sample.Type],
          key =TRUE, key.title="Heatmap Key",
          main="Heatmap of the 147 differentially expressed miRNAs")
################################################################################
################################################################################
## mirTarBase (experimentally validated only)
## Reading the database
hsa_MTI <- read.xlsx("raw_data/miRNA-DBs/hsa_MTI.xlsx")
substring(hsa_MTI$miRNA, 7) = 'r'

## subsetting the DB by the DEMs
hsa_MTI.dems <- hsa_MTI[hsa_MTI$miRNA %in% rownames(mirna.res.TN.dems),]
# rm(hsa_MTI)

## Exploring
length(unique(hsa_MTI.dems$miRNA)) 
# 29 mir out of 147
length(unique(hsa_MTI.dems$Target.Gene)) 
# targeting total of 3,105 genes
sum(unique(hsa_MTI.dems$Target.Gene) %in% repair.genes$symbol) 
# 70 of the 3,105 genes are repair
################################################################################
# ## VennDiagram
# venn.diagram(x=list(rownames(res.degs), repair.genes$symbol, rownames(res.rep.degs), unique(hsa_MTI.dems$Target.Gene), unique(hsa_MTI.dems.rep_degs$Target.Gene)),
#              category.names=c("All DEGs", "All Repair genes", "Repair DEGs", "All genes affected by any of the DEMs", "Repair Genes affected by DEMs"),
#              filename="saved_objects/figures/venn_diagram.tiff",
#              output=TRUE,
#              fill = brewer.pal(5, "Pastel2"))
################################################################################
################################################################################