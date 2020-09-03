################################################################################
### miRNA IDC_Tumor vs. LC_Tumor

library(RColorBrewer)
library(gplots)
library(dplyr)

library(DESeq2)
library(EnhancedVolcano)

source("scripts/proteins.R")
source("scripts/plotting_gene_2.R")
source("scripts/to_tiff.R")

mirna.dds.run <- readRDS("saved_objects/mirna.dds.run.rds")
mirna.dds.vsd <- readRDS("saved_objects/mirna.dds.vsd.rds")
# hsa_MTI <- read.xlsx("raw_data/miRNA-DBs/hsa_MTI.xlsx")
# substring(hsa_MTI$miRNA, 7) = 'r'

mir.gene.pairs <- read.csv("saved_objects/mir.gene.pairs.csv", header=TRUE)
mir.gene.pairs.rep <- read.csv("saved_objects/mir.gene.pairs.rep.csv", header=TRUE)
mir.gene.pairs.merged <- read.csv("saved_objects/mir.gene.pairs.merged.csv", header=TRUE)
################################################################################
mirna.dds.run.idclc <- mirna.dds.run[,mirna.dds.run$group %in% c("IDC_Tumor", "LC_Tumor")]
mirna.dds.run.idclc
# 1881 miRs, 957 Samples

## Creating the mirna.res object from the original mirna.dds
mirna.res.idclc <- results(mirna.dds.run, contrast=c("group", "IDC_Tumor", "LC_Tumor"), alpha=0.05, lfcThreshold=1)
summary(mirna.res.idclc)
## out ot 1,598 mirs, 41 upregulated and 0 downregulated

## Distribution of the LFC and p-adjusted value
summary(mirna.res.idclc$log2FoldChange)
to_tiff(
  hist(mirna.res.idclc$log2FoldChange, main="Log2 Fold Change in miRs between IDC and LC", xlab="Log2 Fold Change"),
  "mirna.idclc.all.l2fc.tiff")
summary(mirna.res.idclc$padj)
to_tiff(
  hist(mirna.res.idclc$padj, main="Adjusted p-value in miRs between IDC and LC", xlab="Adjusted p-value"),
  "mirna.idclc.all.adjpval.tiff")
################################################################################
## complete cases only
mirna.res.idclc.complete <- mirna.res.idclc[complete.cases(mirna.res.idclc),]
################################################################################
## Extracting DEMs
mirna.res.idclc.dems <- mirna.res.idclc.complete[mirna.res.idclc.complete$padj<0.05 & abs(mirna.res.idclc.complete$log2FoldChange)>1,]
summary(mirna.res.idclc.dems)
# 41 DEMs: 41 up, 0 down

## Distribution of the p-adjusted value
summary(mirna.res.idclc.dems$padj)
to_tiff(
  hist(mirna.res.idclc.dems$padj, main="Adjusted p-value in DEMs between IDC and LC", xlab="Adjusted p-value"),
  "mirna.idclc.dems.adjpval.tiff")

################################################################################
## Exporting to file
write.table(as.data.frame(mirna.res.idclc.dems), "saved_objects/mirna.res.idclc.dems.csv", sep=",", quote=FALSE)
################################################################################
## Getting DEMs affecting any repair
mirna.res.idclc.reptar <- mirna.res.idclc.complete[rownames(mirna.res.idclc.complete) %in% mir.gene.pairs.rep$miRNA,]
summary(mirna.res.idclc.reptar)
# out of 116 miRs: 3 up and 0 down

mirna.res.idclc.dems.reptar <- mirna.res.idclc.dems[rownames(mirna.res.idclc.dems) %in% mir.gene.pairs.rep$miRNA,]
summary(mirna.res.idclc.dems.reptar)
# 3 DEMs affecting any Repair: 3 up, 0 down
################################################################################
## Adding the Gene Names from mirTarBase DB
mirna.res.idclc.reptar.df <- cbind(miRs = row.names(mirna.res.idclc.reptar), as.data.frame(mirna.res.idclc.reptar))
mirna.res.idclc.reptar.df <- left_join(mirna.res.idclc.reptar.df, mir.gene.pairs.merged,
                                     by=c("miRs" = "miRNA"))

mirna.res.idclc.dems.reptar.df <- cbind(miRs = row.names(mirna.res.idclc.dems.reptar), as.data.frame(mirna.res.idclc.dems.reptar))
mirna.res.idclc.dems.reptar.df <- left_join(mirna.res.idclc.dems.reptar.df, mir.gene.pairs.merged,
                                          by=c("miRs" = "miRNA"))

## Exporting to file
write.table(mirna.res.idclc.reptar.df, "saved_objects/mirna.res.idclc.reptar.withgenes.csv", sep=",", quote=FALSE, row.names = FALSE)
write.table(mirna.res.idclc.dems.reptar.df, "saved_objects/mirna.res.idclc.dems.reptar.withgenes.csv", sep=",", quote=FALSE, row.names = FALSE)

################################################################################
## Experimentally Validated Only (all miRs) table S6
mirna.res.idclc.dems.df <- cbind(miRs = rownames(mirna.res.idclc.dems), as.data.frame(mirna.res.idclc.dems))
mirna.res.idclc.dems.df <- left_join(mirna.res.idclc.dems.df, mir.gene.pairs.merged,
                                     by=c("miRs" = "miRNA"))
mirna.res.idclc.dems.df <- mirna.res.idclc.dems.df[complete.cases(mirna.res.idclc.dems.df),]
mirna.res.idclc.dems.reptar.df <- mirna.res.idclc.dems.df[mirna.res.idclc.dems.df$miRs %in% rownames(mirna.res.idclc.dems.reptar),]

write.table(mirna.res.idclc.dems.df, "saved_objects/S6-mirna.res.idclc.dems.withgenes.csv", sep=",", quote=FALSE, row.names = FALSE)
write.table(mirna.res.idclc.dems.reptar.df, "saved_objects/S6-mirna.res.idclc.dems.reptar.withgenes.csv", sep=",", quote=FALSE, row.names = FALSE)
################################################################################
## Subset VSD
## All miRs
mirna.dds.vsd.idclc <- mirna.dds.vsd[,mirna.dds.vsd$group %in% c("IDC_Tumor", "LC_Tumor")]
mirna.dds.vsd.idclc

## All miRs affecting any repair
mirna.dds.vsd.idclc.reptar <- mirna.dds.vsd.idclc[rownames(mirna.dds.vsd.idclc) %in% mir.gene.pairs.rep$miRNA,]
mirna.dds.vsd.idclc.reptar

## DEMs
mirna.dds.vsd.idclc.dems <- mirna.dds.vsd.idclc[rownames(mirna.dds.vsd.idclc) %in% rownames(mirna.res.idclc.dems),]
mirna.dds.vsd.idclc.dems

## DEMs affecting any repair
mirna.dds.vsd.idclc.dems.reptar <- mirna.dds.vsd.idclc.dems[rownames(mirna.dds.vsd.idclc.dems) %in% mir.gene.pairs.rep$miRNA,]
mirna.dds.vsd.idclc.dems.reptar
################################################################################
## Plotting: 1. MAplot
to_tiff(plotMA(mirna.res.idclc, alpha=0.05, main="MAplot of the miRs between IDC and LC"), 
        "mirna.idclc.all.plotMA.tiff")
################################################################################
## Plotting: 2. EnhancedVolcano Plot
tiff("final_figures/mirna.idclc.all.volcano.tiff", width = 18, height = 21, units = 'cm', res = 300)
EnhancedVolcano(mirna.res.idclc, 
                lab = rownames(mirna.res.idclc), 
                x = 'log2FoldChange', 
                y = 'padj',
                title = "All miRs between IDC and LC",
                pCutoff = 0.05,
                FCcutoff = 1.0)
dev.off()

tiff("final_figures/mirna.idclc.rep.volcano.tiff", width = 18, height = 21, units = 'cm', res = 300)
EnhancedVolcano(mirna.res.idclc.reptar, 
                lab = rownames(mirna.res.idclc.reptar), 
                x = 'log2FoldChange', 
                y = 'padj',
                title = "miRs of Repair Genes between IDC and LC",
                pCutoff = 0.05,
                FCcutoff = 1.0)
dev.off()
################################################################################
## Plotting: 3. PCA
# PCA all miRs
tiff("final_figures/mirna.idclc.pca.all.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca.idclc <- plotPCA(mirna.dds.vsd.idclc, intgroup="group", ntop = 1881, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca.idclc, "percentVar"))
ggplot(mirna.pca.idclc, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of all miRs (1,881) between IDC and LC")
rm(mirna.pca.idclc, percentVar)
dev.off()

# PCA all dems
tiff("final_figures/mirna.idclc.pca.dems.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca.idclc.dems <- plotPCA(mirna.dds.vsd.idclc.dems, intgroup="group", ntop = 41, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca.idclc.dems, "percentVar"))
ggplot(mirna.pca.idclc.dems, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of all DEMs (41) between IDC and LC")
rm(mirna.pca.idclc.dems, percentVar)
dev.off()

# PCA of all miRs affecting any Repair
tiff("final_figures/mirna.idclc.pca.all.reptar.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca.idclc.reptar <- plotPCA(mirna.dds.vsd.idclc.reptar, intgroup="group", ntop = 586, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca.idclc.reptar, "percentVar"))
ggplot(mirna.pca.idclc.reptar, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of all miRs between IDC and LC Targetiing Repair Genes (586)")
rm(mirna.pca.idclc.reptar, percentVar)
dev.off()

# PCA of DEMs affecting any Repair
tiff("final_figures/mirna.idclc.pca.dems.reptar.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca.idclc.dems.reptar <- plotPCA(mirna.dds.vsd.idclc.dems.reptar, intgroup="group", ntop = 3, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca.idclc.dems.reptar, "percentVar"))
ggplot(mirna.pca.idclc.dems.reptar, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of DEMs between IDC and LC Targetiing Repair Genes (3)")
rm(mirna.pca.idclc.dems.reptar, percentVar)
dev.off()
################################################################################
## Plotting: 8. Heatmap & Dendogram
colors2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)

pdf(file="final_figures/mirna.idclc.heatmap.dems.pdf", onefile=TRUE, paper="a4r", width = 11, height = 8)
heatmap.2(assay(mirna.dds.vsd.idclc.dems), col=colors2,
          scale="row", trace="none", labCol=substring(colnames(mirna.dds.vsd.idclc.dems), 6),
          Colv=order(mirna.dds.vsd.idclc.dems$Sample.Type), Rowv=TRUE,
          ColSideColors = c(IDC="darkgreen", LC="orange")[colData(mirna.dds.vsd.idclc.dems)$primary_diagnosis],
          key =TRUE, key.title="Heatmap Key",
          main="41 DEMs between IDC (darkgreen) and LC (orange)")
dev.off()

pdf(file="final_figures/mirna.idclc.heatmap.dems.reptar.pdf", onefile=TRUE, paper="a4r", width = 11, height = 8)
heatmap.2(assay(mirna.dds.vsd.idclc.dems.reptar), col=colors2,
          scale="row", trace="none", labCol=substring(colnames(mirna.dds.vsd.idclc.dems), 6),
          Colv=order(mirna.dds.vsd.idclc.dems$Sample.Type), Rowv=TRUE,
          ColSideColors = c(IDC="darkgreen", LC="orange")[colData(mirna.dds.vsd.idclc.dems)$primary_diagnosis],
          key =TRUE, key.title="Heatmap Key",
          main="3 DEMs between IDC (darkgreen) and LC (orange) Affecting Repair Genes")
dev.off()
################################################################################
## Plotting: 7. Boxplot with dots
dem.list <- rownames(mirna.res.idclc.dems)
pdf(file="final_figures/mirna.idclc.dems.boxplot.pdf", onefile=TRUE, paper="a4", width = 8, height = 11)
par(mfrow=c(3,3))
for (dem in dem.list){
  plotting_gene_2(mirna.dds.vsd.idclc, dem, "", atr=1)
}
dev.off()
rm(dem, dem.list)

dem.list <- rownames(mirna.res.idclc.dems.reptar)
pdf(file="final_figures/mirna.idclc.dems.reptar.boxplot.pdf", onefile=TRUE, paper="a4", width = 8, height = 11)
par(mfrow=c(3,3))
for (dem in dem.list){
  plotting_gene_2(mirna.dds.vsd.idclc, dem, "", atr=1)
}
dev.off()
rm(dem, dem.list)

par(mfrow=c(1,1))
################################################################################
################################################################################
## mirTarBase (experimentally validated only)
## subsetting the DB by the DEMs
dems.mirtarbase.idclc <- mir.gene.pairs[mir.gene.pairs$miRNA %in% rownames(mirna.res.idclc.dems),]

## Exploring
length(unique(dems.mirtarbase.idclc$miRNA)) 
# 7 mir out of 41
length(unique(dems.mirtarbase.idclc$Target.Gene)) 
# targeting total of 705 genes
sum(unique(dems.mirtarbase.idclc$Target.Gene) %in% repair.genes$symbol) 
# 15 of the 705 genes are repair
################################################################################