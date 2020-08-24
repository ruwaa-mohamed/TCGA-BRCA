################################################################################
### miRNA IDC_Tumor vs. LC_Tumor

library(RColorBrewer)
library(gplots)
library(ggplot2)
library(openxlsx)

library(DESeq2)
library(EnhancedVolcano)

source("scripts/proteins.R")
source("scripts/plotting_gene.R")
source("scripts/to_tiff.R")

mirna.dds.run <- readRDS("saved_objects/mirna.dds.run.rds")
mirna.dds.vsd <- readRDS("saved_objects/mirna.dds.vsd.rds")
hsa_MTI <- read.xlsx("raw_data/miRNA-DBs/hsa_MTI.xlsx")
substring(hsa_MTI$miRNA, 7) = 'r'

mir.gene.pairs <- read.csv("saved_objects/mir.gene.pairs.csv", header=TRUE)
mir.gene.pairs.rep <- read.csv("saved_objects/mir.gene.pairs.rep.csv", header=TRUE)
################################################################################
mirna.dds.run.idclc <- mirna.dds.run[,mirna.dds.run$group %in% c("IDC_Tumor", "LC_Tumor")]
mirna.dds.run.idclc

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

## Distribution of the p-adjusted value
summary(mirna.res.idclc.dems$padj)
to_tiff(
  hist(mirna.res.idclc.dems$padj, main="Adjusted p-value in DEMs between IDC and LC", xlab="Adjusted p-value"),
  "mirna.idclc.dems.adjpval.tiff")

################################################################################
## Exporting to file
write.table(as.data.frame(mirna.res.idclc.dems), "saved_objects/mirna.res.idclc.dems.csv", sep=",", quote=FALSE)
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
                y = 'pvalue',
                title = "Volcano plot of the miRs between IDC and LC",
                pCutoff = 0.05,
                FCcutoff = 1.0)
dev.off()
################################################################################
## Plotting: 3. PCA

## subset vsd
mirna.dds.vsd.idclc <- mirna.dds.vsd[,mirna.dds.vsd$group %in% c("IDC_Tumor", "LC_Tumor")]
mirna.dds.vsd.idclc

## DEMs
mirna.dds.vsd.idclc.dems <- mirna.dds.vsd.idclc[rownames(mirna.dds.vsd.idclc) %in% rownames(mirna.res.idclc.dems),]
mirna.dds.vsd.idclc.dems

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

################################################################################
## Plotting: 8. Heatmap & Dendogram
pdf(file="final_figures/mirna.idclc.heatmap.dems.pdf", onefile=TRUE, paper="a4r", width = 11, height = 8)
colors2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)
heatmap.2(assay(mirna.dds.vsd.idclc.dems), col=colors2,
          scale="row", trace="none", labCol=substring(colnames(mirna.dds.vsd.idclc.dems), 6),
          Colv=order(mirna.dds.vsd.idclc.dems$Sample.Type), Rowv=TRUE,
          ColSideColors = c(IDC="darkgreen", LC="orange")[colData(mirna.dds.vsd.idclc.dems)$primary_diagnosis],
          key =TRUE, key.title="Heatmap Key",
          main="Heatmap of the 41 DEMs between IDC (darkgreen) and LC (orange)")
dev.off()
################################################################################
## Plotting: 7. Boxplot with dots
dem.list <- rownames(mirna.res.idclc.dems)

pdf(file="final_figures/mirna.idclc.dems.boxplot.pdf", onefile=TRUE, paper="a4", width = 8, height = 11)
par(mfrow=c(3,3))
for (dem in dem.list){
  plotting_gene(mirna.dds.run.idclc, dem, "")
}
dev.off()

par(mfrow=c(1,1))
rm(dem, dem.list)
################################################################################
################################################################################
## mirTarBase (experimentally validated only)
## subsetting the DB by the DEMs
hsa_MTI.idclc.dems <- hsa_MTI[hsa_MTI$miRNA %in% rownames(mirna.res.idclc.dems),]

## Exploring
length(unique(hsa_MTI.idclc.dems$miRNA)) 
# 7 mir out of 41
length(unique(hsa_MTI.idclc.dems$Target.Gene)) 
# targeting total of 705 genes
sum(unique(hsa_MTI.idclc.dems$Target.Gene) %in% repair.genes$symbol) 
# 15 of the 705 genes are repair
################################################################################