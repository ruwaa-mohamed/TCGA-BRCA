################################################################################
### IDC_Tumor vs. LC_Tumor

library(gplots)
library(ggplot2)
library(RColorBrewer)

library(DESeq2)
library(EnhancedVolcano)

source("scripts/proteins.R")
source("scripts/plotting_gene.R")
source("scripts/to_tiff.R")

dds.run <- readRDS("saved_objects/dds.run.group.rds")
dds.vsd <- readRDS("saved_objects/dds.vsd.rds")
################################################################################
## Extracting dds.run.idclc from dds.run
dds.run.idclc <- dds.run[,dds.run$group %in% c("IDC_Tumor", "LC_Tumor")]
dds.run.idclc
dds.run.idclc$group <- droplevels(dds.run.idclc$group)
table(dds.run.idclc$group)
################################################################################
## Creating the results (RES) object from the original dds.run
res.idclc <- results(dds.run, contrast=c("group", "IDC_Tumor", "LC_Tumor"), alpha=0.05, lfcThreshold=1)
summary(res.idclc)
# out of 25173 gens: 591 upregulated and 72 downregulated

## Distribution of the LFC and p-adjusted value
summary(res.idclc$log2FoldChange)
to_tiff(
  hist(res.idclc$log2FoldChange, main="Log2 Fold Change of the Genes between IDC and LC", xlab="Log2 Fold Change"),
  "rna.idclc.all.l2fc.hist.tiff")
summary(res.idclc$padj)
to_tiff(
  hist(res.idclc$padj, main="Adjusted p-value of the Genes between IDC and LC", xlab="Adjusted p-value"),
  "rna.idclc.all.adjpval.hist.tiff")

## complete cases only
res.idclc.complete <- res.idclc[complete.cases(res.idclc),]
################################################################################
## Extracting DEGs
res.idclc.degs <- res.idclc.complete[res.idclc.complete$padj<0.05 & abs(res.idclc.complete$log2FoldChange)>1,]
summary(res.idclc.degs)

## Distribution of the p-adjusted value
summary(res.idclc.degs$padj)
to_tiff(hist(
  res.idclc.degs$padj, main="Adjusted p-value of the DEGs between IDC and LC", xlab="Adjusted p-value"),
  "rna.idclc.degs.adjpval.hist.tiff")
################################################################################
## Repair genes subsetting 
res.idclc.rep <- res.idclc.complete[rownames(res.idclc.complete) %in% repair.genes$symbol,]
summary(res.idclc.rep)
## out of 285 genes, 0 up-regulated, 0 down-regulated
## There's no need for any of the rep.degs in the rest of the script.


## Distribution of the LFC and p-adjusted value
summary(res.idclc.rep$log2FoldChange)
to_tiff(
  hist(res.idclc.rep$log2FoldChange, main="Log2 Fold Change of All Repair Genes between IDC and LC", xlab="Log2 Fold Change"),
  "rna.idclc.rep.all.l2fc.hist.tiff")
################################################################################
## Exporting to files
write.table(as.data.frame(res.idclc.degs), "saved_objects/res.idclc.degs.csv", sep=",", quote=FALSE)
################################################################################
## Plotting: 1. MAplot
to_tiff(
  plotMA(res.idclc, alpha=0.05, main="MAplot of All Genes between IDC and LC"),
  "rna.idclc.all.MAplot.tiff")
################################################################################
## Plotting: 2. EnhancedVolcano Plot
tiff("final_figures/rna.idclc.EnhancedVolcano.tiff", width = 18, height = 21, units = 'cm', res = 300)
EnhancedVolcano(res.idclc, 
                lab = rownames(res.idclc), 
                x = 'log2FoldChange', 
                y = 'pvalue',
                title = "All Genes between IDC and LC",
                pCutoff = 0.05,
                FCcutoff = 1.0)
dev.off()
################################################################################
## subset vsd
dds.vsd.idclc <- dds.vsd[,dds.vsd$group %in% c("IDC_Tumor", "LC_Tumor")]
dds.vsd.idclc
dds.vsd.idclc$group <- droplevels(dds.vsd.idclc$group)
table(dds.vsd.idclc$group)

## DEGs
dds.vsd.idclc.degs <- dds.vsd.idclc[rownames(dds.vsd.idclc) %in% rownames(res.idclc.degs),]
dds.vsd.idclc.degs
################################################################################
## Plotting: 3. PCA
# all genes
tiff("final_figures/rna.idclc.pca.all.tiff", width = 18, height = 21, units = 'cm', res = 300)
pca.plot.vsd.idclc <- plotPCA(dds.vsd.idclc, intgroup="group", ntop = 25531, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.idclc, "percentVar"))
ggplot(pca.plot.vsd.idclc, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of All (25,531) Genes between IDC and LC")
rm(percentVar)
dev.off()

# all degs
tiff("final_figures/rna.idclc.pca.degs.tiff", width = 18, height = 21, units = 'cm', res = 300)
pca.plot.vsd.idclc.degs <- plotPCA(dds.vsd.idclc.degs, intgroup="group", ntop = 663, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.idclc.degs, "percentVar"))
ggplot(pca.plot.vsd.idclc.degs, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of all DEGs (663) between IDC and LC")
rm(percentVar)
dev.off()
################################################################################
## Plotting: 4. Dispersion estimate
to_tiff(
  plotDispEsts(dds.run.idclc, main="Disperssion Estimates for all Genes between IDC and LC"),
  "rna.idclc.plotDispEsts.tiff")
################################################################################