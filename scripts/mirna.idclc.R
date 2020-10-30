################################################################################
### miRNA IDC_Tumor vs. LC_Tumor

library(RColorBrewer)
library(gplots)
library(dplyr)

library(DESeq2)
library(EnhancedVolcano)

source("scripts/proteins.R")
source("scripts/plotting_gene.R")
source("scripts/to_tiff.R")

mirna.dds.run <- readRDS("saved_objects/mirna.dds.run.rds")
mirna.dds.vsd <- readRDS("saved_objects/mirna.dds.vsd.rds")

mir.gene.pairs <- read.csv("saved_objects/mir.gene.pairs.csv", header=TRUE)
mir.gene.pairs.rep <- read.csv("saved_objects/mir.gene.pairs.rep.csv", header=TRUE)
mir.gene.pairs.merged <- read.csv("saved_objects/mir.gene.pairs.merged.csv", header=TRUE)
################################################################################
mirna.dds.run.idclc <- mirna.dds.run[,mirna.dds.run$group %in% c("IDC_Tumor", "LC_Tumor")]
mirna.dds.run.idclc
# 1881 miRs, 957 Samples
mirna.dds.run.idclc$group <- droplevels(mirna.dds.run.idclc$group)
levels(mirna.dds.run.idclc$group) <- c("IDC", "LC")
table(mirna.dds.run.idclc$group)

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
write.table(as.data.frame(mirna.res.idclc.dems), "saved_objects/mirna.res.idclc.dems.all.csv", sep=",", quote=FALSE)
################################################################################
## Validated Only from mirTarBase DB
## All miRs
mirna.res.idclc.val <- mirna.res.idclc.complete[rownames(mirna.res.idclc.complete) %in% mir.gene.pairs$miRNA,]
summary(mirna.res.idclc.val)
# only 159 out of 764 complete miR (mirna.res.idclc.complete) 

## All DEMs
mirna.res.idclc.dems.val <- mirna.res.idclc.dems[rownames(mirna.res.idclc.dems) %in% mir.gene.pairs$miRNA,]
summary(mirna.res.idclc.dems.val)
# only 7 out of 40 DEMs 

## miRs affecting any repair
mirna.res.idclc.val.reptar <- mirna.res.idclc.val[rownames(mirna.res.idclc.val) %in% mir.gene.pairs.rep$miRNA,]
summary(mirna.res.idclc.val.reptar)
# 116 miRs (out of 159 validated miR): 3 up and 0 down

## DEMs affecting any repair
mirna.res.idclc.dems.val.reptar <- mirna.res.idclc.dems[rownames(mirna.res.idclc.dems) %in% mir.gene.pairs.rep$miRNA,]
summary(mirna.res.idclc.dems.val.reptar)
# 3 DEMs (out of 7 validated DEMs): 3 up and 0 down
################################################################################
## Adding Target Genes part (for Export)
## All DEMs with Target Genes
mirna.res.idclc.dems.val.df <- cbind(miRs = row.names(mirna.res.idclc.dems.val), as.data.frame(mirna.res.idclc.dems.val))
mirna.res.idclc.dems.val.df <- left_join(mirna.res.idclc.dems.val.df, mir.gene.pairs.merged, by=c("miRs" = "miRNA"))
mirna.res.idclc.dems.val.df <- mirna.res.idclc.dems.val.df [, c("miRs", "baseMean", "log2FoldChange", "pvalue", "padj", "Target.Gene")]

## Exporting to file
write.table(mirna.res.idclc.dems.val.df, "saved_objects/S6-mirna.res.idclc.dems.val.csv", sep=",", quote=FALSE, row.names = FALSE)
################################################################################
################################################################################
## Subsetting VSD
## All Validated miRs in IDC vs. LC
mirna.dds.vsd.idclc.val <- mirna.dds.vsd[rownames(mirna.dds.vsd) %in% rownames(mirna.res.idclc.val), mirna.dds.vsd$group %in% c("IDC_Tumor", "LC_Tumor")]
mirna.dds.vsd.idclc.val
mirna.dds.vsd.idclc.val$group <- droplevels(mirna.dds.vsd.idclc.val$group)
levels(mirna.dds.vsd.idclc.val$group) <- c("IDC", "LC")
table(mirna.dds.vsd.idclc.val$group)

## Validated miRs Targetting any Repair Genes
mirna.dds.vsd.idclc.val.reptar <- mirna.dds.vsd.idclc.val[rownames(mirna.dds.vsd.idclc.val) %in% mir.gene.pairs.rep$miRNA,]
mirna.dds.vsd.idclc.val.reptar

## All Validated DEMs
mirna.dds.vsd.idclc.val.dems <- mirna.dds.vsd.idclc.val[rownames(mirna.dds.vsd.idclc.val) %in% rownames(mirna.res.idclc.dems.val),]
mirna.dds.vsd.idclc.val.dems

## Validated DEMs Targetting any Repair Genes
mirna.dds.vsd.idclc.val.dems.reptar <- mirna.dds.vsd.idclc.val.dems[rownames(mirna.dds.vsd.idclc.val.dems) %in% mir.gene.pairs.rep$miRNA,]
mirna.dds.vsd.idclc.val.dems.reptar
################################################################################## Plotting: 1. MAplot
to_tiff(plotMA(mirna.res.idclc.val, alpha=0.05, main="MAplot of the miRs between IDC and LC"), 
        "mirna/mirna.idclc.all.val.MAplot.tiff")
################################################################################
## Plotting: 2. EnhancedVolcano Plot
tiff("final_figures/mirna/mirna.idclc.val.all.volcano.tiff", width = 18, height = 21, units = 'cm', res = 300)
EnhancedVolcano(mirna.res.idclc.val, 
                lab = rownames(mirna.res.idclc.val), 
                x = 'log2FoldChange', 
                y = 'padj',
                title = "All miRs between IDC and LC",
                pCutoff = 0.05,
                FCcutoff = 1.0)
dev.off()

tiff("final_figures/mirna/mirna.idclc.val.reptar.volcano.tiff", width = 18, height = 21, units = 'cm', res = 300)
EnhancedVolcano(mirna.res.idclc.val.reptar, 
                lab = rownames(mirna.res.idclc.val.reptar), 
                x = 'log2FoldChange', 
                y = 'padj',
                title = "miRs of Repair Genes between IDC and LC",
                pCutoff = 0.05,
                FCcutoff = 1.0)
dev.off()
################################################################################
## Plotting: 3. PCA
## All miRs with Validated Target Genes
tiff("final_figures/mirna/mirna.idclc.val.mirs.all.pca.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca.idclc.val <- plotPCA(mirna.dds.vsd.idclc.val, intgroup="group", ntop = 159, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca.idclc.val, "percentVar"))
ggplot(mirna.pca.idclc.val, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("All Validated miRs (159) between IDC and LC")
rm(mirna.pca.idclc.val, percentVar)
dev.off()

## All DEMs with Validated Target Genes 
tiff("final_figures/mirna/mirna.idclc.val.dems.all.pca.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca.idclc.val.dems <- plotPCA(mirna.dds.vsd.idclc.val.dems, intgroup="group", ntop = 7, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca.idclc.val.dems, "percentVar"))
ggplot(mirna.pca.idclc.val.dems, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("All Validated DEMs (7) between IDC and LC")
rm(mirna.pca.idclc.val.dems, percentVar)
dev.off()

## All miRs with Validated Target Repair Genes
tiff("final_figures/mirna/mirna.idclc.val.mirs.reptar.pca.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca.idclc.val.reptar <- plotPCA(mirna.dds.vsd.idclc.val.reptar, intgroup="group", ntop = 116, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca.idclc.val.reptar, "percentVar"))
ggplot(mirna.pca.idclc.val.reptar, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Validated Repair miRs (116) between IDC and LC")
rm(mirna.pca.idclc.val.reptar, percentVar)
dev.off()

## All DEMs with Validated Target Repair Genes 
tiff("final_figures/mirna/mirna.idclc.val.dems.reptar.pca.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca.idclc.val.dems.reptar <- plotPCA(mirna.dds.vsd.idclc.val.dems.reptar, intgroup="group", ntop = 3, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca.idclc.val.dems.reptar, "percentVar"))
ggplot(mirna.pca.idclc.val.dems.reptar, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Validated Repair DEMs (3) between IDC and LC")
rm(mirna.pca.idclc.val.dems.reptar, percentVar)
dev.off()
################################################################################
## Plotting: 8. Heatmap & Dendogram
colors2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)
plt2c <- c("#92CD2D", "#FC9C30", "#AA5EAC")

## All DEMs with Validated Target Genes 
tiff("final_figures/mirna/mirna.idclc.val.dems.all.heatmap.tiff", width = 25, height = 12, units = 'cm', res = 300)
heatmap.2(assay(mirna.dds.vsd.idclc.val.dems), col=colors2,
          scale="row", trace="none", labCol="", cexRow = 1.2, margins = c(1,8),
          Colv=order(mirna.dds.vsd.idclc.val.dems$group), Rowv=TRUE,
          ColSideColors = c(plt2c[c(2, 3)])[colData(mirna.dds.vsd.idclc.val.dems)$group],
          key =TRUE, key.title="Heatmap Key", lhei = c(1, 1.5),
          main="7 Validated DEMs between IDC and LC")
legend("topleft", legend = levels(mirna.dds.vsd.idclc.val.dems$group), col = plt2c[c(2,3)],
       lty= 1, lwd = 5, cex = 0.8, bty="0", bg="#FEFCF6",
       title="Samples", text.font=4, inset = c(0, 0.17), xjust = 0.5,
       text.width = 0.09)
dev.off()

## Validated DEMs with Target Repair Genes 
tiff("final_figures/mirna/mirna.idclc.val.dems.reptar.heatmap.tiff", width = 25, height = 8, units = 'cm', res = 300)
heatmap.2(assay(mirna.dds.vsd.idclc.val.dems.reptar), col=colors2,
          scale="row", trace="none", labCol="", cexRow = 1.2, margins = c(1,8),
          Colv=order(mirna.dds.vsd.idclc.val.dems.reptar$group), Rowv=TRUE,
          ColSideColors = c(plt2c[c(2,3)])[colData(mirna.dds.vsd.idclc.val.dems.reptar)$group],
          key =TRUE, key.title="Heatmap Key", lhei = c(1.1, 1),
          main="3 Validated DEMs between IDC and LC Targetting Repair Genes")
legend("topleft", legend = levels(mirna.dds.vsd.idclc.val.dems.reptar$group), col = plt2c[c(2,3)],
       lty= 1, lwd = 5, cex = 0.8, bty="0", bg="#FEFCF6",
       title="Samples", text.font=4, inset = c(0, 0.17), xjust = 0.5,
       text.width = 0.09)
dev.off()
################################################################################