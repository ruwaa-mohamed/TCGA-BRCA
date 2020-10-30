################################################################################
### miRNA IDC_Tumor vs. IDC_Normal

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

mir.gene.pairs <- read.csv("saved_objects/mir.gene.pairs.csv", header=TRUE)
mir.gene.pairs.rep <- read.csv("saved_objects/mir.gene.pairs.rep.csv", header=TRUE)
mir.gene.pairs.merged <- read.csv("saved_objects/mir.gene.pairs.merged.csv", header=TRUE)
################################################################################
mirna.dds.run.idc <- mirna.dds.run[,mirna.dds.run$group %in% c("IDC_Tumor", "IDC_Normal")]
mirna.dds.run.idc
# 1881 miRs and 838 samples
mirna.dds.run.idc$group <- droplevels(mirna.dds.run.idc$group)
levels(mirna.dds.run.idc$group) <- c("Normal\nDuctal", "IDC")
table(mirna.dds.run.idc$group)

## Creating the mirna.res object from the original mirna.dds
mirna.res.idc <- results(mirna.dds.run, contrast=c("group", "IDC_Tumor", "IDC_Normal"), alpha=0.05, lfcThreshold=1)
summary(mirna.res.idc)
## out ot 1,598 mirs, 109 upregulated and 62 downregulated

## Distribution of the LFC and p-adjusted value
summary(mirna.res.idc$log2FoldChange)
to_tiff(
  hist(mirna.res.idc$log2FoldChange, main="Log2 Fold Change in miRs between IDC and Normal", xlab="Log2 Fold Change"),
  "mirna.idc.all.l2fc.tiff")
summary(mirna.res.idc$padj)
to_tiff(
  hist(mirna.res.idc$padj, main="Adjusted p-value in miRs between IDC and Normal", xlab="Adjusted p-value"),
  "mirna.idc.all.adjpval.tiff")
################################################################################
## complete cases only
mirna.res.idc.complete <- mirna.res.idc[complete.cases(mirna.res.idc),]
################################################################################
## Extracting DEMs
mirna.res.idc.dems <- mirna.res.idc.complete[mirna.res.idc.complete$padj<0.05 & abs(mirna.res.idc.complete$log2FoldChange)>1,]
summary(mirna.res.idc.dems)
# 171 DEMs: 109 up and 62 down

## Distribution of the p-adjusted value
summary(mirna.res.idc.dems$padj)
to_tiff(
  hist(mirna.res.idc.dems$padj, main="Adjusted p-value in DEMs between IDC and Normal", xlab="Adjusted p-value"),
  "mirna.idc.dems.adjpval.tiff")
################################################################################
## Exporting to file
write.table(as.data.frame(mirna.res.idc.dems), "saved_objects/mirna.res.idc.dems.all.csv", sep=",", quote=FALSE)
################################################################################
## Validated Only from mirTarBase DB
## All miRs
mirna.res.idc.val <- mirna.res.idc.complete[rownames(mirna.res.idc.complete) %in% mir.gene.pairs$miRNA,]
summary(mirna.res.idc.val)
# only 159 out of 764 complete miR (mirna.res.idc.complete) 

## All DEMs
mirna.res.idc.dems.val <- mirna.res.idc.dems[rownames(mirna.res.idc.dems) %in% mir.gene.pairs$miRNA,]
summary(mirna.res.idc.dems.val)
# only 32 out of 171 DEMs 

## miRs affecting any repair
mirna.res.idc.val.reptar <- mirna.res.idc.val[rownames(mirna.res.idc.val) %in% mir.gene.pairs.rep$miRNA,]
summary(mirna.res.idc.val.reptar)
# 116 miRs (out of 159 validated miR): 15 up and 7 down

## DEMs affecting any repair
mirna.res.idc.dems.val.reptar <- mirna.res.idc.dems[rownames(mirna.res.idc.dems) %in% mir.gene.pairs.rep$miRNA,]
summary(mirna.res.idc.dems.val.reptar)
# 22 DEMs (out of 32 validated DEMs): 15 up and 7 down
################################################################################
## Adding Target Genes part (for Export)
## All DEMs with Target Genes
mirna.res.idc.dems.val.df <- cbind(miRs = row.names(mirna.res.idc.dems.val), as.data.frame(mirna.res.idc.dems.val))
mirna.res.idc.dems.val.df <- left_join(mirna.res.idc.dems.val.df, mir.gene.pairs.merged, by=c("miRs" = "miRNA"))
mirna.res.idc.dems.val.df <- mirna.res.idc.dems.val.df [, c("miRs", "baseMean", "log2FoldChange", "pvalue", "padj", "Target.Gene")]

## Exporting to file
write.table(mirna.res.idc.dems.val.df, "saved_objects/S5-mirna.res.idc.dems.val.csv", sep=",", quote=FALSE, row.names = FALSE)
################################################################################
################################################################################
## Subsetting VSD
## All Validated miRs in IDC (Tumor and Normal)
mirna.dds.vsd.idc.val <- mirna.dds.vsd[rownames(mirna.dds.vsd) %in% rownames(mirna.res.idc.val), mirna.dds.vsd$group %in% c("IDC_Tumor", "IDC_Normal")]
mirna.dds.vsd.idc.val
mirna.dds.vsd.idc.val$group <- droplevels(mirna.dds.vsd.idc.val$group)
levels(mirna.dds.vsd.idc.val$group) <- c("Ductal\nNormal", "IDC")
table(mirna.dds.vsd.idc.val$group)

## Validated miRs Targetting any Repair Genes
mirna.dds.vsd.idc.val.reptar <- mirna.dds.vsd.idc.val[rownames(mirna.dds.vsd.idc.val) %in% mir.gene.pairs.rep$miRNA,]
mirna.dds.vsd.idc.val.reptar

## All Validated DEMs
mirna.dds.vsd.idc.val.dems <- mirna.dds.vsd.idc.val[rownames(mirna.dds.vsd.idc.val) %in% rownames(mirna.res.idc.dems.val),]
mirna.dds.vsd.idc.val.dems

## Validated DEMs Targetting any Repair Genes
mirna.dds.vsd.idc.val.dems.reptar <- mirna.dds.vsd.idc.val.dems[rownames(mirna.dds.vsd.idc.val.dems) %in% mir.gene.pairs.rep$miRNA,]
mirna.dds.vsd.idc.val.dems.reptar
################################################################################
## Plotting: 1. MAplot
to_tiff(plotMA(mirna.res.idc.val, alpha=0.05, main="MAplot of All miRs between IDC and Ductal Normal"), 
        "mirna/mirna.idc.val.MAplot.tiff")
################################################################################
## Plotting: 2. EnhancedVolcano Plot
tiff("final_figures/mirna/mirna.idc.val.all.volcano.tiff", width = 18, height = 21, units = 'cm', res = 300)
EnhancedVolcano(mirna.res.idc.val, 
                lab = rownames(mirna.res.idc.val), 
                x = 'log2FoldChange', 
                y = 'padj',
                title = "All miRs between IDC and Ductal Normal",
                pCutoff = 0.05,
                FCcutoff = 1.0)
dev.off()

tiff("final_figures/mirna/mirna.idc.val.rep.volcano.tiff", width = 18, height = 21, units = 'cm', res = 300)
EnhancedVolcano(mirna.res.idc.val.reptar, 
                lab = rownames(mirna.res.idc.val.reptar), 
                x = 'log2FoldChange', 
                y = 'padj',
                title = "Repair miRs between IDC and Ductal Normal",
                pCutoff = 0.05,
                FCcutoff = 1.0)
dev.off()
################################################################################
## Plotting: 3. PCA

## All miRs with Validated Target Genes
tiff("final_figures/mirna/mirna.idc.val.mirs.all.pca.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca.idc.val <- plotPCA(mirna.dds.vsd.idc.val, intgroup="group", ntop = 159, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca.idc.val, "percentVar"))
ggplot(mirna.pca.idc.val, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("All Validated miRs (159) between IDC and Ductal Normal")
rm(mirna.pca.idc.val, percentVar)
dev.off()

## All DEMs with Validated Target Genes 
tiff("final_figures/mirna/mirna.idc.val.dems.all.pca.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca.idc.val.dems <- plotPCA(mirna.dds.vsd.idc.val.dems, intgroup="group", ntop = 32, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca.idc.val.dems, "percentVar"))
ggplot(mirna.pca.idc.val.dems, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("All Validated DEMs (32) between IDC and Ductal Normal")
rm(mirna.pca.idc.val.dems, percentVar)
dev.off()

## All miRs with Validated Target Repair Genes
tiff("final_figures/mirna/mirna.idc.val.mirs.reptar.pca.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca.idc.val.reptar <- plotPCA(mirna.dds.vsd.idc.val.reptar, intgroup="group", ntop = 116, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca.idc.val.reptar, "percentVar"))
ggplot(mirna.pca.idc.val.reptar, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Validated Repair miRs (116) between IDC and Ductal Normal")
rm(mirna.pca.idc.val.reptar, percentVar)
dev.off()

## All DEMs with Validated Target Repair Genes 
tiff("final_figures/mirna/mirna.idc.val.dems.reptar.pca.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca.idc.val.dems.reptar <- plotPCA(mirna.dds.vsd.idc.val.dems.reptar, intgroup="group", ntop = 22, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca.idc.val.dems.reptar, "percentVar"))
ggplot(mirna.pca.idc.val.dems.reptar, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Validated Repair DEMs (22) between IDC and Ductal Normal")
rm(mirna.pca.idc.val.dems.reptar, percentVar)
dev.off()
################################################################################
## Plotting: 8. Heatmap & Dendogram
colors2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)
plt2c <- c("#92CD2D", "#FC9C30", "#AA5EAC")

## All DEMs with Validated Target Genes 
tiff("final_figures/mirna/mirna.idc.val.dems.all.heatmap.tiff", width = 25, height = 18, units = 'cm', res = 300)
heatmap.2(assay(mirna.dds.vsd.idc.val.dems), col=colors2,
          scale="row", trace="none", labCol="", cexRow = 1.2, margins = c(1,8),
          Colv=order(mirna.dds.vsd.idc.val.dems$group), Rowv=TRUE,
          ColSideColors = c(plt2c[c(1,2)])[colData(mirna.dds.vsd.idc.val.dems)$group],
          key =TRUE, key.title="Heatmap Key",
          main="32 Validated DEMs between IDC and Ductal Normal")
legend("topleft", legend = levels(mirna.dds.vsd.idc.val.dems$group), col = plt2c[c(1,2)],
       lty= 1, lwd = 5, cex=0.75, bty="0", bg="#FEFCF6",
       title="Samples", text.font=4, inset = c(0, 0.17), xjust = 0.5,
       text.width = 0.09)
dev.off()

## Validated DEMs with Target Repair Genes 
tiff("final_figures/mirna/mirna.idc.val.dems.reptar.heatmap.tiff", width = 25, height = 18, units = 'cm', res = 300)
heatmap.2(assay(mirna.dds.vsd.idc.val.dems.reptar), col=colors2,
          scale="row", trace="none", labCol="", cexRow = 1.2, margins = c(1,8),
          Colv=order(mirna.dds.vsd.idc.val.dems.reptar$group), Rowv=TRUE,
          ColSideColors = c(plt2c[c(1,2)])[colData(mirna.dds.vsd.idc.val.dems.reptar)$group],
          key =TRUE, key.title="Heatmap Key",
          main="22 Repair DEMs between IDC and Ductal Normal")
legend("topleft", legend = levels(mirna.dds.vsd.idc.val.dems.reptar$group), col = plt2c[c(1,2)],
       lty= 1, lwd = 5, cex=0.75, bty="0", bg="#FEFCF6",
       title="Samples", text.font=4, inset = c(0, 0.17), xjust = 0.5,
       text.width = 0.09)
dev.off()
################################################################################
## Plotting: 7. Boxplot with dots
dem.list <- rownames(mirna.res.idc.dems.val.reptar)
dem.list <- unique(c(rownames(mirna.res.idc.dems.val), rownames(mirna.res.idclc.dems.val)))

mirna.dds.vsd.bp <- mirna.dds.vsd[rownames(mirna.dds.vsd) %in% rownames(mirna.dds.vsd.idc.val), mirna.dds.vsd$group %in% c("IDC_Normal", "IDC_Tumor", "LC_Tumor")]
mirna.dds.vsd.bp$group <- droplevels(mirna.dds.vsd.bp$group)
levels(mirna.dds.vsd.bp$group) <- c("Ductal\nNormal", "IDC", "LC")
table(mirna.dds.vsd.bp$group)

source("scripts/plotting_gene_2.R")

## PDF
pdf(file="final_figures/mirna/boxplots/notreptar/mirna.idc.idclc.dems.boxplot.ylim.pdf", onefile=TRUE, paper="a4", width = 8, height = 11)
par(mfrow=c(3,3))
for (dem in dem.list){
  plotting_gene_2(mirna.dds.vsd.bp, dem, l="", atr=1, colorslist = plt2c)
}
dev.off()
par(mfrow=c(1,1))

## TIFF
for (dem in dem.list){
  tiff(paste0("final_figures/mirna/boxplots/notreptar/mirna.dems.", dem, 'ylim.tiff'),
       width = 15, height = 21, units = 'cm', res = 300)
  plotting_gene_2(mirna.dds.vsd.bp, dem, l="", atr=1, colorslist = plt2c)
  dev.off()
  
  # to_tiff(
  #   plotting_gene_2(mirna.dds.vsd.bp, dem, l="", atr=1, colorslist = plt2c),
  #   paste0("mirna/boxplots/mirna.dems.reptar.", dem, '.tiff'))
}

rm(dem, dem.list)
################################################################################