################################################################################
### miRNA IDC_Tumor vs. IDC_Normal

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
mirna.dds.run.idc <- mirna.dds.run[,mirna.dds.run$group %in% c("IDC_Normal", "IDC_Tumor")]
mirna.dds.run.idc

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

## Distribution of the p-adjusted value
summary(mirna.res.idc.dems$padj)
to_tiff(
  hist(mirna.res.idc.dems$padj, main="Adjusted p-value in DEMs between IDC and Normal", xlab="Adjusted p-value"),
  "mirna.idc.dems.adjpval.tiff")

################################################################################
## Exporting to file
write.table(as.data.frame(mirna.res.idc.dems), "saved_objects/mirna.res.idc.dems.csv", sep=",", quote=FALSE)
################################################################################
## Plotting: 1. MAplot
to_tiff(plotMA(mirna.res.idc, alpha=0.05, main="MAplot of the miRs between IDC and Normal"), 
        "mirna.idc.all.plotMA.tiff")
################################################################################
## Plotting: 2. EnhancedVolcano Plot
tiff("final_figures/mirna.idc.all.volcano.tiff", width = 18, height = 21, units = 'cm', res = 300)
EnhancedVolcano(mirna.res.idc, 
                lab = rownames(mirna.res.idc), 
                x = 'log2FoldChange', 
                y = 'pvalue',
                title = "Volcano plot of the miRs between IDC and Normal",
                pCutoff = 0.05,
                FCcutoff = 1.0)
dev.off()
################################################################################
## Plotting: 3. PCA

## subset vsd
mirna.dds.vsd.idc <- mirna.dds.vsd[,mirna.dds.vsd$group %in% c("IDC_Tumor", "IDC_Normal")]
mirna.dds.vsd.idc

## DEMs
mirna.dds.vsd.idc.dems <- mirna.dds.vsd.idc[rownames(mirna.dds.vsd.idc) %in% rownames(mirna.res.idc.dems),]
mirna.dds.vsd.idc.dems

# PCA all miRs
tiff("final_figures/mirna.idc.pca.all.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca.idc <- plotPCA(mirna.dds.vsd.idc, intgroup="group", ntop = 1881, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca.idc, "percentVar"))
ggplot(mirna.pca.idc, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of all miRs (1,881) between IDC and Normal")
rm(mirna.pca.idc, percentVar)
dev.off()

# PCA all dems
tiff("final_figures/mirna.idc.pca.dems.tiff", width = 18, height = 21, units = 'cm', res = 300)
mirna.pca.idc.dems <- plotPCA(mirna.dds.vsd.idc.dems, intgroup="group", ntop = 171, returnData=TRUE)
percentVar <- round(100 * attr(mirna.pca.idc.dems, "percentVar"))
ggplot(mirna.pca.idc.dems, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of all DEMs (171) between IDC and Normal")
rm(mirna.pca.idc.dems, percentVar)
dev.off()

################################################################################
## Plotting: 8. Heatmap & Dendogram
pdf(file="final_figures/mirna.idc.heatmap.dems.pdf", onefile=TRUE, paper="a4r", width = 11, height = 8)
colors2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)
heatmap.2(assay(mirna.dds.vsd.idc.dems), col=colors2,
          scale="row", trace="none", labCol=substring(colnames(mirna.dds.vsd.idc.dems), 6),
          Colv=order(mirna.dds.vsd.idc.dems$Sample.Type), Rowv=TRUE,
          ColSideColors = c(Normal="darkgreen", Tumor="orange")[colData(mirna.dds.vsd.idc.dems)$Sample.Type],
          key =TRUE, key.title="Heatmap Key",
          main="Heatmap of the 171 DEMs between IDC (orange) and Normal (darkgreen)")
dev.off()
################################################################################
## Plotting: 7. Boxplot with dots
dem.list <- rownames(mirna.res.idc.dems)

pdf(file="final_figures/mirna.idc.dems.boxplot.pdf", onefile=TRUE, paper="a4", width = 8, height = 11)
par(mfrow=c(3,3))
for (dem in dem.list){
  plotting_gene(mirna.dds.run.idc, dem, "")
}
dev.off()

par(mfrow=c(1,1))
rm(dem, dem.list)
################################################################################
################################################################################
## mirTarBase (experimentally validated only)
## subsetting the DB by the DEMs
hsa_MTI.idc.dems <- hsa_MTI[hsa_MTI$miRNA %in% rownames(mirna.res.idc.dems),]

## Exploring
length(unique(hsa_MTI.idc.dems$miRNA)) 
# 32 mir out of 171
length(unique(hsa_MTI.idc.dems$Target.Gene)) 
# targeting total of 3,590 genes
sum(unique(hsa_MTI.idc.dems$Target.Gene) %in% repair.genes$symbol) 
# 79 of the 3,590 genes are repair
################################################################################