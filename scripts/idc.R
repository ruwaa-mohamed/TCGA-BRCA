### RNA - IDC
dds.run <- readRDS("saved_objects/dds.run.group.rds")
dds.vsd <- readRDS("saved_objects/dds.vsd.rds")
################################################################################
## Libraries
library(DESeq2)
library(EnhancedVolcano)
library(biomaRt)
library(dplyr)
library(RColorBrewer)
library(gplots)

source("scripts/proteins.R")
source("scripts/plotting_gene.R")
source("scripts/to_tiff.R")
################################################################################
## Extracting dds.run.idc from dds.run
dds.run.idc <- dds.run[,dds.run$group %in% c("IDC_Tumor", "IDC_Normal")]
dds.run.idc
# 25,531 gene, 856 samples

dds.run.idc$group <- droplevels(dds.run.idc$group)
table(dds.run.idc$group)
## 856 samples = 89 Normal + 767 Tumor
################################################################################
## Creating the results (RES) object from the original dds.run
res.idc <- results(dds.run, contrast=c("group", "IDC_Tumor", "IDC_Normal"), alpha=0.05, lfcThreshold=1)
summary(res.idc)
# out of 25,173 gens: 2,236 upregulated and 1,618 downregulated

## Distribution of the LFC and p-adjusted value
summary(res.idc$log2FoldChange)
to_tiff(
  hist(res.idc$log2FoldChange, main="Log2 Fold Change of the Genes between IDC and Normal", xlab="Log2 Fold Change"),
  "rna.idc.all.l2fc.hist.tiff")
summary(res.idc$padj)
to_tiff(
  hist(res.idc$padj, main="Adjusted p-value of the Genes between IDC and Normal", xlab="Adjusted p-value"),
  "rna.idc.all.adjpval.hist.tiff")
################################################################################
## complete cases only
res.idc.complete <- res.idc[complete.cases(res.idc),]
################################################################################
## Extracting DEGs
res.idc.degs <- res.idc.complete[res.idc.complete$padj<0.05 & abs(res.idc.complete$log2FoldChange)>1,]
summary(res.idc.degs)
# 3,854 DEGs: 2,236 up and 1,618 down

## Distribution of the p-adjusted value
summary(res.idc.degs$padj)
to_tiff(hist(
  res.idc.degs$padj, main="Adjusted p-value of the DEGs between IDC and Normal", xlab="Adjusted p-value"),
  "rna.idc.degs.adjpval.hist.tiff")
################################################################################
## Repair genes subsetting 
res.idc.rep <- res.idc.complete[rownames(res.idc.complete) %in% repair.genes$symbol,]
summary(res.idc.rep)
## out of 285 genes, 36 up-regulated, 0 down-regulated
res.idc.degs.rep <- res.idc.degs[rownames(res.idc.degs) %in% repair.genes$symbol,]
summary(res.idc.degs.rep)

## Distribution of the LFC and p-adjusted value
summary(res.idc.rep$log2FoldChange)
to_tiff(
  hist(res.idc.rep$log2FoldChange, main="Log2 Fold Change of All Repair Genes between IDC and Normal", xlab="Log2 Fold Change"),
  "rna.idc.rep.all.l2fc.hist.tiff")

summary(res.idc.rep$padj)
to_tiff(
  hist(res.idc.rep$padj, main="Adjusted p-value of All Repair Genes between IDC and Normal", xlab="Adjusted p-value"),
  "rna.idc.rep.all.adjpval.hist.tiff")

summary(res.idc.degs.rep$padj)
to_tiff(
  hist(res.idc.degs.rep$padj, main="Adjusted p-value of Repair DEGs between IDC and Normal", xlab="Adjusted p-value"),
  "rna.idc.rep.degs.adjpval.hist.tiff")
################################################################################
## Exporting to files
write.table(as.data.frame(res.idc.degs), "saved_objects/res.idc.degs.csv", sep=",", quote=FALSE)
write.table(as.data.frame(res.idc.degs.rep), "saved_objects/res.idc.degs.rep.csv", sep=",", quote=FALSE)
################################################################################
## Biomart part
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# filters <- listFilters(ensembl)
# attributes <- listAttributes(ensembl)
mart <- getBM(attributes = c("external_gene_name", "description"),
              filters = "external_gene_name",
              values = rownames(res.idc.degs), 
              mart = ensembl)

# ## Using the online Biomart
# # Ensembl Genes 101
# # Human Genes (GRCH38.p13)
# # Date: Aug. 30, 2020

# mart.online.rep <- read.csv("saved_objects/mart_export.rep.degs.txt", header=TRUE, na.strings = "")
# mart.online.rep <- mart.online.rep[complete.cases(mart.online.rep),]
# length(unique(mart.online.rep$Gene.name)) # 36 gene
# 
# mart.online <- read.csv("saved_objects/mart_export.degs.txt", header=TRUE, na.strings = "")
# mart.online <- mart.online[complete.cases(mart.online),]
# length(unique(mart.online$Gene.name)) # 3,116


## Adding mart to the DEGs
res.idc.degs.mart <- cbind(Gene.Name = rownames(res.idc.degs), as.data.frame(res.idc.degs))
res.idc.degs.mart <- left_join(res.idc.degs.mart, mart,
                               by=c("Gene.Name" = "external_gene_name"))
# res.idc.degs.mart <- left_join(res.idc.degs.mart, mart.online, 
#                                by = c("Gene.Name" = "Gene.name"))

## Adding mart to Repair DEGs
res.idc.degs.rep.mart <- cbind(Gene.Name = rownames(res.idc.degs.rep), as.data.frame(res.idc.degs.rep))
res.idc.degs.rep.mart <- left_join(res.idc.degs.rep.mart, mart, 
                                   by = c("Gene.Name" = "external_gene_name"))

## Rordering for Export
res.idc.degs.mart <- res.idc.degs.mart [, c("Gene.Name", "description", "baseMean", "log2FoldChange", "pvalue", "padj")]
res.idc.degs.mart <- res.idc.degs.mart %>% rename(Description = description)

res.idc.degs.rep.mart <- res.idc.degs.rep.mart [, c("Gene.Name", "description", "baseMean", "log2FoldChange", "pvalue", "padj")]
res.idc.degs.rep.mart <- res.idc.degs.rep.mart %>% rename(Description = description)

## Exporting the results
write.table(res.idc.degs.mart, "saved_objects/S1-res.idc.degs.csv", sep=",", quote=FALSE, row.names = FALSE)
write.table(res.idc.degs.rep.mart, "saved_objects/S4-res.idc.degs.rep.csv", sep=",", quote=FALSE, row.names = FALSE)
################################################################################
## Plotting: 1. MAplot
to_tiff(
  plotMA(res.idc.rep, alpha=0.05, main="MAplot of Repair Genes between IDC and Normal"),
  "rna.idc.rep.MAplot.tiff")
################################################################################
## Plotting: 2. EnhancedVolcano Plot
tiff("final_figures/rna.idc.rep.EnhancedVolcano.tiff", width = 18, height = 21, units = 'cm', res = 300)
EnhancedVolcano(res.idc.rep, 
                lab = rownames(res.idc.rep), 
                x = 'log2FoldChange', 
                y = 'padj',
                title = "Repair Genes between IDC and Normal",
                pCutoff = 0.05,
                FCcutoff = 1.0)
dev.off()
################################################################################
## subset vsd
dds.vsd.idc <- dds.vsd[,dds.vsd$group %in% c("IDC_Tumor", "IDC_Normal")]
dds.vsd.idc
dds.vsd.idc$group <- droplevels(dds.vsd.idc$group)
table(dds.vsd.idc$group)

## All Repair
dds.vsd.idc.rep <- dds.vsd.idc[rownames(dds.vsd.idc) %in% repair.genes$symbol,]
dds.vsd.idc.rep

## DEGs
dds.vsd.idc.degs <- dds.vsd.idc[rownames(dds.vsd.idc) %in% rownames(res.idc.degs),]
dds.vsd.idc.degs

## Repair DEGs
dds.vsd.idc.degs.rep <- dds.vsd.idc.degs[rownames(dds.vsd.idc.degs) %in% repair.genes$symbol ,]
dds.vsd.idc.degs.rep
################################################################################
## Plotting: 3. PCA
# All Genes
tiff("final_figures/rna.idc.pca.all.tiff", width = 18, height = 21, units = 'cm', res = 300)
pca.plot.vsd.idc <- plotPCA(dds.vsd.idc, intgroup="group", ntop = 25531, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.idc, "percentVar"))
ggplot(pca.plot.vsd.idc, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of All (25,531) Genes between IDC and Normal")
rm(pca.plot.vsd.idc, percentVar)
dev.off()

# All Repair 
tiff("final_figures/rna.idc.pca.rep.all.tiff", width = 18, height = 21, units = 'cm', res = 300)
pca.plot.vsd.idc.rep <- plotPCA(dds.vsd.idc.rep, intgroup="group", ntop = 285, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.idc.rep, "percentVar"))
ggplot(pca.plot.vsd.idc.rep, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of All Repair Genes (285) between IDC and Normal")
rm(pca.plot.vsd.idc.rep, percentVar)
dev.off()

# All DEGs
tiff("final_figures/rna.idc.pca.degs.tiff", width = 18, height = 21, units = 'cm', res = 300)
pca.plot.vsd.idc.degs <- plotPCA(dds.vsd.idc.degs, intgroup="group", ntop = 3854, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.idc.degs, "percentVar"))
ggplot(pca.plot.vsd.idc.degs, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of all DEGs (3,854) between IDC and Normal")
rm(pca.plot.vsd.idc.degs, percentVar)
dev.off()

# Repair DEGs
tiff("final_figures/rna.idc.pca.rep.degs.tiff", width = 18, height = 21, units = 'cm', res = 300)
pca.plot.vsd.idc.degs.rep <- plotPCA(dds.vsd.idc.degs.rep, intgroup="group", ntop = 36, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.idc.degs.rep, "percentVar"))
ggplot(pca.plot.vsd.idc.degs.rep, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of Repair DEGs (36) between IDC and Normal")
rm(pca.plot.vsd.idc.degs.rep, percentVar)
dev.off()
################################################################################
## Plotting: 4. Dispersion estimate
to_tiff(
  plotDispEsts(dds.run.idc, main="Disperssion Estimates for all Genes between IDC and Normal"),
  "rna.idc.plotDispEsts.tiff")
################################################################################
## Colors for Boxplot and Heatmap clusters
# plt2 <- brewer.pal(12,"Set3")[c(7, 6, 10)]
plt2 <- c("#B3DE69", "#FDB462", "#BC80BD")
plt2c <- c("#92CD2D", "#FC9C30", "#AA5EAC")
# library(scales)
# show_col(c(plt2, plt2c))
################################################################################
## Plotting: 7. Boxplot with dots
dds.vsd.bp <- dds.vsd[,dds.vsd$group %in% c("IDC_Tumor", "IDC_Normal", "LC_Tumor")]
dds.vsd.bp$group <- droplevels(dds.vsd.bp$group)
levels(dds.vsd.bp$group) <- c("Normal\nDuctal", "IDC", "LC")

genes.list <- rownames(res.idc.degs.rep)
source("scripts/plotting_gene.R")

for (gen in genes.list){
  tiff(paste0("final_figures/boxplots/rna.degs.rep.", gen, ".tiff"),
       width = 15, height = 21, units = 'cm', res = 300)
  plotting_gene(dds.vsd.bp, gen, l="", atr=1, colorslist = plt2c)
  dev.off()

  # to_tiff(
  #   plotting_gene(dds.vsd.bp, gen, l="", atr=1, colorslist = plt2c),
  #   paste0("boxplots/rna.degs.rep.", gen, ".tiff"))
}


pdf(file="final_figures/rna.vsd.rep.degs.boxplot.pdf", onefile=TRUE, paper="a4", width = 8, height = 11)
par(mfrow=c(3,3))
for (gen in genes.list){
  plotting_gene(dds.vsd.bp, gen, l="", atr=1, colorslist = plt2c)
}
dev.off()
par(mfrow=c(1,1))

# rm(gen, genes.list)


for (gen in c("BRCA1", "TP53")){
  tiff(paste0("final_figures/boxplots/rna.nondegs.rep.", gen, ".tiff"),
       width = 15, height = 21, units = 'cm', res = 300)
  plotting_gene(dds.vsd.bp, gen, l="", atr=1, colorslist = plt2c)
  dev.off()
  # to_tiff(
  #   plotting_gene(dds.vsd.bp, gen, l="", atr=1, colorslist = plt2c),
  #   paste0("boxplots/rna.rep.nondeg.", gen, ".tiff"))
}


pdf(file="final_figures/rna.vsd.rep.nondegs.boxplot.pdf", onefile=TRUE, paper="a4", width = 8, height = 11)
par(mfrow=c(3,3))
for (gen in c("BRCA1", "TP53")){
  plotting_gene(dds.vsd.bp, gen, l="", atr=1, colorslist = plt2c)
}
dev.off()
par(mfrow=c(1,1))

################################################################################
## Plotting: 8. Heatmap & Dendogram
colors2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)
levels(dds.vsd.idc.degs.rep$group) <- c("Normal\nDuctal", "IDC")

# pdf(file="final_figures/rna.idc.heatmap.rep.degs.pdf", onefile=TRUE, paper="a4r", width = 11, height = 8)
tiff("final_figures/rna.idc.heatmap.rep.degs.tiff", width = 25, height = 18, units = 'cm', res = 300)
heatmap.2(assay(dds.vsd.idc.degs.rep), col=colors2,
          scale="row", trace="none", labCol="", #dendrogram = "row",
          Colv=order(dds.vsd.idc.degs.rep$Sample.Type), Rowv=TRUE,
          ColSideColors = plt2c[c(1,2)][colData(dds.vsd.idc.degs.rep)$group],
          key =TRUE, key.title="Heatmap Key",
          main="36 Repair DEGs between IDC and Normal")
legend("topleft", legend = levels(dds.vsd.idc.degs.rep$group), col = plt2c[c(1,2)],
       lty= 1, lwd = 5, cex=0.75, bty="0", bg="#FEFCF6",
       title="Samples", text.font=4, inset = c(0, 0.17), xjust = 0.5,
       text.width = 0.09)
dev.off()
################################################################################