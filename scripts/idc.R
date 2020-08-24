################################################################################
################################################################################
## R version 4.0.1 (2020-06-06) -- "See Things Now"
## Copyright (C) 2020 The R Foundation for Statistical Computing
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04 LTS

## BiocManager version: 1.30.10
## Bioconductor version: 3.11
## Date: Wednesday, July 15, 2020.
## Created by Ruwaa I. Mohamed
################################################################################
################################################################################
### IDC
################################################################################
## Extracting dds.run.idc from dds.run
dds.run.idc <- dds.run[,dds.run$group %in% c("IDC_Normal", "IDC_Tumor")]
dds.run.idc
dds.run.idc$group <- factor(dds.run.idc$group)
dds.run.idc$group <- droplevels(dds.run.idc$group)
table(dds.run.idc$group)
################################################################################
## Creating the results (RES) object from the original dds.run
res.idc <- results(dds.run, contrast=c("group", "IDC_Tumor", "IDC_Normal"), alpha=0.05, lfcThreshold=1)
summary(res.idc)
# out of 25173 gens: 2,236 upregulated and 1,618 downregulated

## Distribution of the LFC and p-adjusted value
summary(res.idc$log2FoldChange)
to_tiff(
  hist(res.idc$log2FoldChange, main="Log2 Fold Change of the Genes between IDC and Normal", xlab="Log2 Fold Change"),
  "rna.idc.all.l2fc.hist.tiff")
summary(res.idc$padj)
to_tiff(
  hist(res.idc$padj, main="Adjusted p-value of the Genes between IDC and Normal", xlab="Adjusted p-value"),
  "rna.idc.all.adjpval.hist.tiff")

## complete cases only
res.idc.complete <- res.idc[complete.cases(res.idc),]
################################################################################
# ## Converting the results to DF
# res.idc.df <- as.data.frame(res.idc)
# res.idc.df <- res.idc.df[order(res.idc.df$padj),]
# write.csv(res.idc.df, "saved_objects/IDC.results.sorted.csv")
################################################################################
## Extracting DEGs
res.idc.degs <- res.idc.complete[res.idc.complete$padj<0.05 & abs(res.idc.complete$log2FoldChange)>1,]
summary(res.idc.degs)

# write.table(as.data.frame(res.idc.degs), "saved_objects/res.idc.degs.csv", sep=",", quote=FALSE)

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
## Plotting: 1. MAplot
to_tiff(
  plotMA(res.idc, alpha=0.05, main="MAplot of All Genes between IDC and Normal"),
  "rna.idc.all.MAplot.tiff")
################################################################################
## Plotting: 2. EnhancedVolcano Plot
tiff("final_figures/rna.idc.EnhancedVolcano.tiff", width = 18, height = 21, units = 'cm', res = 300)
EnhancedVolcano(res.idc, 
                lab = rownames(res.idc), 
                x = 'log2FoldChange', 
                y = 'pvalue',
                title = "All Genes between IDC and Normal",
                pCutoff = 0.05,
                FCcutoff = 1.0)
dev.off()
################################################################################
## subset vsd
dds.vsd.idc <- dds.vsd[,dds.vsd$group %in% c("IDC_Normal", "IDC_Tumor")]
dds.vsd.idc
dds.vsd.idc$group <- factor(dds.vsd.idc$group)
dds.vsd.idc$group <- droplevels(dds.vsd.idc$group)
table(dds.vsd.idc$group)

## DEGs
dds.vsd.idc.degs <- dds.vsd.idc[rownames(dds.vsd.idc) %in% rownames(res.idc.degs),]
dds.vsd.idc.degs

## Repair DEGs
dds.vsd.idc.degs.rep <- dds.vsd.idc.degs[rownames(dds.vsd.idc.degs) %in% repair.genes$symbol ,]
dds.vsd.idc.degs.rep
################################################################################
## Plotting: 3. PCA
# all genes
tiff("final_figures/rna.idc.pca.all.tiff", width = 18, height = 21, units = 'cm', res = 300)
pca.plot.vsd.idc <- plotPCA(dds.vsd.idc, intgroup="group", ntop = 25531, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.idc, "percentVar"))
ggplot(pca.plot.vsd.idc, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of All (25,531) Genes between IDC and Normal")
rm(percentVar)
dev.off()

# all degs
tiff("final_figures/rna.idc.pca.degs.tiff", width = 18, height = 21, units = 'cm', res = 300)
pca.plot.vsd.idc.degs <- plotPCA(dds.vsd.idc.degs, intgroup="group", ntop = 3854, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.idc.degs, "percentVar"))
ggplot(pca.plot.vsd.idc.degs, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of all DEGs (3,854) between IDC and Normal")
rm(percentVar)
dev.off()

# repair degs
tiff("final_figures/rna.idc.pca.rep.degs.tiff", width = 18, height = 21, units = 'cm', res = 300)
pca.plot.vsd.idc.degs.rep <- plotPCA(dds.vsd.idc.degs.rep, intgroup="group", ntop = 36, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.idc.degs.rep, "percentVar"))
ggplot(pca.plot.vsd.idc.degs.rep, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of Repair DEGs (36) between IDC and Normal")
rm(percentVar)
dev.off()
################################################################################
## Plotting: 4. Dispersion estimate
to_tiff(
  plotDispEsts(dds.run.idc, main="Disperssion Estimates for all Genes between IDC and Normal"),
  "rna.idc.plotDispEsts.tiff")
################################################################################
## Plotting: 7. Boxplot with dots
genes.list <- rownames(res.idc.degs.rep)
# "H4C1", "NEIL3"

pdf(file="final_figures/rna.idc.rep.degs.boxplot.pdf", onefile=TRUE, paper="a4", width = 8, height = 11)
par(mfrow=c(3,3))
for (gen in genes.list){
  if (gen %in% c("H4C1", "NEIL3")){
    plotting_gene(dds.run.idc, gen, "")
  } else {
    plotting_gene(dds.run.idc, gen)
  }
}
dev.off()

par(mfrow=c(1,1))
################################################################################
## Plotting: 8. Heatmap & Dendogram
pdf(file="final_figures/rna.idc.heatmap.rep.degs.pdf", onefile=TRUE, paper="a4r", width = 11, height = 8)

colors2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)
heatmap.2(assay(dds.vsd.idc.degs.rep), col=colors2,
          scale="row", trace="none", labCol=substring(colnames(dds.vsd.idc.degs.rep), 6),
          # dendrogram = "row",
          Colv=order(dds.vsd.idc.degs.rep$Sample.Type), Rowv=TRUE,
          ColSideColors = c(Normal="darkgreen", Tumor="orange")[colData(dds.vsd.idc.degs.rep)$Sample.Type],
          key =TRUE, key.title="Heatmap Key",
          main="Heatmap of the 36 Repair DEGS between IDC (orange) and Normal (darkgreen)")
dev.off()
################################################################################