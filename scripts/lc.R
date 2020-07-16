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
### LC
################################################################################
## Extracting dds.run.lc from dds.run
dds.run.lc <- dds.run[,dds.run$group %in% c("LC_Normal", "LC_Tumor")]
dds.run.lc
dds.run.lc$group <- factor(dds.run.lc$group)
table(dds.run.lc$group)
################################################################################
## Creating the results (RES) object from the original dds.run
res.lc <- results(dds.run, contrast=c("group", "LC_Tumor", "LC_Normal"), alpha=0.05, lfcThreshold=1)
summary(res.lc)
# out of 25175 gens: 162 upregulated and 119 downregulated

## Distribution of the LFC and p-adjusted value
summary(res.lc$log2FoldChange)
hist(res.lc$log2FoldChange, main="Distribution of the Log2 fold change in the results of LC", xlab="Log2 Fold Change")
summary(res.lc$padj)
hist(res.lc$padj, main="Distribution of the adjusted p-value in the results of LC", xlab="Adjusted p-value")

## complete cases only
res.lc.complete <- res.lc[complete.cases(res.lc),]
################################################################################
## Converting the results to DF
res.lc.df <- as.data.frame(res.lc)
res.lc.df <- res.lc.df[order(res.lc.df$padj),]
write.csv(res.lc.df, "saved_objects/LC.results.sorted.csv")
################################################################################
## Extracting DEGs
res.lc.df.degs <- res.lc.df[res.lc.df$padj<0.05 & abs(res.lc.df$log2FoldChange)>1,]
write.csv(res.lc.df.degs, "saved_objects/lc-degs-p.05-LFC1.csv")

res.lc.degs <- res.lc.complete[res.lc.complete$padj<0.05 & abs(res.lc.complete$log2FoldChange)>1,]
summary(res.lc.degs)

## Distribution of the p-adjusted value
summary(res.lc.degs$padj)
hist(res.lc.degs$padj, main="Distribution of the adjusted p-value in the DEGs of LC", xlab="Adjusted p-value")
################################################################################
## Repair genes subsetting 
res.lc.rep <- res.lc[rownames(res.lc) %in% repair.genes$symbol,]
summary(res.lc.rep)
## out of 285 genes, 4 up-regulated, 0 down-regulated
res.lc.degs.rep <- res.lc.degs[rownames(res.lc.degs) %in% repair.genes$symbol,]
summary(res.lc.degs.rep)

## Distribution of the LFC and p-adjusted value
summary(res.lc.rep$log2FoldChange)
hist(res.lc.rep$log2FoldChange, main="Distribution of the Log2 fold change in the results (Repair Genes Only) of LC", xlab="Log2 Fold Change")
summary(res.lc.rep$padj)
hist(res.lc.rep$padj, main="Distribution of the adjusted p-value in the results (Repair Genes Only) of LC", xlab="Adjusted p-value")

summary(res.lc.degs.rep$padj)
hist(res.lc.degs.rep$padj, main="Distribution of the adjusted p-value in the Repair DEGs of LC", xlab="Adjusted p-value")
################################################################################
## Plotting: 1. MAplot
plotMA(res.lc, alpha=0.05, main="MA plot of the not-mormalized DESeq Results of LC")
################################################################################
## Plotting: 2. EnhancedVolcano Plot
EnhancedVolcano(res.lc, 
                lab = rownames(res.lc), 
                x = 'log2FoldChange', 
                y = 'pvalue',
                title = "Volcano plot of the not-normalized DESeq results",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                legendPosition = "right")
################################################################################
## subset vsd
dds.vsd.lc <- dds.vsd[,dds.vsd$group %in% c("LC_Normal", "LC_Tumor")]
dds.vsd.lc
dds.vsd.lc$group <- factor(dds.vsd.lc$group)
table(dds.vsd.lc$group)

## DEGs
dds.vsd.lc.degs <- dds.vsd.lc[rownames(dds.vsd.lc) %in% rownames(res.lc.degs),]
dds.vsd.lc.degs

## Repair DEGs
dds.vsd.lc.degs.rep <- dds.vsd.lc.degs[rownames(dds.vsd.lc.degs) %in% repair.genes$symbol ,]
dds.vsd.lc.degs.rep
################################################################################
## Plotting: 3. PCA
# all genes
pca.plot.vsd.lc <- plotPCA(dds.vsd.lc, intgroup="group", ntop = 25531, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.lc, "percentVar"))
ggplot(pca.plot.vsd.lc, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vst-normalized DESeq DataSet for all genes (ntop = 25531)")
rm(percentVar)

# all degs
pca.plot.vsd.lc.degs <- plotPCA(dds.vsd.lc.degs, intgroup="group", ntop = 281, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.lc.degs, "percentVar"))
ggplot(pca.plot.vsd.lc.degs, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vst-normalized DESeq DataSet for all DEGs (ntop = 281)")
rm(percentVar)

# repair degs
pca.plot.vsd.lc.degs.rep <- plotPCA(dds.vsd.lc.degs.rep, intgroup="group", ntop = 4 , returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.lc.degs.rep, "percentVar"))
ggplot(pca.plot.vsd.lc.degs.rep, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vst-normalized DESeq DataSet for repair DEGs (ntop = 4)")
rm(percentVar)
################################################################################
## Plotting: 4. Dispersion estimate
plotDispEsts(dds.run.lc)
################################################################################
## Plotting: 7. Boxplot with dots
genes.list <- rownames(res.lc.degs.rep)
par(mfrow=c(2,2))
for (gen in genes.list){
  plotting_gene(dds.run.lc, gen)
}
par(mfrow=c(1,1))
################################################################################
## Plotting: 8. Heatmap & Dendogram
colors2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)
heatmap.2(assay(dds.vsd.lc.degs.rep), col=colors2,
          scale="row", trace="none", labCol=substring(colnames(dds.vsd.lc.degs.rep), 6),
          # dendrogram = "row",
          Colv=order(dds.vsd.lc.degs.rep$Sample.Type), Rowv=TRUE,
          ColSideColors = c(Primary.Tumor="darkgreen", Solid.Tissue.Normal="orange")[colData(dds.vsd.lc.degs.rep)$Sample.Type],
          key =TRUE, key.title="Heatmap Key",
          main="Heatmap of the top 4 differentially expressed repair genes in LC")
################################################################################