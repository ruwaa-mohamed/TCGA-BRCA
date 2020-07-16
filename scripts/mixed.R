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
### Mixed
################################################################################
## Extracting dds.run.mixed from dds.run
dds.run.mixed <- dds.run[,dds.run$group %in% c("Mixed_Normal", "Mixed_Tumor")]
dds.run.mixed
dds.run.mixed$group <- factor(dds.run.mixed$group)
table(dds.run.mixed$group)
################################################################################
## Creating the results (RES) object from the original dds.run
res.mixed <- results(dds.run, contrast=c("group", "Mixed_Tumor", "Mixed_Normal"), alpha=0.05,lfcThreshold=1)
summary(res.mixed)
# out of 25175 gens: 387 upregulated and 620 downregulated

## Distribution of the LFC and p-adjusted value
summary(res.mixed$log2FoldChange)
hist(res.mixed$log2FoldChange, main="Distribution of the Log2 fold change in the results of Mixed", xlab="Log2 Fold Change")
summary(res.mixed$padj)
hist(res.mixed$padj, main="Distribution of the adjusted p-value in the results of Mixed", xlab="Adjusted p-value")

## complete cases only
res.mixed.complete <- res.mixed[complete.cases(res.mixed),]
################################################################################
## Converting the results to DF
res.mixed.df <- as.data.frame(res.mixed)
res.mixed.df <- res.mixed.df[order(res.mixed.df$padj),]
write.csv(res.mixed.df, "saved_objects/Mixed.results.sorted.csv")
################################################################################
## Extracting DEGs
res.mixed.df.degs <- res.mixed.df[res.mixed.df$padj<0.05 & abs(res.mixed.df$log2FoldChange)>1,]
write.csv(res.mixed.df.degs, "saved_objects/mixed-degs-p.05-LFC1.csv")

res.mixed.degs <- res.mixed.complete[res.mixed.complete$padj<0.05 & abs(res.mixed.complete$log2FoldChange)>1,]
summary(res.mixed.degs)

## Distribution of the p-adjusted value
summary(res.mixed.degs$padj)
hist(res.mixed.degs$padj, main="Distribution of the adjusted p-value in the DEGs of Mixed", xlab="Adjusted p-value")
################################################################################
## Repair genes subsetting 
res.mixed.rep <- res.mixed[rownames(res.mixed) %in% repair.genes$symbol,]
summary(res.mixed.rep)
## out of 285 genes, 6 up-regulated, 0 down-regulated
res.mixed.degs.rep <- res.mixed.degs[rownames(res.mixed.degs) %in% repair.genes$symbol,]
summary(res.mixed.degs.rep)

## Distribution of the LFC and p-adjusted value
summary(res.mixed.rep$log2FoldChange)
hist(res.mixed.rep$log2FoldChange, main="Distribution of the Log2 fold change in the results (Repair Genes Only) of Mixed", xlab="Log2 Fold Change")
summary(res.mixed.rep$padj)
hist(res.mixed.rep$padj, main="Distribution of the adjusted p-value in the results (Repair Genes Only) of Mixed", xlab="Adjusted p-value")

summary(res.mixed.degs.rep$padj)
hist(res.mixed.degs.rep$padj, main="Distribution of the adjusted p-value in the Repair DEGs of Mixed", xlab="Adjusted p-value")
################################################################################
## Plotting: 1. MAplot
plotMA(res.mixed, alpha=0.05, main="MA plot of the not-mormalized DESeq Results of Mixed")
################################################################################
## Plotting: 2. EnhancedVomixedano Plot
EnhancedVomixedano(res.mixed, 
                lab = rownames(res.mixed), 
                x = 'log2FoldChange', 
                y = 'pvalue',
                title = "Vomixedano plot of the not-normalized DESeq results",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                legendPosition = "right")
################################################################################
## subset vsd
dds.vsd.mixed <- dds.vsd[,dds.vsd$group %in% c("Mixed_Normal", "Mixed_Tumor")]
dds.vsd.mixed
dds.vsd.mixed$group <- factor(dds.vsd.mixed$group)
table(dds.vsd.mixed$group)

## DEGs
dds.vsd.mixed.degs <- dds.vsd.mixed[rownames(dds.vsd.mixed) %in% rownames(res.mixed.degs),]
dds.vsd.mixed.degs

## Repair DEGs
dds.vsd.mixed.degs.rep <- dds.vsd.mixed.degs[rownames(dds.vsd.mixed.degs) %in% repair.genes$symbol ,]
dds.vsd.mixed.degs.rep
################################################################################
## Plotting: 3. PCA
# all genes
pca.plot.vsd.mixed <- plotPCA(dds.vsd.mixed, intgroup="group", ntop = 25531, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.mixed, "percentVar"))
ggplot(pca.plot.vsd.mixed, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vst-normalized DESeq DataSet for all genes (ntop = 25531)")
rm(percentVar)

# all degs
pca.plot.vsd.mixed.degs <- plotPCA(dds.vsd.mixed.degs, intgroup="group", ntop = 1007, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.mixed.degs, "percentVar"))
ggplot(pca.plot.vsd.mixed.degs, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vst-normalized DESeq DataSet for all DEGs (ntop = 1007)")
rm(percentVar)

# repair degs
pca.plot.vsd.mixed.degs.rep <- plotPCA(dds.vsd.mixed.degs.rep, intgroup="group", ntop = 6 , returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.mixed.degs.rep, "percentVar"))
ggplot(pca.plot.vsd.mixed.degs.rep, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vst-normalized DESeq DataSet for repair DEGs (ntop = 4)")
rm(percentVar)
################################################################################
## Plotting: 4. Dispersion estimate
plotDispEsts(dds.run.mixed)
################################################################################
## Plotting: 7. Boxplot with dots
genes.list <- rownames(res.mixed.degs.rep)
par(mfrow=c(2,3))
for (gen in genes.list){
  plotting_gene(dds.run.mixed, gen)
}
par(mfrow=c(1,1))
################################################################################
## Plotting: 8. Heatmap & Dendogram
colors2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)
heatmap.2(assay(dds.vsd.mixed.degs.rep), col=colors2,
          scale="row", trace="none", labCol=substring(colnames(dds.vsd.mixed.degs.rep), 6),
          # dendrogram = "row",
          Colv=order(dds.vsd.mixed.degs.rep$Sample.Type), Rowv=TRUE,
          ColSideColors = c(Primary.Tumor="darkgreen", Solid.Tissue.Normal="orange")[colData(dds.vsd.mixed.degs.rep)$Sample.Type],
          key =TRUE, key.title="Heatmap Key",
          main="Heatmap of the top 6 differentially expressed repair genes in Mixed")
################################################################################