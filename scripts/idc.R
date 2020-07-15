### IDC
################################################################################
## Extracting dds.run.idc from dds.run
dds.run.idc <- dds.run[,dds.run$group %in% c("IDC_Normal", "IDC_Tumor")]
dds.run.idc
table(dds.run.idc$group)
dds.run.idc$group <- factor(dds.run.idc$group)
table(dds.run.idc$group)
################################################################################
## Creating the results (RES) object from the original dds.run
res.idc <- results(dds.run, contrast=c("group", "IDC_Tumor", "IDC_Normal"), alpha=0.05,lfcThreshold=1)
summary(res.idc)

## Distribution of the LFC and p-adjusted value
summary(res.idc$log2FoldChange)
hist(res.idc$log2FoldChange, main="Distribution of the Log2 fold change in the results of IDC", xlab="Log2 Fold Change")
summary(res.idc$padj)
hist(res.idc$padj, main="Distribution of the adjusted p-value in the results of IDC", xlab="Adjusted p-value")

## complete cases only
res.idc.complete <- res.idc[complete.cases(res.idc),]
################################################################################
## Converting the results to DF
res.idc.df <- as.data.frame(res.idc)
res.idc.df <- res.idc.df[order(res.idc.df$padj),]
write.csv(res.idc.df, "saved_objects/IDC.results.sorted.csv")
################################################################################
## Extracting DEGs
res.idc.df.degs <- res.idc.df[res.idc.df$padj<0.05 & abs(res.idc.df$log2FoldChange)>1,]
write.csv(res.idc.df.degs, "saved_objects/idc-degs-p.05-LFC1.csv")

res.idc.degs <- res.idc.complete[res.idc.complete$padj<0.05 & abs(res.idc.complete$log2FoldChange)>1,]
summary(res.idc.degs)

## Distribution of the p-adjusted value
summary(res.idc.degs$padj)
hist(res.idc.degs$padj, main="Distribution of the adjusted p-value in the DEGs of IDC", xlab="Adjusted p-value")
################################################################################
## Repair genes subsetting 
res.idc.rep <- res.idc[rownames(res.idc) %in% repair.genes$symbol,]
summary(res.idc.rep)
## out of 285 genes, 0 up-regulated, 36 down-regulated
res.idc.degs.rep <- res.idc.degs[rownames(res.idc.degs) %in% repair.genes$symbol,]
summary(res.idc.degs.rep)

## Distribution of the LFC and p-adjusted value
summary(res.idc.rep$log2FoldChange)
hist(res.idc.rep$log2FoldChange, main="Distribution of the Log2 fold change in the results (Repair Genes Only) of IDC", xlab="Log2 Fold Change")
summary(res.idc.rep$padj)
hist(res.idc.rep$padj, main="Distribution of the adjusted p-value in the results (Repair Genes Only) of IDC", xlab="Adjusted p-value")

summary(res.idc.degs.rep$padj)
hist(res.idc.degs.rep$padj, main="Distribution of the adjusted p-value in the Repair DEGs of IDC", xlab="Adjusted p-value")
################################################################################
## Plotting: 1. MAplot
plotMA(res.idc, alpha=0.05, main="MA plot of the not-mormalized DESeq Results of IDC")
################################################################################
## Plotting: 2. EnhancedVolcano Plot
EnhancedVolcano(res.idc, 
                lab = rownames(res.idc), 
                x = 'log2FoldChange', 
                y = 'pvalue',
                title = "Volcano plot of the not-normalized DESeq results",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                legendPosition = "right")
################################################################################
## subset vsd
dds.vsd.idc <- dds.vsd[,dds.vsd$group %in% c("IDC_Normal", "IDC_Tumor")]
dds.vsd.idc
table(dds.vsd.idc$group)
dds.vsd.idc$group <- factor(dds.vsd.idc$group)
table(dds.vsd.idc$group)

## DEGs
dds.vsd.idc.degs <- dds.vsd.idc[rownames(dds.vsd.idc) %in% rownames(res.idc.degs),]
dds.vsd.idc.degs

## Repair DEGs
dds.vsd.idc.degs.rep <- dds.vsd.idc.degs[rownames(dds.vsd.idc.degs) %in% repair.genes$symbol ,]
dds.vsd.idc.degs.rep
################################################################################
## Plotting: 3. PCA
pca.plot.vsd.idc <- plotPCA(dds.vsd.idc, intgroup="group", ntop = 25531, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.idc, "percentVar"))
ggplot(pca.plot.vsd.idc, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vst-normalized DESeq DataSet (ntop = 25531)")
rm(percentVar)
################################################################################
## Plotting: 4. Dispersion estimate
plotDispEsts(dds.run.idc)
################################################################################
## Plotting: 7. Boxplot with dots
source('scripts/plotting_gene.R')
genes.list <- rownames(res.idc.degs.rep)
par(mfrow=c(3,3))
for (gen in genes.list){
  plotting_gene(dds.run.idc, gen)
}
par(mfrow=c(1,1))
################################################################################
## Plotting: 8. Heatmap & Dendogram
colors2 <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
heatmap.2(assay(dds.vsd.idc.degs.rep), col=colors2,
          scale="row", trace="none", labCol=substring(colnames(dds.vsd.idc.degs.rep), 6),
          # dendrogram = "row",
          Colv=order(dds.vsd.idc.degs.rep$Sample.Type), Rowv=TRUE,
          ColSideColors = c(Primary.Tumor="darkgreen", Solid.Tissue.Normal="orange")[colData(dds.vsd.idc.degs.rep)$Sample.Type],
          key =TRUE, key.title="Heatmap Key",
          main="Heatmap of the top 36 differentially expressed repair genes in IDC")
