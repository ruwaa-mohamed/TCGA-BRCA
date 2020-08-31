################################################################################
### RNA - IDC_Tumor vs. LC_Tumor
dds.run <- readRDS("saved_objects/dds.run.group.rds")
dds.vsd <- readRDS("saved_objects/dds.vsd.rds")
res.idc.degs.rep <- read.csv("saved_objects/res.idc.degs.rep.csv", header = TRUE)
################################################################################
## Libraries
library(DESeq2)
library(EnhancedVolcano)
library(biomaRt)
library(dplyr)
library(RColorBrewer)
library(gplots)

source("scripts/proteins.R")
source("scripts/plotting_gene_2.R")
source("scripts/to_tiff.R")
################################################################################
## Extracting dds.run.idclc from dds.run
dds.run.idclc <- dds.run[,dds.run$group %in% c("IDC_Tumor", "LC_Tumor")]
dds.run.idclc
# 25,531 gene, 969 samples

dds.run.idclc$group <- droplevels(dds.run.idclc$group)
table(dds.run.idclc$group)
# 969 samples = 767 IDC + 202 LC
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
################################################################################
## complete cases only
res.idclc.complete <- res.idclc[complete.cases(res.idclc),]
################################################################################
## Extracting DEGs
res.idclc.degs <- res.idclc.complete[res.idclc.complete$padj<0.05 & abs(res.idclc.complete$log2FoldChange)>1,]
summary(res.idclc.degs)
# 663 DEGs: 591 up, 72 down

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

res.idclc.idc.degs.rep <- res.idclc.complete[rownames(res.idclc.complete) %in% rownames(res.idc.degs.rep),]
summary(res.idclc.idc.degs.rep)

## Distribution of the LFC and p-adjusted value
summary(res.idclc.rep$log2FoldChange)
to_tiff(
  hist(res.idclc.rep$log2FoldChange, main="Log2 Fold Change of All Repair Genes between IDC and LC", xlab="Log2 Fold Change"),
  "rna.idclc.rep.all.l2fc.hist.tiff")
################################################################################
## Exporting to files
write.table(as.data.frame(res.idclc.degs), "saved_objects/res.idclc.degs.csv", sep=",", quote=FALSE)
write.table(as.data.frame(res.idclc.idc.degs.rep), "saved_objects/res.idclc.idc.degs.rep.csv", sep=",", quote=FALSE)
################################################################################
## Biomart part
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# filters <- listFilters(ensembl)
# attributes <- listAttributes(ensembl)
mart <- getBM(attributes = c("external_gene_name", "description", "name_1006"),
              filters = "external_gene_name",
              values = c(rownames(res.idclc.degs), rownames(res.idc.degs.rep)), 
              mart = ensembl)
mart$name_1006[mart$name_1006 == ""] <- NA
mart <- mart[complete.cases(mart),]
mart.merged <- mart %>% 
  group_by(external_gene_name, description) %>%
  summarize(name_1006 = paste(name_1006, collapse="; "))

## Adding mart to the DEGs
res.idclc.degs.desc <- cbind(Gene.Name = rownames(res.idclc.degs), as.data.frame(res.idclc.degs))
res.idclc.degs.desc <- left_join(res.idclc.degs.desc, mart.merged,
                               by=c("Gene.Name" = "external_gene_name"))
res.idclc.degs.desc <- res.idclc.degs.desc %>% rename(Go.term.name = name_1006)

## Adding mart to Repair DEGs
res.idclc.degs.rep.desc <- cbind(Gene.Name = rownames(res.idc.degs.rep), as.data.frame(res.idc.degs.rep))
res.idclc.degs.rep.desc <- left_join(res.idclc.degs.rep.desc, mart.merged, 
                                   by = c("Gene.Name" = "external_gene_name"))
res.idclc.degs.rep.desc <- res.idclc.degs.rep.desc %>% rename(Go.term.name = name_1006)

## Exporting the results
write.table(res.idclc.degs.desc, "saved_objects/res.idclc.degs.description.csv", sep=",", quote=FALSE)
write.table(res.idclc.degs.rep.desc, "saved_objects/res.idclc.degs.rep.description.csv", sep=",", quote=FALSE)
################################################################################
## Plotting: 1. MAplot
to_tiff(
  plotMA(res.idclc.rep, alpha=0.05, main="MAplot of Repair Genes between IDC and LC"),
  "rna.idclc.rep.MAplot.tiff")
################################################################################
## Plotting: 2. EnhancedVolcano Plot
tiff("final_figures/rna.idclc.rep.EnhancedVolcano.tiff", width = 18, height = 21, units = 'cm', res = 300)
EnhancedVolcano(res.idclc.rep, 
                lab = rownames(res.idclc.rep), 
                x = 'log2FoldChange', 
                y = 'padj',
                title = "Repair Genes between IDC and LC",
                pCutoff = 0.05,
                FCcutoff = 1.0)
dev.off()
################################################################################
## subset vsd
dds.vsd.idclc <- dds.vsd[,dds.vsd$group %in% c("IDC_Tumor", "LC_Tumor")]
dds.vsd.idclc
dds.vsd.idclc$group <- droplevels(dds.vsd.idclc$group)
table(dds.vsd.idclc$group)

## All Repair
dds.vsd.idclc.rep <- dds.vsd.idclc[rownames(dds.vsd.idclc) %in% repair.genes$symbol,]
dds.vsd.idclc.rep

## DEGs
dds.vsd.idclc.degs <- dds.vsd.idclc[rownames(dds.vsd.idclc) %in% rownames(res.idclc.degs),]
dds.vsd.idclc.degs

## Repair IDC DEGs
dds.vsd.idclc.degs.rep <- dds.vsd.idclc[rownames(dds.vsd.idclc) %in% rownames(res.idc.degs.rep),]
dds.vsd.idclc.degs.rep
################################################################################
## Plotting: 3. PCA
# All Genes
tiff("final_figures/rna.idclc.pca.all.tiff", width = 18, height = 21, units = 'cm', res = 300)
pca.plot.vsd.idclc <- plotPCA(dds.vsd.idclc, intgroup="group", ntop = 25531, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.idclc, "percentVar"))
ggplot(pca.plot.vsd.idclc, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of All (25,531) Genes between IDC and LC")
rm(pca.plot.vsd.idclc, percentVar)
dev.off()

# All Repair
tiff("final_figures/rna.idclc.pca.all.rep.tiff", width = 18, height = 21, units = 'cm', res = 300)
pca.plot.vsd.idclc.rep <- plotPCA(dds.vsd.idclc.rep, intgroup="group", ntop = 285, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.idclc.rep, "percentVar"))
ggplot(pca.plot.vsd.idclc.rep, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of Repair (285) Genes between IDC and LC")
rm(pca.plot.vsd.idclc.rep, percentVar)
dev.off()

# All DEs
tiff("final_figures/rna.idclc.pca.degs.tiff", width = 18, height = 21, units = 'cm', res = 300)
pca.plot.vsd.idclc.degs <- plotPCA(dds.vsd.idclc.degs, intgroup="group", ntop = 663, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.idclc.degs, "percentVar"))
ggplot(pca.plot.vsd.idclc.degs, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of all DEGs (663) between IDC and LC")
rm(pca.plot.vsd.idclc.degs, percentVar)
dev.off()

# Repair DEGs
tiff("final_figures/rna.idclc.pca.degs.rep.tiff", width = 18, height = 21, units = 'cm', res = 300)
pca.plot.vsd.idclc.degs.rep <- plotPCA(dds.vsd.idclc.degs.rep, intgroup="group", ntop = 36, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.idclc.degs.rep, "percentVar"))
ggplot(pca.plot.vsd.idclc.degs.rep, aes(PC1, PC2, color=group)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of IDC Repair DEGs (36) between IDC and LC")
rm(pca.plot.vsd.idclc.degs.rep, percentVar)
dev.off()
################################################################################
## Plotting: 4. Dispersion estimate
to_tiff(
  plotDispEsts(dds.run.idclc, main="Disperssion Estimates for all Genes between IDC and LC"),
  "rna.idclc.plotDispEsts.tiff")
################################################################################
## Plotting: 7. Boxplot with dots
genes.list <- rownames(res.idc.degs.rep)

pdf(file="final_figures/rna.idclc.vsd.rep.degs.idc.boxplot.pdf", onefile=TRUE, paper="a4", width = 8, height = 11)
par(mfrow=c(3,3))
for (gen in genes.list){
  plotting_gene_2(dds.vsd.idclc.degs.rep, gen, "", atr=1)
}
dev.off()

par(mfrow=c(1,1))
rm(gen, genes.list)
################################################################################
## Plotting: 8. Heatmap & Dendogram
pdf(file="final_figures/rna.idclc.heatmap.rep.degs.idc.pdf", onefile=TRUE, paper="a4r", width = 11, height = 8)

colors2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)
heatmap.2(assay(dds.vsd.idclc.degs.rep), col=colors2,
          scale="row", trace="none", labCol=substring(colnames(dds.vsd.idclc.degs.rep), 6),
          # dendrogram = "row",
          Colv=order(dds.vsd.idclc.degs.rep$Sample.Type), Rowv=TRUE,
          ColSideColors = c(IDC="darkgreen", LC="orange")[colData(dds.vsd.idclc.degs.rep)$primary_diagnosis],
          key =TRUE, key.title="Heatmap Key",
          main="Heatmap of the 36 Repair IDC DEGs between IDC (darkgreen) and LC (orange)")
dev.off()
################################################################################