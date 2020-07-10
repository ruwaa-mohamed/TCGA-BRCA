################################################################################
################################ Session Info. #################################
## R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
## Copyright (C) 2020 The R Foundation for Statistical Computing
## Platform: x86_64-pc-linux-gnu (64-bit)

## Created by Ruwaa I. Mohamed
## Date: Saturday, May 30, 2020.
################################################################################
## Installing Bioconductor and the required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("DESeq2")              # packageVersion: 1.28.1
BiocManager::install('EnhancedVolcano')     # packageVersion: 1.6.0

install.packages("RColorBrewer")  # packageVersion: 1.1.2
install.packages("gplots")        # packageVersion: 3.0.3
################################################################################
## Loading the required libraries/scripts
source('scripts/get_df_from_gdc_files.R')

library(DESeq2)
library(EnhancedVolcano)

library(RColorBrewer)
library(gplots)
################################################################################
## Reading the sample sheet
#mirna.sample_sheet <- read.csv("saved_objects/mirna.sample_sheet.csv", header=TRUE)
mirna.sample_sheet <- read.csv("raw_data/gdc_sample_sheet.2020-07-10.tsv", sep="\t", header=TRUE)
mirna.sample_sheet <- mirna.sample_sheet[! mirna.sample_sheet$Sample.Type == "Metastatic",]
mirna.sample_sheet$Sample.Type <- as.factor(mirna.sample_sheet$Sample.Type)
mirna.sample_sheet$Sample.Type <- relevel(mirna.sample_sheet$Sample.Type, "Solid Tissue Normal")
levels(mirna.sample_sheet$Sample.Type) <- c("Normal", "Tumor")
table(mirna.sample_sheet$Sample.Type)
################################################################################
## Reading the files
mirna.exp.df <- get_df_from_gdc_mirna(mirna.path="raw_data/miRNA-seq", file_ext=".mirnas.quantification.txt$")
mirna.exp.df.sub <- mirna.exp.df[, colnames(mirna.exp.df) %in% mirna.sample_sheet$File.ID]

# comparing with the RNA-seq data
summary(mirna.sample_sheet$Case.ID %in% sample_sheet$Case.ID)
file.ids <- mirna.sample_sheet[mirna.sample_sheet$Case.ID %in% sample_sheet$Case.ID,]$File.ID
mirna.exp.df.sub <- mirna.exp.df.sub[, colnames(mirna.exp.df.sub) %in% file.ids]
mirna.sample_sheet <- mirna.sample_sheet[mirna.sample_sheet$File.ID %in% file.ids,]
################################################################################
## Preparing DESeq2 Data
count.data.mirna <- mirna.exp.df.sub

col.data.mirna <- mirna.sample_sheet[, colnames(mirna.sample_sheet) %in% c("File.ID", "Sample.Type", "Sample.ID")]

# reording 
new.order <- match(colnames(count.data.mirna), col.data.mirna$File.ID)
count.data.mirna <- count.data.mirna[order(new.order)]
rm(new.order)

# convert to matrix (numbers only)
count.data.mirna <- apply (count.data.mirna, 2, as.integer)
rownames(count.data.mirna) <- rownames(mirna.exp.df.sub)

# check the final order match
all(colnames(count.data.mirna) == col.data.mirna$File.ID)
################################################################################
## Building DESeqDataSet
dds.mirna <- DESeqDataSetFromMatrix(countData=count.data.mirna , colData=col.data.mirna , design=~Sample.Type)
dds.mirna
saveRDS(dds.mirna, "saved_objects/dds.mirna.rds")
dds.mirna <- readRDS("saved_objects/dds.mirna.rds")

## Collapse Technical Replicates
ddsCollapsed.mirna <- collapseReplicates(dds.mirna, groupby=dds.mirna$Sample.ID)
ddsCollapsed.mirna
saveRDS(ddsCollapsed.mirna, "saved_objects/ddsCollapsed.mirna.rds")
ddsCollapsed.mirna <- readRDS("saved_objects/ddsCollapsed.mirna.rds")

## to check that collapsing worked!
original <- rowSums(counts(dds.mirna)[, dds.mirna$Sample.ID == "TCGA-A7-A0DB-01A"])
all(original == counts(ddsCollapsed.mirna)[,"TCGA-A7-A0DB-01A"])
rm(original)
################################################################################
## Run DESEQ2
dds.mirna.run <- DESeq(ddsCollapsed.mirna)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 109 genes
# -- DESeq argument 'minReplicatesForReplace' = 7
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
dds.mirna.run

saveRDS(dds.mirna.run, file="saved_objects/dds.mirna.run.rds")
dds.mirna.run <- readRDS("saved_objects/dds.mirna.run.rds")
################################################################################
## Creating the results (RES) object
res.mirna <- results(dds.mirna.run, contrast=c("Sample.Type", "Tumor", "Normal"), alpha=0.05,lfcThreshold=1)
summary(res.mirna)
# out of 1600 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 1.00 (up)    : 87, 5.4%
# LFC < -1.00 (down) : 64, 4%
# outliers [1]       : 0, 0%
# low counts [2]     : 676, 42%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

## Removing incomplete genes (No need for this step!)
res.mirna.2 <- res.mirna[complete.cases(res.mirna),]
summary(res.mirna.2)

## saving to and reading from RDS object
saveRDS(res.mirna, file="saved_objects/res.mirna.rds")
res.mirna <- readRDS("saved_objects/res.mirna.rds")

## Chkeckin the distribution of the LFC and p-adjusted value
summary(res.mirna$log2FoldChange)
hist(res.mirna$log2FoldChange, main="Distribution of the Log2 fold change in the miRNA results", xlab="Log2 Fold Change")
summary(res.mirna$padj)
hist(res.mirna$padj, main="Distribution of the adjusted p-value in the miRNA results", xlab="Adjusted p-value")
################################################################################
## Getting the DEMs
res.mirna.dems <- res.mirna.2[abs(res.mirna.2$log2FoldChange)>1 & res.mirna.2$padj<0.05,]
summary(res.mirna.dems)

summary(res.mirna.dems$padj)
hist(res.mirna.dems$padj, main="Distribution of the adjusted p-value in the miRNA DEGs", xlab="Adjusted p-value")

write.csv(res.mirna.dems, file="saved_objects/mirna.dems.all.csv", quote=FALSE)
write.table(rownames(res.mirna.dems), file="saved_objects/res.mirna.dems.tsv", row.names = FALSE, col.names=FALSE, quote=FALSE)

res.mirna.dems.df <- as.data.frame(res.mirna.dems)
res.mirna.dems.df <- res.mirna.dems.df[order(res.mirna.dems.df$padj),]
################################################################################
## Data Normalization for plotting: 1. VST Normalization
dds.mirna.vsd <- varianceStabilizingTransformation(dds.mirna.run, blind=FALSE)
dds.mirna.vsd

dds.mirna.vsd.dems <- dds.mirna.vsd[rownames(dds.mirna.vsd) %in% rownames(res.mirna.dems),]
dds.mirna.vsd.dems
################################################################################
## Plotting: 1. MA plot
plotMA(res.mirna, alpha=0.05, main="MA plot of the not-mormalized miRNA DESeq Results")
################################################################################
## Plotting: 2. EnhancedVolcano
EnhancedVolcano(res.mirna, 
                lab = rownames(res.mirna), 
                x = 'log2FoldChange', 
                y = 'pvalue',
                title = "Volcano plot of the not-normalized miRNA DESeq results",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                legendPosition = "right")
################################################################################
## Plotting: 3. PCA
pca.plot.mirna.vsd <- plotPCA(dds.mirna.vsd, intgroup="Sample.Type", ntop = 1881, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.mirna.vsd, "percentVar"))
ggplot(pca.plot.mirna.vsd, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vsd-normalized DESeq DataSet of the miRNA data (ntop = 1881)")
rm(percentVar)

pca.plot.mirna.vsd.degs <- plotPCA(dds.mirna.vsd.dems, intgroup="Sample.Type", ntop = 151, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.mirna.vsd.degs, "percentVar"))
ggplot(pca.plot.mirna.vsd.degs, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of the vsd-normalized DESeq DEGs of the miRNA data (ntop = 151)")
rm(percentVar)
################################################################################
## Plotting: 3. Heatmap
colors2 <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
heatmap.2(assay(dds.mirna.vsd.dems), col=colors2,
          scale="row", trace="none", labCol=substring(colnames(dds.mirna.vsd.dems), 6),
          # dendrogram = "row",
          Colv=order(dds.mirna.vsd.dems$Sample.Type), Rowv=TRUE,
          ColSideColors = c(Tumor="darkgreen", Normal="orange")[colData(dds.mirna.vsd.dems)$Sample.Type],
          key =TRUE, key.title="Heatmap Key",
          main="Heatmap of the 151 differentially expressed miRNAs")
################################################################################
## Plotting: 4. Dispersion estimate
plotDispEsts(dds.mirna.run)
################################################################################
################################################################################
# ## miRDB Data Download (http://mirdb.org/download.html)
# miRDB <- read.table("raw_data/miRNA-DBs/miRDB_v6.0_prediction_result.txt", header=FALSE, sep="\t")
# miRDB.hsa <- miRDB[startsWith(miRDB$V1, "hsa-"),]
# 
# miRDB.mirna.dems <- miRDB.hsa[startsWith(miRDB$V1, rownames(res.mirna.dems)[1]),]
# for (i in 2:length(rownames(res.mirna.dems))){
#   miRDB.temp <- miRDB.hsa[startsWith(miRDB$V1, rownames(res.mirna.dems)[i]),]
#   miRDB.mirna.dems <- rbind(miRDB.mirna.dems, miRDB.temp)
# }
# 
# ## only 25 hsa were matched miRWalk_miRNA_Targets!
# miRWalk_miRNA_Targets <- read.csv("raw_data/miRNA-DBs/miRWalk_miRNA_Targets.csv")
# rna.res.degs.mirna.dems.miRWalk <- res.degs[rownames(res.degs) %in% miRWalk_miRNA_Targets$genesymbol,]
# summary(rna.res.degs.mirna.dems.miRWalk)
# rna.res.degs.mirna.dems.rep.miRWalk <- rna.res.degs.mirna.dems.miRWalk[rownames(rna.res.degs.mirna.dems.miRWalk) %in% repair.genes$symbol,]
# summary(rna.res.degs.mirna.dems.rep.miRWalk)
# 
# rna.res.degs.rep.mirna.dems.rep.miRWalk <- rna.res.degs.mirna.dems.miRWalk[rownames(rna.res.degs.mirna.dems.miRWalk) %in% rownames(res.rep.degs),]
# summary(rna.res.degs.rep.mirna.dems.rep.miRWalk)
################################################################################
## mirDIP : microRNA Data Integration Portal
mirDIP.inputVerification <- read.table("raw_data/miRNA-DBs/mirDIP_E_2020_07_09_19_23_48.txt", sep="\t", skip=163, nrow=124, header=TRUE)
mirDIP.Results <- read.table("raw_data/miRNA-DBs/mirDIP_E_2020_07_09_19_23_48.txt", sep="\t", skip=293, header=TRUE)

# mir is written with capital R in the results and small r in the dems!
substring(mirDIP.Results$MicroRNA, 7) = 'r'

## subsetting DEMs
mirDIP.dems.up <- mirDIP.Results[mirDIP.Results$MicroRNA %in% rownames(res.mirna.dems)[res.mirna.dems$log2FoldChange>0],]
mirDIP.dems.dwn <- mirDIP.Results[mirDIP.Results$MicroRNA %in% rownames(res.mirna.dems)[res.mirna.dems$log2FoldChange<0],]

## subsetting DEMs on degs 
mirDIP.dems.up.degs.up <- mirDIP.dems.up[mirDIP.dems.up$Gene.Symbol %in% rownames(res.degs)[res.degs$log2FoldChange>0],]
mirDIP.dems.up.degs.dwn <- mirDIP.dems.up[mirDIP.dems.up$Gene.Symbol %in% rownames(res.degs)[res.degs$log2FoldChange<0],]

mirDIP.dems.dwn.degs.up <- mirDIP.dems.dwn[mirDIP.dems.dwn$Gene.Symbol %in% rownames(res.degs)[res.degs$log2FoldChange>0],]
mirDIP.dems.dwn.degs.dwn <- mirDIP.dems.dwn[mirDIP.dems.dwn$Gene.Symbol %in% rownames(res.degs)[res.degs$log2FoldChange<0],]

## getting DEMs on repair genes only
mirDIP.Results.rep.all <- mirDIP.Results[mirDIP.Results$Gene.Symbol %in% repair.genes$symbol,]
mirDIP.Results.rep.degs <- mirDIP.Results[mirDIP.Results$Gene.Symbol %in% rownames(res.rep.degs),]

## chechking for negative and positive correlations
table(unique(mirDIP.Results.rep.degs$MicroRNA) %in% mirDIP.dems.dwn$MicroRNA)
table(unique(mirDIP.Results.rep.degs$MicroRNA) %in% mirDIP.dems.up$MicroRNA)
# there are 8 miRNA with negative correlation and 12 with positive correlation! (we need the negative!)
# on gene level, 17 repair DEGs with negative correlation and 22 with positive correlation
################################################################################