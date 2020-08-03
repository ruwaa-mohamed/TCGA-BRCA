################################################################################
################################################################################
## R version 4.0.1 (2020-06-06) -- "See Things Now"
## Copyright (C) 2020 The R Foundation for Statistical Computing
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04 LTS

## BiocManager version: 1.30.10
## Bioconductor version: 3.11
## Date: Wednesday, June 17, 2020.
## Created by Ruwaa I. Mohamed
################################################################################
################################################################################
## Installing Bioconductor and the required libraries
# use packageVersion('') to check the versio of all installed packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::version()
install.packages("BiocManager")

BiocManager::install("DESeq2")              # packageVersion: 1.28.1
BiocManager::install("org.Hs.eg.db")        # packageVersion: 3.11.4
BiocManager::install('EnhancedVolcano')     # packageVersion: 1.6.0
# BiocManager::install("genefilter")        # packageVersion: 1.70.0    # already installed with DESeq2
BiocManager::install("apeglm")              # packageVersion: 1.10.0

install.packages("stringr")       # packageVersion: 1.4.0
install.packages("gplots")        # packageVersion: 3.0.3
install.packages("ggplot2")       # packageVersion: 3.3.2
install.packages("RColorBrewer")  # packageVersion: 1.1.2
# install.packages("cluster")     # packageVersion: 2.1.0
# install.packages("biclust")     # packageVersion: 2.0.2
# install.packages("ggdendro")    # packageVersion: 0.1.20
################################################################################
## Loading the required libraries/scripts
source('scripts/get_df_from_gdc.R')
source('scripts/ENSG_to_symbol.R')
source('scripts/aggregate_rows.R')
source('scripts/proteins.R')
source('scripts/plotting_gene.R')
source('scripts/to_tiff.R')

library(gplots)
library(ggplot2)
library(RColorBrewer)
# library(cluster)
# library(biclust)
# library(ggdendro)

library(DESeq2)
library(EnhancedVolcano)
# library(apeglm)
# library(genefilter)
################################################################################
## Reading the RNA-seq files (1,222 files)
rna.data.path <- "raw_data/RNA-seq/"

## reading the expression values from the raw data and saving them to RDS object (Don't use this is unless you don't have the RDS object)!
rna.exp.df <- get_df_from_gdc(rna.data.path, 'htseq.counts.gz$')
rna.exp.df <- rna.exp.df[1:(nrow(rna.exp.df)-5),]
saveRDS(rna.exp.df, file = "saved_objects/exp.rna.rds")

## reading the expression values from the previously saved RDS object (use this instead)!
rna.exp.df <- readRDS(file = "saved_objects/exp.rna.rds")
View(rna.exp.df) # DON'T!
################################################################################
## Reading the sample sheet
sample.sheet <- read.csv("raw_data/females_only.ductal-and-lobular-only.gdc_sample_sheet.2020-06-17.tsv", sep="\t", header=TRUE)
table(sample.sheet$Sample.Type)

## dropping metastatic samples from the sample sheet
sample.sheet <- sample.sheet[! sample.sheet$Sample.Type == "Metastatic",]
table(sample.sheet$Sample.Type)

## dropping metastatic samples from the RNA-seq expression dataframe
rna.exp.df <- rna.exp.df[, colnames(rna.exp.df) %in% sample.sheet$File.ID]

saveRDS(rna.exp.df, file = "saved_objects/exp.rna.subsetted.rds")
rna.exp.df <- readRDS("saved_objects/exp.rna.subsetted.rds")
################################################################################
## Reading the clinical data 
clinical <- read.csv("raw_data/clinical.project-TCGA-BRCA.2020-03-05/clinical.tsv", header=TRUE, sep="\t", na="--")
clinical <- clinical[, colSums(is.na(clinical)) != nrow(clinical)]
clinical <- clinical[clinical$submitter_id %in% sample.sheet$Case.ID ,]
write.csv(clinical, 'saved_objects/clinical.all.csv', row.names=FALSE)
write.csv(clinical[clinical$submitter_id %in% sample.sheet$Case.ID,], 'saved_objects/clinical_filtered.csv', row.names=FALSE)
################################################################################
## subsetting from clinical and mixing with sample sheet
clinical <- clinical[clinical$primary_diagnosis %in% c("Infiltrating duct and lobular carcinoma", "Infiltrating duct carcinoma, NOS", "Lobular carcinoma, NOS"),]
clinical <- unique(clinical[,c("submitter_id", "primary_diagnosis")])
mapperIDs <- match(sample.sheet$Case.ID, clinical$submitter_id)
sample.sheet.clinical <- cbind(sample.sheet, primary_diagnosis=clinical$primary_diagnosis[mapperIDs])
sample.sheet.clinical <- sample.sheet.clinical[! is.na(sample.sheet.clinical$primary_diagnosis),]
rm(mapperIDs)
sample.sheet <- sample.sheet.clinical

write.csv(sample.sheet, "saved_objects/sample.sheet.csv", row.names=FALSE)
sample.sheet <- read.csv("saved_objects/sample.sheet.csv", header=TRUE)
write.csv(clinical, "saved_objects/clinical.csv", row.names=FALSE)
clinical <- read.csv("saved_objects/clinical.csv", header=TRUE)
## 992 individual, 2 columns

# drop from the RNA-exp DF too
rna.exp.df <- rna.exp.df[, colnames(rna.exp.df) %in% sample.sheet$File.ID]
################################################################################
## Changing ENSG to SYMBOL
rna.exp.df.sym <- ENSG_to_symbol(rna.exp.df)
saveRDS(rna.exp.df.sym, file = "saved_objects/rna.exp.df.sym.rds")
rna.exp.df.sym <- readRDS(file = "saved_objects/rna.exp.df.sym.rds")
## 60,483 ENSG where reduced to 25,596 gene symbol (org.Hs.eg.db version 3.11.4)

rna.exp.df.sym.agg <- aggregate_rows(rna.exp.df.sym, agg.var=rna.exp.df.sym$symbol) # This step takes too long to run!
saveRDS(rna.exp.df.sym.agg, file = "saved_objects/rna.exp.df.sym.agg.rds")
rna.exp.df.sym.agg <- readRDS(file="saved_objects/rna.exp.df.sym.agg.rds")
## They're aggregated to 25,531 genes
rm(rna.exp.df.sym)
################################################################################
## Preparing the inputs of DESEQ2
count.data <- rna.exp.df.sym.agg

col.data <- sample.sheet[, colnames(sample.sheet) %in% c("File.ID", "Sample.ID", "Sample.Type", "primary_diagnosis")]

col.data$Sample.Type <- as.factor(col.data$Sample.Type)
col.data$Sample.Type <- relevel(col.data$Sample.Type, "Solid Tissue Normal")
levels(col.data$Sample.Type) <- c("Normal", "Tumor")
table(col.data$Sample.Type)
## 105 Normal and 1,002 Tumor

col.data$primary_diagnosis <- as.factor(col.data$primary_diagnosis)
levels(col.data$primary_diagnosis) <- c("Mixed", "IDC", "LC")
table(col.data$primary_diagnosis)
## 37 Mixed, 86 IDC, and 210 LC

# reording 
new.order <- match(colnames(count.data), col.data$File.ID)
count.data <- count.data[order(new.order)]
rm(new.order)

# convert to matrix (numbers only)
count.data <- apply (count.data, 2, as.integer)
rownames(count.data) <- rownames(rna.exp.df.sym.agg)

# check the final order match
all(colnames(count.data) == col.data$File.ID)
################################################################################
## Building DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=count.data , colData=col.data , design=~Sample.Type) #primary_diagnosis+Sample.Type+primary_diagnosis:Sample.Type)
dds

dds$group <- factor(paste(dds$primary_diagnosis, dds$Sample.Type, sep="_"))
design(dds) <- ~ group
dds

table(dds$group)
## IDC_ Normal: 89, IDC_Tumor: 771
## LC_Normal: 7, LC_Tumor:203
## Mixed_Normal: 9, Mixed_Tumor: 28

saveRDS(dds, "saved_objects/dds.rds")
dds <- readRDS("saved_objects/dds.rds")

## Collapse Technical Replicates
ddsCollapsed <- collapseReplicates(dds, groupby=dds$Sample.ID)
ddsCollapsed
saveRDS(ddsCollapsed, "saved_objects/ddsCollapsed.rds")
ddsCollapsed <- readRDS("saved_objects/ddsCollapsed.rds")

## to check that collapsing worked!
original <- rowSums(counts(dds)[, dds$Sample.ID == "TCGA-A7-A26J-01A"])
all(original == counts(ddsCollapsed)[,"TCGA-A7-A26J-01A"])
rm(original)

## testing
table(ddsCollapsed$group)
## IDC_ Normal: 89, IDC_Tumor: 767 (-4)
## LC_Normal: 7, LC_Tumor:202 (-1)
## Mixed_Normal: 9, Mixed_Tumor: 28
################################################################################
################################################################################
## Run DESEQ2: on the 1107 samples (on the new 'group' column)
dds.run <- DESeq(ddsCollapsed)
# The following steps were run:
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 3267 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
dds.run
# resultsNames(dds.run)
saveRDS(dds.run, file="saved_objects/dds.run.rds")
dds.run <- readRDS("saved_objects/dds.run.rds")
################################################################################
################################################################################
## To subset based on BC subtype, run the following section

## Data Normalization for plotting: 1. VST Normalization
dds.vsd <- varianceStabilizingTransformation(dds.run, blind=FALSE)
dds.vsd

saveRDS(dds.vsd, "saved_objects/dds.vsd.rds")
dds.vsd <- readRDS("saved_objects/dds.vsd.rds")

source('scripts/idc.R')
source('scripts/lc.R')
source('scripts/mixed.R')

## For P-val 0.05 and LFC 1, out of 25,175 genes
## IDC: 2,261 up-regulated, 1,619 down-regulated
## LC: 162 up-regulated, 119 down-regulated
## Mixed: 387 up-regulated, 620 down-regulated
################################################################################
################################################################################
## continue the rest of this scipt only if you want to ignore BC subtyps and 
## deal with tumor vs. normal only

## Adjusting Design and Rerunning
design(ddsCollapsed) <- ~Sample.Type
dds.run.TN <- DESeq(ddsCollapsed)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 3906 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing

saveRDS(dds.run.TN, "saved_objects/dds.run.TN.rds")
dds.run.TN <- readRDS("saved_objects/dds.run.TN.rds")

## Creating the results (RES) object
res <- results(dds.run.TN, contrast=c("Sample.Type", "Tumor", "Normal"), alpha=0.05,lfcThreshold=1)
summary(res)
# out of 25,175 genes: 2,220 up-regulated and 1,608 down-regulated

## Removing incomplete genes (No need for this step!)
res.complete <- res[complete.cases(res),]
summary(res.complete)

## saving to and reading from RDS object
saveRDS(res, file="saved_objects/res.rds")
res <- readRDS("saved_objects/res.rds")

## Chkeckin the distribution of the LFC and p-adjusted value
summary(res$log2FoldChange)
# pdf(paste(figures_dir, "res_log2FoldChange.pdf", sep="/"))
to_tiff(
  hist(res$log2FoldChange, main="Distribution of the Log2 Fold Change", xlab="Log2 Fold Change"),
  "res_log2FoldChange.tiff")

summary(res$padj)
to_tiff(
  hist(res$padj, main="Distribution of the Adjusted p-value", xlab="Adjusted p-value"),
  "res_Adjusted_p-value.tiff")
################################################################################
## Converting the results to DF
res.df <- as.data.frame(res.complete)
res.df <- res.df[order(res.df$padj),]
write.csv(res.df, "saved_objects/all.results.sorted.csv")
rm(res.df)
################################################################################
## Extracting DEGs
res.degs <- res.complete[res.complete$padj<0.05 & abs(res.complete$log2FoldChange)>1,]
summary(res.degs)
write.csv(as.data.frame(res.degs), "saved_objects/degs.csv")

summary(res.degs$padj)
to_tiff(
  hist(res.degs$padj, main="Distribution of the adjusted p-value in the DEGs", xlab="Adjusted p-value"),
  "res.degs.Adjusted_p-value.tiff")
################################################################################
## Repair genes subsetting
res.rep <- res[rownames(res) %in% repair.genes$symbol,]
summary(res.rep)
## 285 genes, 34 up-regulated, 0 down-regulated

## Extracting DEGs
res.degs.rep <- res.degs[rownames(res.degs) %in% repair.genes$symbol,]
summary(res.degs.rep)

## Chkeckin the distribution of the LFC and p-adjusted value
summary(res.rep$log2FoldChange)
to_tiff(
  hist(res.rep$log2FoldChange, main="Distribution of the Log2 Fold Change (Repair Genes Only)", xlab="Log2 Fold Change"),
  "res.rep.Hist_L2FC.tiff")
summary(res.rep$padj)
to_tiff(
  hist(res.rep$padj, main="Distribution of the Adjusted p-value (Repair Genes Only)", xlab="Adjusted p-value"),
  "res.rep.Hist_Adjusted_p-value.tiff")

summary(res.degs.rep$log2FoldChange)
# hist(res.degs.rep$log2FoldChange, main="Distribution of the Log2 fold change in the DEGs (Repair Genes Only)", xlab="Log2 Fold Change")
summary(res.degs.rep$padj)
to_tiff(
  hist(res.degs.rep$padj, main="Distribution of the Adjusted p-value in Repair DEGs.", xlab="Adjusted p-value"),
  "res.rep.degs.Hist_Adjusted_p-value.tiff")
################################################################################
## Normalization Methods Available:
## 1. rlog()  --> takes too long and not applicable for our dataset (not as sensetive as vst and takes too long and fails)
## 2. vst() or varianceStabilizingTransformation() --> the first is just a wrapper for the second. 
## 3. normTransform() --> log2 transformation or any other desired fuction. Will not be used!
## 4. lfcShrink
################################################################################
## Data Normalization for plotting: 1. VST Normalization
dds.vsd.TN <- varianceStabilizingTransformation(dds.run, blind=FALSE)
dds.vsd.TN
# the matrix of transformed values is stored in assay(dds.vsd.TN)

saveRDS(dds.vsd.TN, "saved_objects/dds.vsd.TN.rds")
dds.vsd.TN <- readRDS("saved_objects/dds.vsd.TN.rds")

# Extracting all DEGs
dds.vsd.TN.degs <- dds.vsd.TN[rownames(dds.vsd.TN) %in% rownames(res.degs),]
dds.vsd.TN.degs

# Extracting all repair and repair DEGs
dds.vsd.TN.rep <- dds.vsd.TN[rownames(dds.vsd.TN) %in% repair.genes$symbol,]
dds.vsd.TN.rep.degs <- dds.vsd.TN[row.names(dds.vsd.TN) %in% rownames(res.degs.rep)]
################################################################################
## Data Normalization for plotting: 2. lfcShrink Normalization
resultsNames(dds.run.TN)
res.lfc <- lfcShrink(dds.run.TN, coef=2, res=res, type = "apeglm")
summary(res.lfc)

## Chkeckin the distribution of the p-adjusted value
summary(res.lfc$padj)
hist(res.lfc$padj, main="Distribution of the adjusted p-value in the lfc-shrinked (apeglm) results", xlab="Adjusted p-value")
################################################################################
## Data Normalization for plotting: 3. normTransform
# dds.nrd <- normTransform(dds.run)
# dds.nrd
# saveRDS(dds.nrd, "saved_objects/dds.nrd.rds")
# dds.nrd <- readRDS("saved_objects/dds.nrd.rds")
################################################################################
## Plotting: 1. MAplot
## not-normalized data
to_tiff(plotMA(res, alpha=0.05, main="MA plot"), "MAplot.tiff")

## normalized data
# plotMA(res.lfc, alpha=0.05, main="MA plot of the lfcschrink-normalized DESeq Results")

# ## repair genes (normalized and not normalized)
# plotMA(res.rep, alpha=0.05, main="MA plot of the not-normalized DESeq Results for repair genes only")
# plotMA(res.lfc.rep, alpha=0.05, main="MA plot of the lfcschrink-normalized DESeq Results for repair genes only")
# 
# ## DEGs
# plotMA(res.degs, alpha=0.05, main="MA plot of the not-normalized all DEGs")
# plotMA(res.lfc.degs, alpha=0.05, main="MA plot of the lfc-shrunk all DEGs")
# 
# plotMA(res.degs.rep, alpha=0.05, main="MA plot of the not-normalized repair DEGs")
# plotMA(res.lfc.rep.degs, alpha=0.05, main="MA plot of the lfc-shrunk repair DEGs")
################################################################################
## Plotting: 2. EnhancedVolcano Plot
to_tiff(
  EnhancedVolcano(res, 
                  lab = rownames(res), 
                  x = 'log2FoldChange', 
                  y = 'pvalue',
                  title = NULL,
                  pCutoff = 0.05,
                  FCcutoff = 1.0,
                  legendPosition = "right", 
                  selectLab=rownames(res.degs.rep),
                  labCol = 'black'),
  "volcanoplot.tiff")
################################################################################
## Plotting: 3. PCA
# should we use normalized data instead? The manual says yes!
# what is the best ntop value? size of the DEGs?

## All genes
tiff(paste(figures_dir, "pca.plot.vsd.tiff", sep="/"),
     width = 18, height = 21, units = 'cm', res = 300)
pca.plot.vsd <- plotPCA(dds.vsd.TN, intgroup="Sample.Type", ntop = 25531, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd, "percentVar"))
ggplot(pca.plot.vsd, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of DESeq DataSet (ntop = 25531)")
rm(percentVar)
rm(pca.plot.vsd)
dev.off()

# all repair
pca.plot.vsd.rep <- plotPCA(dds.vsd.TN.rep, intgroup="Sample.Type", ntop = 285, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.rep, "percentVar"))
ggplot(pca.plot.vsd.rep, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of DESeq DataSet of the repair genes only (ntop = 285)")
rm(percentVar)
rm(pca.plot.vsd.rep)

## all DEGs
pca.plot.vsd.degs <- plotPCA(dds.vsd.TN.degs, intgroup="Sample.Type", ntop = 3828, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsd.degs, "percentVar"))
ggplot(pca.plot.vsd.degs, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of DEGs (ntop = 3828)")
rm(percentVar)
rm(pca.plot.vsd.degs)

## PCA for Repair DEGs
pca.plot.vsdrep.degs <- plotPCA(dds.vsd.TN.rep.degs, intgroup="Sample.Type", ntop = 34, returnData=TRUE)
percentVar <- round(100 * attr(pca.plot.vsdrep.degs, "percentVar"))
ggplot(pca.plot.vsdrep.degs, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  # xlab(percentage[1]) + ylab(percentage[2]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot of Repair DEGs (ntop = 34)")
rm(percentVar)
rm(pca.plot.vsdrep.degs)
################################################################################
## Plotting: 4. Dispersion estimate
# estimateDispersions(dds)
to_tiff(plotDispEsts(dds.run.TN), "plotDispEsts.tiff")
################################################################################
## Plotting: 5. meanSdPlot
# meanSdPlot()
################################################################################
## Plotting: 6. Counts plot
# genes.list <- c("PCLAF", "ISG15", "EXO1")
# genes.list <- rownames(res.degs.rep)
# for (gen in genes.list){
#   plotCounts(dds.run, gene=gen, intgroup = "Sample.Type") # returnData = FALSE
# }
# rm(genes.list)
# plotCounts(dds.run, gene=gen, intgroup = "Sample.Type")
################################################################################
## Plotting: 7. Boxplot with dots
list(assays(dds.run)) # counts mu H cooks replaceCounts replaceCooks
genes.list <- rownames(res.degs.rep)

par(mfrow=c(3,3))
for (gen in genes.list){
  plotting_gene(dds.run.TN, gen)
}
par(mfrow=c(1,1))
rm(genes.list)
# list(assays(dds.vsd.TN)) # counts mu H cooks replaceCounts replaceCooks
# boxplot(t(assays(dds.vsd.TN["EXO1"])[[1]])~dds.vsd.TN$Sample.Type, 
#         range=0, las=1, log='y', 
#         xlab="Groups", ylab="Counts", main="EXO1 (vsd-normalized)",
#         col=c("darksalmon", "darkred"))
# stripchart(t(assays(dds.vsd.TN["EXO1"])[[1]])~dds.vsd.TN$Sample.Type, 
#            vertical=TRUE, method='jitter', add=TRUE, pch=23, col="black", cex=0.5)
################################################################################
## Plotting: 8. Heatmap & Dendogram
colors2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)
colors.paired <- brewer.pal(12,"Paired")

cc <- c(IDC_Normal=colors.paired[3], IDC_Tumor=colors.paired[4], 
        LC_Normal=colors.paired[5], LC_Tumor=colors.paired[6],
        Mixed_Normal=colors.paired[7], Mixed_Tumor=colors.paired[8])[colData(dds.vsd.TN.rep.degs)$group]

# for repair degs
heatmap.2(assay(dds.vsd.TN.rep.degs), 
          col=colors2,
          labCol=substring(colnames(dds.vsd.TN.rep.degs), 6),
          scale="row", trace="none",
          # dendrogram = "row",
          
          Rowv=TRUE,
          Colv=order(dds.vsd.TN.rep.degs$Sample.Type),
          # Colv=order(dds.vsd.TN.rep.degs$group),
          # split = dds.vsd.TN.rep.degs$group,
          
          # ColSideColors = c(Primary.Tumor="darkgreen", Solid.Tissue.Normal="orange")[colData(dds.vsd.TN.rep.degs)$Sample.Type],
          ColSideColors = cc, colCol=cc,
          key =TRUE, key.title="Heatmap Key",
          main="Heatmap of the 34 Differentially Expressed Repair Genes")
################################################################################
################################################################################
## Including Interaction Term
design(ddsCollapsed) <- ~primary_diagnosis+Sample.Type+primary_diagnosis:Sample.Type
dds.run.ia <- DESeq(ddsCollapsed)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 3267 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing

saveRDS(dds.run.ia, "saved_objects/dds.run.ia.rds")
dds.run.ia <- readRDS("saved_objects/dds.run.ia.rds")

## Extracting Results
## https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html#example-4-two-conditionss-three-genotpes-with-interaction-terms
resultsNames(dds.run.ia)
levels(dds.run.ia$primary_diagnosis)

## Tumor vs. Normal for mixed subtype
res.ia.tn.mixed <- results(dds.run.ia, alpha=0.05,lfcThreshold=1, 
                           contrast = c("Sample.Type", "Tumor", "Normal"))
res.ia.tn.mixed <- res.ia.tn.mixed[order(res.ia.tn.mixed$padj),] # ordering
summary(res.ia.tn.mixed)
# out of 25,175 genes: 387 up-regulated and 620 down-regulated

## Tumor vs. Normal in IDC subtype
res.ia.tn.idc <- results(dds.run.ia, alpha=0.05,lfcThreshold=1, 
                         contrast=list("Sample.Type_Tumor_vs_Normal","primary_diagnosisIDC.Sample.TypeTumor"))
res.ia.tn.idc <- res.ia.tn.idc[order(res.ia.tn.idc$padj),] # ordering
summary(res.ia.tn.idc)
# out of 25,175 genes: 42 up-regulated and 116 down-regulated

## Tumor vs. Normal in LC subtype
res.ia.tn.lc <- results(dds.run.ia, alpha=0.05,lfcThreshold=1, 
                        contrast=list("Sample.Type_Tumor_vs_Normal","primary_diagnosisLC.Sample.TypeTumor"))
res.ia.tn.lc <- res.ia.tn.lc[order(res.ia.tn.lc$padj),] # ordering
summary(res.ia.tn.lc)
# out of 25,175 genes: 37 up-regulated and 159 down-regulated

## The effect of LC vs. IDC, with Normal.
res.ia.idc.lc.normal <- results(dds.run.ia, alpha=0.05,lfcThreshold=1, 
                         contrast = c(0,-1,1,0,0,0))
summary(res.ia.idc.lc.normal)
# out of 25,175 genes: 0 up-regulated and 5 down-regulated

## The interaction term for Tumor vs Normal in IDC vs genotype mixed
res.ia.idc.ia <- results(dds.run.ia, alpha=0.05,lfcThreshold=1, 
                         name="primary_diagnosisIDC.Sample.TypeTumor")
summary(res.ia.idc.ia)
# out of 25,175 genes: 2 up-regulated and 1 down-regulated

## The interaction term for Tumor vs Normal in LC vs genotype mixed
res.ia.lc.ia <- results(dds.run.ia, alpha=0.05,lfcThreshold=1, 
                        name="primary_diagnosisLC.Sample.TypeTumor")
summary(res.ia.lc.ia)
# out of 25,175 genes: 5 up-regulated and 1 down-regulated

## The interaction term for tumor vs normal in IDC vs LC.
res.ia.idc.lc.tumor <- results(dds.run.ia, alpha=0.05,lfcThreshold=1, 
                               contrast=list("primary_diagnosisIDC.Sample.TypeTumor", "primary_diagnosisLC.Sample.TypeTumor"))
summary(res.ia.idc.lc.tumor)
# out of 25,175 genes: 0 up-regulated and 2 down-regulated
################################################################################
################################################################################
### Saving all figures from the plots tab at once 
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="saved_objects/figures/")
################################################################################