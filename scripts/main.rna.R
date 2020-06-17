################################################################################
################################ Session Info. #################################
## R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
## Copyright (C) 2020 The R Foundation for Statistical Computing
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 18.04.4 LTS
## BiocManager version: 1.30.10

## Created by Ruwaa I. Mohamed
## Date: Saturday, May 30, 2020.
sessionInfo()

################################################################################
################################################################################
## Second Device 
## R version 4.0.1 (2020-06-06) -- "See Things Now"
## Copyright (C) 2020 The R Foundation for Statistical Computing
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04 LTS

## BiocManager version: 1.30.10
## Bioconductor version: 3.11
################################################################################
## Installing Bioconductor and the required libraries
# use packageVersion('') to check the versio of all installed packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
## BiocManager::install(version = "1.30.10") # <-- BiocManager version
BiocManager::version() # 3.9 <-- bioconductor version

BiocManager::install("org.Hs.eg.db")        # packageVersion: 3.8.2
BiocManager::install('EnhancedVolcano')     # packageVersion: 1.2.0

install.packages("gplots")        # packageVersion: 3.0.3
install.packages("ggplot2")       # packageVersion: 3.3.1
# install.packages("cluster")     # packageVersion: 2.1.0
install.packages("RColorBrewer")  # packageVersion: 1.1.2
install.packages("genefilter")    # packageVersion: 1.66.0
# install.packages("biclust")     # packageVersion: 2.0.2
# install.packages("ggdendro")    # packageVersion: 0.1.20
################################################################################
## Loading the required libraries/scripts
#scripts
source('scripts/get_df_from_gdc.R')
source('scripts/ENSG_to_symbol.R')
source('scripts/aggregate_rows.R')

# libraries
library(DESeq2)
library(gplots)
library(ggplot2)
library(EnhancedVolcano)
# library(cluster)
library(RColorBrewer)
library(genefilter)
library(biclust)
# library(ggdendro)
################################################################################
## Reading the RNA-seq files (1,222 files)
rna.data.path <- "./raw_data/RNA-seq/"

## reading the expression values from the raw data and saving them to RDS object (Don't use this is unless you don't have the RDS object)!
rna.exp.df <- get_df_from_gdc(rna.data.path, 'htseq.counts.gz$')
rna.exp.df <- rna.exp.df[1:(nrow(rna.exp.df)-5),]
saveRDS(rna.exp.df, file = "./saved_objects/exp.rna.rds")

## reading the expression values from the previously saved RDS object (use this instead)!
rna.exp.df <- readRDS(file = "./saved_objects/exp.rna.rds")
View(rna.exp.df) # DON'T!
################################################################################
## Reading the sample sheet
sample_sheet <- read.table('./raw_data/gdc_sample_sheet.2020-05-29.tsv', sep='\t', header=TRUE)
## subsetting only RNA-seq related records.
rna.sample_sheet <- sample_sheet[sample_sheet$File.ID %in% colnames(rna.exp.df),]
table(rna.sample_sheet$Sample.Type)
mirna.sample_sheet <- sample_sheet[! sample_sheet$File.ID %in% colnames(rna.exp.df),]
table(mirna.sample_sheet$Sample.Type)
write.csv(mirna.sample_sheet, "saved_objects/mirna.sample_sheet.csv")
################################################################################
## Changing ENSG to SYMBOL
rna.exp.df.sym <- ENSG_to_symbol(rna.exp.df)
saveRDS(rna.exp.df.sym, file = "./saved_objects/rna.exp.df.sym.rds")
rna.exp.df.sym <- readRDS(file = "./saved_objects/rna.exp.df.sym.rds")
## 60,483 ENSG where reduced to 25,393 gene symbol

rna.exp.df.sym.agg <- aggregate_rows(rna.exp.df.sym, agg.var=rna.exp.df.sym$symbol) # This step takes too long to run!
saveRDS(rna.exp.df.sym.agg, file = "./saved_objects/rna.exp.df.sym.agg.rds")
rna.exp.df.sym.agg <- readRDS(file="./saved_objects/rna.exp.df.sym.agg.rds")
rm(rna.exp.df.sym)
## After aggregation, we have 25,330 gene symbol after aggregation
################################################################################
## Dividing the df into normal, tumor, and metastatic
normal.FileIDs <- rna.sample_sheet[rna.sample_sheet$Sample.Type == 'Solid Tissue Normal',]$File.ID
normal.rna.exp <- rna.exp.df.sym.agg[, names(rna.exp.df.sym.agg) %in% normal.FileIDs]
rm(normal.FileIDs)
saveRDS(normal.rna.exp, file='./saved_objects/normal.rna.exp.all.rds')
normal.rna.exp <- readRDS(file='./saved_objects/normal.rna.exp.all.rds')

tumor.FileIDs <- rna.sample_sheet[rna.sample_sheet$Sample.Type == 'Primary Tumor',]$File.ID
tumor.rna.exp <- rna.exp.df.sym.agg[, names(rna.exp.df.sym.agg) %in% tumor.FileIDs]
rm(tumor.FileIDs)
saveRDS(tumor.rna.exp, file='./saved_objects/tumor.rna.exp.all.rds')
tumor.rna.exp <- readRDS(file='./saved_objects/tumor.rna.exp.all.rds')

metastatic.FileIDs <- rna.sample_sheet[rna.sample_sheet$Sample.Type == 'Metastatic',]$File.ID
metastatic.rna.exp <- rna.exp.df.sym.agg[, names(rna.exp.df.sym.agg) %in% metastatic.FileIDs]
rm(metastatic.FileIDs)
saveRDS(metastatic.rna.exp, file='./saved_objects/metastatic.rna.exp.all.rds')
metastatic.rna.exp <- readRDS(file='./saved_objects/metastatic.rna.exp.all.rds')

################################################################################
## Preparing the inputs of DESEQ2
count.data <- cbind(normal.rna.exp, tumor.rna.exp)
col.data <- rna.sample_sheet[, colnames(rna.sample_sheet) %in% c("File.ID", "Sample.Type", "Sample.ID")]
col.data <- col.data[! col.data$Sample.Type == "Metastatic",]
col.data <- droplevels(col.data)
col.data$Sample.Type <- relevel(col.data$Sample.Type, "Solid Tissue Normal")
levels(col.data$Sample.Type) <- c("Solid.Tissue.Normal", "Primary.Tumor")
table(col.data$Sample.Type)

new.order <- match(colnames(count.data), col.data$File.ID)
count.data.sorted <- count.data[order(new.order)]
rm(new.order)

count.data.sorted <- apply (count.data.sorted, 2, as.integer)
rownames(count.data.sorted) <- rownames(normal.rna.exp)

all(colnames(count.data.sorted) == col.data$File.ID)
################################################################################
## Run DESEQ2
## This tutorial was used mainly >> https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf
## Prepare the DDS object
dds <- DESeqDataSetFromMatrix(countData=count.data.sorted , colData=col.data , design=~Sample.Type)
dds

## Collapse Technocal Replicates
ddsCollapsed <- collapseReplicates(dds, groupby=dds$Sample.ID)
ddsCollapsed

## to check that collapsing worked!
original <- rowSums(counts(dds)[, dds$Sample.ID == "TCGA-A7-A26J-01A"])
all(original == counts(ddsCollapsed)[,"TCGA-A7-A26J-01A"])

## Run DESEQ2 
dds.run <- DESeq(ddsCollapsed)     ## This step takes almost 4 hours to run
saveRDS(dds.run, file="saved_objects/dds.run.rds")
dds.run <- readRDS("saved_objects/dds.run.rds")
################################################################################
### Create res object based on cases vs. controls (normal vs. tumor)
# res <- results(dds.run)
res.default <- results(dds.run, contrast=c("Sample.Type", "Primary.Tumor", "Solid.Tissue.Normal"))
summary(res.default)

res <- results(dds.run, contrast=c("Sample.Type", "Primary.Tumor", "Solid.Tissue.Normal"), alpha=0.05, lfcThreshold=1)
summary(res)
## no. of genes: 24,981

## Removing NULLs from the res object
res.default <- res.default[complete.cases(res.default),]
summary(res.default)

res <- res[complete.cases(res),]
summary(res)
# Nothing changed for the res obkect in this case!
## no. of genes: 24,981

saveRDS(res, file="saved_objects/res.rds")
res <- readRDS("saved_objects/res.rds")

## Converting the results to DF
res.df <- as.data.frame(res)
res.df.sorted <- res.df[order(res.df$padj),]
write.csv(res.df.sorted, "saved_objects/all.results.sorted.csv")
# res.ups <- res.df[res.df$log2FoldChange > 0,]
# res.dwns <- res.df[res.df$log2FoldChange < 0,]

## Select p-value and FC
res.degs <- res.df.sorted[res.df.sorted$padj<0.05 & abs(res.df.sorted$log2FoldChange)>log2(2),]
write.csv(res.degs, "saved_objects/degs-p.05-FC1.csv")
## p-value 0.05 & log-FC 1 --> 3,877 genes
# dim(res.degs[res.degs$log2FoldChange > 0,])
# dim(res.degs[res.degs$log2FoldChange < 0,])
## (up-regulated: 2,273 and downregulated: 1,604)
## This is matching summary(res)
################################################################################
## Comparing me and Salma
ruwaa <- res.degs
salma <- read.csv("saved_objects/res.sym.degs.salma.may30.csv", header=TRUE)
salma <- salma[,-1]
salma.uniq <- salma[! salma$symbol %in% rownames(ruwaa),]$symbol
salma.uniq <- droplevels(salma.uniq)
ruwaa.uniq <- rownames(ruwaa[! rownames(ruwaa) %in% salma$symbol,])
common <- rownames(ruwaa[rownames(ruwaa) %in% salma$symbol,])
################################################################################
## MA Plot
plotMA(res, main="MA plot of res with LFC=1, p-val=0.05, filtered on res (LFC=1, pvalue=0.05)")
plotMA(dds.run, main="MA plot of dds.run")
# plotMA(res.default, main="MA plot of res with default vales (LFC=0, p-val=0.1), alpha here is 0.05", alpha=0.05)     ## Please consider that the p-value is 0.1 here and the LFC is 0.
################################################################################
## Volcano plot
# volcanoPlot(res)
vp <- EnhancedVolcano(res, 
                      lab = rownames(res), 
                      x = 'log2FoldChange', 
                      y = 'pvalue',
                      title = "Volcano plot of the results (p-val=0.05, LFC=1.0)",
                      pCutoff = 0.05,
                      FCcutoff = 1.0,
                      legendPosition = "right")
vp
################################################################################
## Histogram of adjusted p-values distribution
pval.his <- hist(res$padj, main="Distribution of the adjusted value, res with LFC=1, p-val=0.05")
# pval.his_res.default <- hist(res.default$padj, main="Distribution of the adjusted value, default res object (LFC=0, p-val=0.1)")
################################################################################
## Count plot of top DEGs

################################################################################
### Heatmap Sample clustering 

## regularized logarithm transformation
# rld <- rlog(dds.run)     # R session aborted!
# head(assay(rld))

## variance stabilizing transformations (VST)
vsd <- vst(dds.run)      # This step takes too long to run. Don't run this again unless you re-run the DeSEQ2
vsd <- varianceStabilizingTransformation(dds.run)
saveRDS(vsd, file = "./saved_objects/vsd.rds")
vsd <- readRDS(file="./saved_objects/vsd.rds")

#  Sample distances
sampleDists <- dist(t(assay(vsd)))     # either rld or vsd
saveRDS(sampleDists, file = "./saved_objects/sampleDists.rds")
sampleDists <- readRDS(file="./saved_objects/sampleDists.rds")

# plot into heat map
sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(vsd$Sample.Type, sd$Sample.ID, sep="-" )
# colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
hm.sam <- heatmap.2(sampleDistMatrix, trace="none", col=colors, main="Heatmap of the Samples", dendrogram="row")
hm.sam
################################################################################
# Dendogram
dendo <- hm.sam$rowDendrogram
plot(dendo, main="Dendogram of the Samples from the Heatmap")
################################################################################
# plot Principal components analysis (PCA) 
pca.plot <- plotPCA(vsd, intgroup="Sample.Type", ntop = 25330, returnData=TRUE, )
ggplot(pca.plot, aes(PC1, PC2, color=Sample.Type)) + geom_point(size=1) + stat_ellipse(type = "norm")
# pca.plot$frame <- TRUE
# pca.plot$frame.type <- 'norm'
# pca.plot$dotSize <- 2
# ggplot_build(pca.plot)
# plot(pca.plot)
################################################################################
##Heatmap with gene clustering
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing=TRUE ), 50)
colors2 <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
hm.gen <- heatmap.2(assay(vsd)[topVarGenes,], scale="row", trace="none", col=colors2, key.title="Heatmap of the genes")
hm.gen
################################################################################
## Biclustering Trial
# bics <- biclust(as.matrix(dist(assay(vsd))), BCSpectral)
# drawHeatmap(as.matrix(assay(vsd)[topVarGenes,]), bicResult=bics)

################################################################################
################################################################################
################################################################################
## subsetting all repair-related genes.
repair.genes <- read.csv("saved_objects/repair.genes.csv", row.names = "X")
# res.repair.genes.df <- res.df.sorted[rownames(res.df.sorted) %in% repair.genes$x,]

dds.run.rep <- dds.run[rownames(dds.run) %in% repair.genes$x,]
dds.run.rep

res.rep <- results(dds.run.rep, contrast=c("Sample.Type", "Primary.Tumor", "Solid.Tissue.Normal"), alpha=0.05, lfcThreshold=1)
summary(res.rep)

res.rep <- res.rep[complete.cases(res.rep),]
summary(res.rep)

res.rep.sorted <- res.rep[order(res.rep$padj),]
res.rep.sorted.df <- as.data.frame(res.rep.sorted)
write.csv(res.rep.sorted.df, "saved_objects/repair_genes_results.csv")

res.rep.degs <- res.rep.sorted.df[res.rep.sorted.df$padj < 0.05 & abs(res.rep.sorted.df$log2FoldChange) > 1,]
write.csv(res.rep.degs, "saved_objects/repair_genes_results_degs-only.csv")

saveRDS(res.rep, file="saved_objects/res.rep.rds")
res.rep <- readRDS("saved_objects/res.rep.rds")
################################################################################
################################################################################
################################################################################
### Saving all figures from the plots tab at once 
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="saved_objects/figures/")
################################################################################
## 1. Data + Experiment design (male and female)  Done
## 2. figures: Dispersion estimate and counts (and their data normalization)
## 3. Heatmap
