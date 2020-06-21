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

install.packages("stringr")       # packageVersion: 1.4.0
install.packages("gplots")        # packageVersion: 3.0.3
install.packages("ggplot2")       # packageVersion: 3.3.1
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

library(DESeq2)
library(EnhancedVolcano)
# library(genefilter)

library(gplots)
library(ggplot2)
library(RColorBrewer)
# library(cluster)
# library(biclust)
# library(ggdendro)
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
sample_sheet <- read.csv("raw_data/females_only.ductal-and-lobular-only.gdc_sample_sheet.2020-06-17.tsv", sep="\t", header=TRUE)
table(sample_sheet$Sample.Type)

## dropping metastatic samples from the sample sheet
sample_sheet <- sample_sheet[! sample_sheet$Sample.Type == "Metastatic",]
table(sample_sheet$Sample.Type)

## dropping metastatic samples from the RNA-seq expression dataframe
rna.exp.df.subsetted <- rna.exp.df[, colnames(rna.exp.df) %in% sample_sheet$File.ID]

saveRDS(rna.exp.df.subsetted, file = "saved_objects/exp.rna.subsetted.rds")
rna.exp.df.subsetted <- readRDS("saved_objects/exp.rna.subsetted.rds")
################################################################################
## Changing ENSG to SYMBOL
rna.exp.df.sym <- ENSG_to_symbol(rna.exp.df.subsetted)
saveRDS(rna.exp.df.sym, file = "saved_objects/rna.exp.df.sym.rds")
rna.exp.df.sym <- readRDS(file = "saved_objects/rna.exp.df.sym.rds")
## 60,483 ENSG where reduced to 25,596 gene symbol (org.Hs.eg.db version 3.11.4)

rna.exp.df.sym.agg <- aggregate_rows(rna.exp.df.sym, agg.var=rna.exp.df.sym$symbol) # This step takes too long to run!
saveRDS(rna.exp.df.sym.agg, file = "saved_objects/rna.exp.df.sym.agg.rds")
rna.exp.df.sym.agg <- readRDS(file="saved_objects/rna.exp.df.sym.agg.rds")
rm(rna.exp.df.sym)
################################################################################
## Preparing the inputs of DESEQ2
count.data <- rna.exp.df.sym.agg

col.data <- sample_sheet[, colnames(sample_sheet) %in% c("File.ID", "Sample.Type", "Sample.ID")]
col.data$Sample.Type <- as.factor(col.data$Sample.Type)
col.data$Sample.Type <- relevel(col.data$Sample.Type, "Solid Tissue Normal")
levels(col.data$Sample.Type) <- c("Solid.Tissue.Normal", "Primary.Tumor")
table(col.data$Sample.Type)

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
dds <- DESeqDataSetFromMatrix(countData=count.data , colData=col.data , design=~Sample.Type)
dds
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
################################################################################
## Run DESEQ2 
dds.run <- DESeq(ddsCollapsed)     ## This step takes almost 4 hours to run
saveRDS(dds.run, file="saved_objects/dds.run.rds")
dds.run <- readRDS("saved_objects/dds.run.rds")
################################################################################
## Creating the results (RES) object
res <- results(dds.run, contrast=c("Sample.Type", "Primary.Tumor", "Solid.Tissue.Normal"), alpha=0.05,lfcThreshold=1)
summary(res)

## Removing incomplete genes
res <- res[complete.cases(res),]
summary(res)
##  genes,  up-regulated,  down-regulated

## saving to and reading from RDS object
saveRDS(res, file="saved_objects/res.rds")
res <- readRDS("saved_objects/res.rds")

## Converting the results to DF
res.df <- as.data.frame(res)
res.df <- res.df[order(res.df$padj),]
write.csv(res.df, "saved_objects/all.results.sorted.csv")

## Extracting DEGs
res.degs <- res.df[res.df$padj<0.05 & abs(res.df$log2FoldChange)>log2(2),]
write.csv(res.degs, "saved_objects/degs-p.05-FC1.csv")
################################################################################
## Repair genes subsetting (only one of the two following sets!)

dds.run.rep <- dds.run.rep <- dds.run[rownames(dds.run) %in% repair.genes$symbol,]
res.rep <- results(dds.run.rep, contrast=c("Sample.Type", "Primary.Tumor", "Solid.Tissue.Normal"), alpha=0.05,lfcThreshold=1)
summary(res.rep)
res.rep <- res.rep[complete.cases(res.rep),]
summary(res.rep)

res.rep <- res[, colnames(res) %in% repair.genes$symbol] # not sure of the code!
summary(res.rep)
################################################################################
## Normalization Methods Available:
## 1. rlog()  --> takes too long and not applicable for our dataset (not as sensetive as vst and takes too long and fails)
## 2. vst() or varianceStabilizingTransformation() --> the first is just a wrapper for the second. 
## 3. normTransform()
################################################################################
## Data Normalization for plotting: 1. VST Normalization
dds.vsd <- varianceStabilizingTransformation(dds.run)
dds.vsd.res <- results(dds.vsd, contrast=c("Sample.Type", "Primary.Tumor", "Solid.Tissue.Normal"), alpha=0.05,lfcThreshold=1)
summary(dds.vsd.res)
# the matrix of transformed values is stored in assay(vsd)
################################################################################
## Data Normalization for plotting: 2. lfcShrink Normalization
resultsNames(dds.run)
res.lfc <- lfcShrink(dds.run, coef=2, res=res, type = "apeglm") # change the coef with resultsNames(dds.run)
# should we specify lfcThreshold and get s-value instead of p-value or not?
################################################################################
## Data Normalization for plotting: 1. normTransform
nrd <- normTransform(dds.run)
################################################################################
## Plotting: 1. MAplot
plotMA(res, main="MA plot of the DESeq Results (LFC=1.0, p-val=0.05)")
plotMA(dds.run, main="MA plot of the DESeq DataSet")

# should we use normalized data instead?
plotMA(res.lfc, main="MA plot of the lfcschrink-normalized DESeq Results (LFC=1.0, p-val=0.05)")
plotMA(dds.vsd, main="MA plot of the vsd-normalized DESeq DataSet")
plotMA(dds.vsd.res, main="MA plot of the vsd-normalized DESeq Results (LFC=1.0, p-val=0.05)")

# what about the repair genes?
################################################################################
## Plotting: 2. EnhancedVolcano Plot
EnhancedVolcano(res, 
                lab = rownames(res), 
                x = 'log2FoldChange', 
                y = 'pvalue',
                title = "Volcano plot of the results (LFC=1.0, p-val=0.05)",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                legendPosition = "right")

# should we use normalized data instead?
EnhancedVolcano(res.lfc, 
                lab = rownames(res.lfc), 
                x = 'log2FoldChange', 
                y = 'pvalue',
                title = "Volcano plot of the lfcschrink-normalized DESeq Results (LFC=1.0, p-val=0.05)",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                legendPosition = "right")
EnhancedVolcano(dds.vsd.res, 
                lab = rownames(dds.vsd.res), 
                x = 'log2FoldChange', 
                y = 'pvalue',
                title = "Volcano plot of the vsd-normalized DESeq Results (LFC=1.0, p-val=0.05)",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                legendPosition = "right")


# what about the repair genes?
EnhancedVolcano(res.rep, 
                lab = rownames(res.rep), 
                x = 'log2FoldChange', 
                y = 'pvalue',
                title = "Volcano plot of the repair genes results (LFC=1.0, p-val=0.05)",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                legendPosition = "right")
################################################################################
## Plotting: 3. PCA
# should we use normalized data instead? The manual says yes!
# what is the best ntop value? size ofthe DEGs?
pca.plot <- plotPCA(dds.vsd, intgroup="Sample.Type", ntop = 25531, returnData=TRUE, )
ggplot(pca.plot, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(percentage[1]) + ylab(percentage[2]) +
  ggtitle("PCA plot of the vsd-normalized DESeq DataSet")

# if not normalized 
pca.plot <- plotPCA(dds.run, intgroup="Sample.Type", ntop = 25531, returnData=TRUE, )
ggplot(pca.plot, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(percentage[1]) + ylab(percentage[2]) + 
  ggtitle("PCA plot of the not-normalized DESeq DataSet")

# what about the repair genes?
pca.plot <- plotPCA(dds.run.rep, intgroup="Sample.Type", ntop = 310, returnData=TRUE, )
ggplot(pca.plot, aes(PC1, PC2, color=Sample.Type)) + 
  geom_point(size=1) + stat_ellipse(type = "norm") + 
  xlab(percentage[1]) + ylab(percentage[2]) + 
  ggtitle("PCA plot of the not-normalized DESeq DataSet of the repair genes only")
################################################################################
## Plotting: 4. Dispersion estimate
# estimateDispersions(dds)
plotDispEsts(dds) # dds not dds.run!
plotDispEsts(dds.run) # just out of curiosity!
################################################################################
## Plotting: 5. meanSdPlot
################################################################################
################################################################################
### Saving all figures from the plots tab at once 
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="saved_objects/figures/")
################################################################################
