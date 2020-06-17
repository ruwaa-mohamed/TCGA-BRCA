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
################################################################################
## Loading the required libraries/scripts
source('scripts/get_df_from_gdc_files.R')
################################################################################
mirna.sample_sheet <- read.csv("saved_objects/mirna.sample_sheet.csv", header=TRUE)
table(mirna.sample_sheet$Sample.Type)
mirna.exp.df <- get_df_from_gdc_mirna(data.path="raw_data/miRNA-seq/", file_ext=".mirnas.quantification.txt")
