# TCGA-BRCA

### Aim
The aim of this project is to analyze the change of gene expression levels in lobular and ductal neoplasms breast cancer subtypes, with focus on personalized approach for all genes involved in repair pathways. 

### Methodology:

#### 1. Selecting the dataset (Cases)
The dataset selected was from the [TCGA-BRCA Project](https://portal.gdc.cancer.gov/projects/TCGA-BRCA "TCGA-BRCA Project Page"). The study was made on 1,098 cases: 1,085 females and 12 males. After excluding the males, 1,041 of the females are with ductal and lobular neoplasms. Of them, only **1,036 cases** had RNA-seq Gene Expression Quantification available. 

#### 2. Description of the Samples
The **1,036** female cases with ductal and lobular neoplasms breast cancer contributed with **1,164** RNA-seq, HTSeq Counts (raw counts) samples/files (1,046 samples from primary tumor, 111 samples from solid tissue normal, and 7 metastatic samples). 
The filters applied in the TCGA data repository to get those files were: 

1. Data Category: Transcriptome Profiling
2. Data Type: Gene Expression Quantification
3. Experimental Strategy: RNA-Seq
4. Workflow Type: HTSeq - Counts.

miRNA data of the same cases were downloaded. The filters applied in the TCGA data repository to get those files were: 

1. Data Category: Transcriptome Profiling
2. Data Type: miRNA Expression Quantification
3. Experimental Strategy: miRNA-Seq
4. Workflow Type: BCGSC miRNA Profiling

miRNA-seq data available were taken from only **1,022** females not 1,036 as in the case of the RNA-seq data. They, also, contributed with **1,164** BCGSC miRNA Profiling files/samples (1,037 samples from primary tumor, 102 samples from solid tissue normal, and 7 metastatic samples).

#### 3. Downloading the dataset
On March 05, 2020, the files selected were added to the cart on GDC portal. The sample sheet, manifest file, clinical data, and metadata files were downloaded. After that, the counts' files were downloaded via The GDC Data Transfer Tool ([gdc-client](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Getting_Started/ "The GDC Data Transfer Tool")) on the same day using the following `Bash` commands.

```bash
# create a working environment
mkdir ~/TCGA-BRCA
cd ~/TCGA-BRCA
mkdir raw_data scripts saved_objects
mkdir raw_data/RNA-seq raw_data/miRNA-seq
# move the files downloaded from TCGA GDC data portal to ~/TCGA-BRCA/raw_data

# Downloading the gdc-client tool.
wget https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.5.0_Ubuntu_x64.zip 
unzip gdc-client_v1.5.0_Ubuntu_x64.zip
alias gdc-client="./gdc-client"
gdc-client download -h

# Downloading the RNA-seq Data
gdc-client download -m raw_data/rna-seq_gdc_manifest.2020-03-05.txt -d raw_data/RNA-seq/

# Downloading the miRNA-seq Data
gdc-client download -m raw_data/mirna-seq_gdc_manifest.2020-03-05.txt -d raw_data/miRNA-seq/
```

#### 4. Genes in Repair Pathways
A list of all genes involved in the following 12 pathways were downloaded from [Reactome](https://reactome.org/ "Reactome database") database.
1. BER Participating Molecules [[R-HSA-73884](https://reactome.org/content/detail/R-HSA-73884)]
2. DNA BypassParticipating Molecules [[R-HSA-73893](https://reactome.org/content/detail/R-HSA-73893)]
3. DNA damage Reversal Participating Molecules [[R-HSA-73942](https://reactome.org/content/detail/R-HSA-73942)]
4. DNA double strand break response Participating Molecules [[R-HSA-5693606](https://reactome.org/content/detail/R-HSA-5693606)]
5. Fanconi Anemia Participating Molecules [[R-HSA-6783310](https://reactome.org/content/detail/R-HSA-6783310)]
6. HDR through HRR alone Participating Molecules [[R-HSA-5693567](https://reactome.org/content/detail/R-HSA-5693567)]
7. HDR through HRR and SSA Participating Molecules [[R-HSA-5693567](https://reactome.org/content/detail/R-HSA-5693567)]
8. HDR through MMEJ Participating Molecules [[R-HSA-5685939](https://reactome.org/content/detail/R-HSA-5685939)]
9. HDR through SSA Participating Molecules [[R-HSA-5693567](https://reactome.org/content/detail/R-HSA-5693567)]
10. MMR Participating Molecules [[R-HSA-5358508](https://reactome.org/content/detail/R-HSA-5358508)]
11. NER Participating Molecules [[R-HSA-5696398](https://reactome.org/content/detail/R-HSA-5696398)]
12. Non homologus end joining Participating Molecules [[R-HSA-5693571](https://reactome.org/content/detail/R-HSA-5693571)]

Because some genes are involved in multiple pathways, the R-script `proteins.R` was used to get the unique gene symbol from the **12 pathways: 310 genes**; they were saved in `repair.genes.csv`

#### 5. Analysis Pipeline
1. All RNA-seq expression files were read into a single dataframe and metastatic samples were dropped (60,483 ENSG transcripts x 1,157 samples). The samples were 1,046 from Primary Tumor and 111 from Solid Tissue Normal >> saved as `exp.rna.rds` object.

2. ENSG IDs were converted to gene symbol using `org.Hs.eg.db (3.11.4)` R package, then duplicate gene symbol entries were aggregated by averaging the values of the duplicate rows (25,531 genes x 1,157 samples) >> saved as `exp.rna.sym.agg.rds` object.

3. Standard `DESeq2 (1.28.1)` analysis was run as illustrated in the common package vignettes. 
   1. A DESeqDataSet object was created from the expression data and the sample sheet. The order of the samples in the colnames of the expression matrix was matched with the order of the samples column in the sample sheet.
   2. Technical replicates (5 samples) were collapsed.
   3. DESeq2 was run on the DESeqDataSet object created (25,531 genes x 1,152 samples) with study design of comparing the primary tumor to the solid tissue normal samples (it performs all needed data transformation).
   4. A DESeqResults object was retried using `results()` function. The `DESeqResults` object had **25,177** genes. Of them, only **285** out of the 310 repair genes were found. The p-value and LFC of the study were specified at this step (lfcThreshold=1; p-value=0.05).
   
4. Significant DEGs were collected with LFC > 1 or less than < -1 and with adjusted p-value < 0.05:  **2,233** genes are up-regulated and **1,639** genes down-regulated. Only **35** genes out of the 285 repair genes were significant on these conditions; all are up regulated.

5. Data were normalized using `varianceStabilizingTransformation()` on the `DESeqDataSet` object level or `lfcShrink()` on the `DESeqResults` object for data visualization.

7. Data visualization and plotting:
   
   1. Histogram of distribution of p-value and LFC in the results:
   
      Two histogram figures of the distribution of the LFC of all results and of repair genes only were plotted as in the two figures below.
   
      ---------- | -----------
   
      ![LFC Distribution of the results](saved_objects/figures/6982ec97-2155-468f-8a2d-af27f8a7e03a.png) |  ![LFC Distribution of the results (repair genes only)](saved_objects/figures/9e4ad53c-aa58-4949-aa24-b9e7c27e7028.png)
   
      Another two  histogram figures of the distribution of the adjusted p-value of all results and of repair genes only were plotted as in the two figures below.
   
      
   
   2. MA/Volcano plot:
   
      
   
      1. MA plot
      2. Volcano plot
   
   3. PCA plot
   
   4. Dispersion Estimates plot
   
   5. Counts/Box plot
   
      1. Count plot
      2. Boxplot
   
   6. Heatmap and genes clustering
