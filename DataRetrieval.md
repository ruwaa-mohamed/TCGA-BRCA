## Data Retrieval

### 1. Selecting the dataset (Cases)

The dataset selected was from the [TCGA-BRCA Project](https://portal.gdc.cancer.gov/projects/TCGA-BRCA "TCGA-BRCA Project Page"). The study was made on 1,098 cases: 1,085 females and 12 males. After excluding the males, 1,041 of the females are with ductal and lobular neoplasms. Of them, only **1,036 cases** had RNA-seq Gene Expression Quantification available, of them, only 1,018 cases had both. 

### 2. Description of the Samples

The **1,036** female cases with ductal and lobular neoplasms breast cancer contributed with **1,164** RNA-seq, HTSeq Counts (raw counts) samples/files (1,046 samples from primary tumor, 111 samples from solid tissue normal, and 7 metastatic samples). 
The filters applied in the GDC data repository to get those files were: 

1. Data Category: Transcriptome Profiling
2. Data Type: Gene Expression Quantification
3. Experimental Strategy: RNA-Seq
4. Workflow Type: HTSeq - Counts.

miRNA data of the same cases filter was used in this study. Out of the 1,041 females with ductal and lobular neoplasms, only **1,022** females had miRNA-seq files (not 1,036 as in the case of the RNA-seq data). However, 4 of them were not included in the cases with RNA-seq data available and were excluded later in our analysis. The filters applied in the GDC data repository to get the miRNA-seq files were: 

1. Data Category: Transcriptome Profiling
2. Data Type: miRNA Expression Quantification
3. Experimental Strategy: miRNA-Seq
4. Workflow Type: BCGSC miRNA Profiling

The 1,022 cases contributed with **1,164** BCGSC miRNA Profiling files/samples (1,037 samples from primary tumor, 102 samples from solid tissue normal, and 7 metastatic samples). After removing the metastatic entries and the 4 cases not included in the RNA-seq data, we started our analysis with **1,018 cases** with **1,135 files** (1,033 primary tumor and 102 solid tissue normal samples).

### 3. Downloading the dataset

On March 05, 2020, the files selected were added to the cart on GDC portal. The sample sheet, manifest file, clinical data, and metadata files were downloaded from the GDC data portal. After that, the counts' files were downloaded via The GDC Data Transfer Tool ([gdc-client](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Getting_Started/ "The GDC Data Transfer Tool")) on the same day using the following `Bash` commands.

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

### 4. Genes in Repair Pathways

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

Because some genes are involved in multiple pathways, the [Proteins](./scripts/proteins.R) R-script was used to get the unique gene symbol from the **12 pathways: 310 genes**; they were saved in [Repair Genes](./saved_objects/repair.genes.csv) CSV file.
