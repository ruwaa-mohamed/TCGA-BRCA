# TCGA-BRCA

### Aim
The aim of this project is to analyze the change of gene expression levels in lobular and ductal neoplasms breast cancer subtypes, with focus on personalized approach for all genes involved in repair pathways. 

### Methodology:

#### Selecting the dataset
The dataset selected was from the [TCGA-BRCA Project](https://portal.gdc.cancer.gov/projects/TCGA-BRCA "TCGA-BRCA Project Page"). The study was made on 1,098 cases: 1,085 emales and 12 males. After excluding the males, 1,041 of the females are with ductal and lobular neoplasms. Of them, only 1,036 cases had RNA-seq Gene Expression Quantification available. 


#### Downloading the dataset
The dataset **TCGA-BRCA** was downloaded from TCGA database as manifest files on March 05, 2020.
After that, GDC-data transfer tool was used to download the data using the following `Bash` commands:

```bash
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

#### Genes in Repair Pathways

#### Analysis Pipeline
1. Read the raw data of all RNA-seq expression files into a single dataframe >> save it as `exp.rna.rds` object.
2. Collect all protein UniProt IDs of the DNA repair pathways of interest in one array >> save it as `proteins.ids.rds` object.
3. convert the transcript IDs to uniprot IDs and aggregate dublicate by averaging.
4. Select the rows of the proteins of interest only >> save it as `exp.rna.repair.agg.rds` object.
5. convert file IDs to Sample IDs in the working dataframe based on the sample datasheet.
6. subsetting the normal, tumor, and metastatic samples to different dfs based on the sample datasheet.
7. using DESeq2 for DGE analysis.
