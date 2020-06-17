# TCGA-BRCA

### Aim
The aim of this project is to analyze the change of gene expression levels in lobular and ductal neoplasms breast cancer subtypes, with focus on personalized approach for all genes involved in repair pathways. 

### Methodology:

#### 1. Selecting the dataset (Cases)
The dataset selected was from the [TCGA-BRCA Project](https://portal.gdc.cancer.gov/projects/TCGA-BRCA "TCGA-BRCA Project Page"). The study was made on 1,098 cases: 1,085 females and 12 males. After excluding the males, 1,041 of the females are with ductal and lobular neoplasms. Of them, only 1,036 cases had RNA-seq Gene Expression Quantification available. 

#### 2. Description of the Samples
The 1,036 female cases with ductal and lobular neoplasms breast cancer contributed with 1,164 RNA-seq HTSeq Counts (raw counts) file. 
The filters applied to get those files were: 
1. Data Category: Transcriptome Profiling
2. Data Type: Gene Expression Quantification
3. Experimental Strategy: RNA-Seq
4. Workflow Type: HTSeq - Counts.

#### 3. Downloading the dataset
On March 05, 2020, the files selected were added to the cart on GDC portal. The sample sheet, manifest file, clinical data, and metadata files were downloaded. After that, The GDC Data Transfer Tool [gdc-client](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Getting_Started/ "The GDC Data Transfer Tool") on the same day using the following `Bash` commands.

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

#### 4. Genes in Repair Pathways
A list of all genes invloved in the foloowing 12 pathways were downloaded from Reactome database.
1. BER Participating Molecules [R-HSA-73884]
2. DNA BypassParticipating Molecules [R-HSA-73893]
3. DNA damage Reversal Participating Molecules [R-HSA-73942]
4. DNA double strand break response Participating Molecules [R-HSA-5693606]
5. Fanconi Anemia Participating Molecules [R-HSA-6783310]
6. HDR through HRR alone Participating Molecules [R-HSA-5693567]
7. HDR through HRR and SSA Participating Molecules [R-HSA-5693567]
8. HDR through MMEJ Participating Molecules [R-HSA-5685939]
9. HDR through SSA Participating Molecules [R-HSA-5693567]
10. MMR Participating Molecules [R-HSA-5358508]
11. NER Participating Molecules [R-HSA-5696398]
12. Non homologus end joining Participating Molecules [R-HSA-5693571]

Because some genes are involved in multiple pathways, the Rscript `proteins.R` was used to get the unique gene symbol from the 12 pathwys: 310 genes; they were saved in `repair.genes.csv`

#### 5. Analysis Pipeline
1. Read the raw data of all RNA-seq expression files into a single dataframe >> save it as `exp.rna.rds` object.
2. Collect all protein UniProt IDs of the DNA repair pathways of interest in one array >> save it as `proteins.ids.rds` object.
3. convert the transcript IDs to uniprot IDs and aggregate dublicate by averaging.
4. Select the rows of the proteins of interest only >> save it as `exp.rna.repair.agg.rds` object.
5. convert file IDs to Sample IDs in the working dataframe based on the sample datasheet.
6. subsetting the normal, tumor, and metastatic samples to different dfs based on the sample datasheet.
7. using DESeq2 for DGE analysis.
