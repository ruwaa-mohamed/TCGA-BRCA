# TCGA-BRCA

The aim of this project is to analyze the change of expression levels of proteins that are involved in DNA repair.

### Methodology:
#### Downloading the dataset
The dataset **TCGA-BRCA** was downloaded from TCGA database as manifest files on March 05, 2020.
After that, GDC-data transfer tool was used to download the data using the following `Bash` commands:

```bash
cwd=$(pwd)
cd ~/Downloads/gdc-client_v1.5.0_Ubuntu_x64
./gdc-client download -h

./gdc-client download -m $cwd/raw_data/rna-seq_gdc_manifest.2020-03-05.txt -d $cwd/raw_data/RNA-seq/

./gdc-client download -m $cwd/raw_data/mirna-seq_gdc_manifest.2020-03-05.txt -d $cwd/raw_data/miRNA-seq/

./gdc-client download -m $cwd/raw_data/meth_gdc_manifest.2020-03-05.txt -d $cwd/raw_data/Methylation-Array/
```
#### Analysis Pipeline
1. Read the raw data of all RNA-seq expression files into a single dataframe >> save it as `exp.rna.rds` object.
2. Collect all protein UniProt IDs of the DNA repair pathways of interest in one array >> save it as `proteins.ids.rds` object.
3. convert the transcript IDs to uniprot IDs and aggregate dublicate by averaging.
4. Select the rows of the proteins of interest only >> save it as `exp.rna.repair.agg.rds` object.
5. convert file IDs to Sample IDs in the working dataframe based on the sample datasheet.
6. subsetting the normal, tumor, and metastatic samples to different dfs based on the sample datasheet.
7. using DESeq2 for DGE analysis.
