# TCGA-BRCA Dataset Analysis Code in R

## The Overexpression of DNA Repair Genes in Invasive Ductal and Lobular Breast Cancers: Implications on Precision Medicine

### Abstract

In the era of precision medicine, analyzing the transcriptomic profile is essential to tailor the appropriate therapy for the patients. In this study, we analyzed the differential gene expression between two invasive breast cancer subtypes; infiltrating ductal carcinoma (IDC) and lobular carcinoma (LC) using RNA-Seq data from the TCGA-BRCA project (|log2FC| > 1 and adjusted P-value < 0.05). We found that 3854 are differentially expressed between normal ductal tissues (normal_ductal) and IDC. In addition, upon comparing IDC to LC, 663 genes were also differentially expressed. Since DNA repair players are known to be associated with breast cancer and their expression levels can affect the patients’ response to therapy dramatically, we focused on their analysis in the same samples. We here report that 36 DNA repair players are overexpressed in a significant number of both IDC and LC patient samples. The majority of these players are known to be associated with patients’ response or resistance to therapy. Despite the overexpression of the DNA repair genes in a significant number of patients, we observed a noticeable variation between expression levels in individual patients. Therefore, we believe that the expression of the DNA repair genes should be analyzed carefully in individual patients rather than a cohort. Upon performing a differential expression analysis of miRNAs in the same samples, we could also observe remarkable variations between patients of the same cancer subtype.

In conclusion, we report that a significant number of IDC and LC patients overexpress 36 DNA repair genes, which could affect their response to different treatments and confer resistance to certain drugs. We finally propose that the future of cancer diagnostics and therapy will inevitably depend on high-throughput genomic and transcriptomic data analysis. However, the analysis will have to be performed on individual patients rather than a big set of patients’ samples to ensure that the best treatment is determined and therapy resistance is reduced. 

### Keywords: 

Invasive ductal carcinoma, invasive lobular carcinoma, DNA repair, precision medicine and RNA-Seq.

### Aim

The aim of this project is to analyze the change of gene expression levels in lobular and ductal neoplasms breast cancer subtypes, with focus on personalized approach for all genes involved in repair pathways. 

### Data Retrieval

Data were obtained from The Cancer Genome Atlas Database, [TCGA-BRCA Project](https://portal.gdc.cancer.gov/projects/TCGA-BRCA "TCGA-BRCA Project Page") as described in [DataRetrieval](./DataRetrieval.md) in detail. 

### Data Selection

The Samples' sheet and clinical data sheet were augmented for samples selection. Code is in [scripts/sample.sheet_clinical_prep R file](./scripts/sample.sheet_clinical_prep.R).

### RNA-seq Analysis Pipeline

The following scripts in order:

1. [Main RNA-seq Analysis File](./scripts/main.rna.R)
2. [Comparing IDC and Ductal Normal for RNA-seq Results](./scripts/idc.R)
3. [Comparing IDC to LC for RNA-seq Results](./scripts/idc.lc.R)

### miRNA-seq Analysis Pipeline

The following scripts in order:

1. [Main miRNA-seq Analysis File](./scripts/main.mirna.R)
2. [Comparing IDC and Ductal Normal for miRNA-seq Results](./scripts/mirna-idc.R)
3. [Comparing IDC to LC for miRNA-seq Results](./scripts/mirna.idclc.R)

