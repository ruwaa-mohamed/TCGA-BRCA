## For Sample sheets and Clinical data
################################################################################
## RNA-seq Sample Sheet
sample.sheet <- read.csv("raw_data/females_only.ductal-and-lobular-only.gdc_sample_sheet.2020-06-17.tsv", sep="\t", header=TRUE)
table(sample.sheet$Sample.Type)

## Drop Metastatic from sample.sheet
sample.sheet <- sample.sheet[! sample.sheet$Sample.Type == "Metastatic",]
table(sample.sheet$Sample.Type)

## As factor and Releveling
sample.sheet$Sample.Type <- as.factor(sample.sheet$Sample.Type)
sample.sheet$Sample.Type <- relevel(sample.sheet$Sample.Type, "Solid Tissue Normal")
levels(sample.sheet$Sample.Type) <- c("Normal", "Tumor")
table(sample.sheet$Sample.Type)
################################################################################
## miRNAs-seq Sample Sheet
mirna.sample.sheet <- read.csv("raw_data/gdc_sample_sheet.2020-07-10.tsv", sep="\t", header=TRUE)
table(mirna.sample.sheet$Sample.Type)

## Drop Metastatic from sample.sheet
mirna.sample.sheet <- mirna.sample.sheet[! mirna.sample.sheet$Sample.Type == "Metastatic",]

## As factor and Releveling
mirna.sample.sheet$Sample.Type <- as.factor(mirna.sample.sheet$Sample.Type)
mirna.sample.sheet$Sample.Type <- relevel(mirna.sample.sheet$Sample.Type, "Solid Tissue Normal")
levels(mirna.sample.sheet$Sample.Type) <- c("Normal", "Tumor")
table(mirna.sample.sheet$Sample.Type)
################################################################################
## Clinical Data
clinical <- read.csv("raw_data/clinical.project-TCGA-BRCA.2020-03-05/clinical.tsv", header=TRUE, sep="\t", na="--")
clinical <- clinical[, colSums(is.na(clinical)) != nrow(clinical)]
length(unique(clinical$submitter_id))

## Subsetting by cases
clinical <- clinical[clinical$submitter_id %in% c(sample.sheet$Case.ID, mirna.sample.sheet$Case.ID) ,]
length(unique(clinical$submitter_id))

## Saving for later
write.csv(clinical, 'saved_objects/clinical.all.csv', row.names=FALSE)

## Subsetting by Primary_diagnosis
clinical <- clinical[clinical$primary_diagnosis %in% c("Infiltrating duct carcinoma, NOS", "Lobular carcinoma, NOS"),]
clinical <- unique(clinical[,c("submitter_id", "primary_diagnosis")])

## As factor and Releveling
clinical$primary_diagnosis <- as.factor(clinical$primary_diagnosis)
levels(clinical$primary_diagnosis) <- c("IDC", "LC")
table(clinical$primary_diagnosis)

## Saving for later
write.csv(clinical, 'saved_objects/clinical.csv', row.names=FALSE)
clinical <- read.csv('saved_objects/clinical.csv', header=TRUE)
################################################################################
## Sync Samples Sheets with the new clinical data
## with clinical
sample.sheet <- sample.sheet[sample.sheet$Case.ID %in% clinical$submitter_id,]
mirna.sample.sheet <- mirna.sample.sheet[mirna.sample.sheet$Case.ID %in% clinical$submitter_id,]

## with each other
# sample.sheet <- sample.sheet[sample.sheet$Case.ID %in% mirna.sample.sheet$Case.ID,]
mirna.sample.sheet <- mirna.sample.sheet[mirna.sample.sheet$Case.ID %in% sample.sheet$Case.ID,]

## check
all(sample.sheet$Case.ID %in% mirna.sample.sheet$Case.ID)
all(mirna.sample.sheet$Case.ID %in% sample.sheet$Case.ID)
################################################################################
## adding the primary diagnosis column
mapperIDs <- match(sample.sheet$Case.ID, clinical$submitter_id)
sample.sheet <- cbind(sample.sheet, primary_diagnosis=clinical$primary_diagnosis[mapperIDs])
rm(mapperIDs)

mapperIDs <- match(mirna.sample.sheet$Case.ID, clinical$submitter_id)
mirna.sample.sheet <- cbind(mirna.sample.sheet, primary_diagnosis=clinical$primary_diagnosis[mapperIDs])
rm(mapperIDs)
################################################################################
## save for importing in other scripts
write.csv(sample.sheet, 'saved_objects/sample.sheet.csv', row.names=FALSE)
sample.sheet <- read.csv('saved_objects/sample.sheet.csv', header=TRUE)

write.csv(mirna.sample.sheet, 'saved_objects/mirna.sample.sheet.csv', row.names=FALSE)
mirna.sample.sheet <- read.csv('saved_objects/mirna.sample.sheet.csv', header=TRUE)
################################################################################