## Reading the list of proteins of interest
folder.path <- "C:/Users/Salma Bargal/Desktop/Zewail City/Bioinformatics Projects 2020/TCGA-BRCA/Menna/Poteins in pathways/"
pathways.files <- list.files(path = folder.path, pattern = '.tsv')
 
pathways.files.temp <- read.table(file = paste(folder.path, pathways.files[1], sep = ""), sep = "\t", header = TRUE)

pathways.files.df <- pathways.files.temp
for (i in 2:length(pathways.files)) {
  pathways.files.temp <- read.table(file = paste(folder.path, pathways.files[i], sep = ""), sep = "\t", header = TRUE)
  pathways.files.df <- rbind(pathways.files.df, pathways.files.temp)
}
pathways.files.df.uniq <- unique(pathways.files.df)

library(stringr)

prot.df.new <- str_split_fixed(pathways.files.df.uniq$MoleculeName, " ", 2)
prot.df.new <- as.data.frame(prot.df.new)
colnames(prot.df.new) <- c("Uniprot", "Symbol")
proteins.ids <- unique(prot.df.new["Symbol"])

