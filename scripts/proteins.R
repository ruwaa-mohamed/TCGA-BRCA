library(stringr)
prot.dir.path <- "proteinsinthepathway/"
pathways.files <- list.files(path = prot.dir.path, pattern = '.tsv')
 
proteins.temp <- read.table(file = paste(prot.dir.path, pathways.files[1], sep = ""), sep = "\t", header = TRUE)

proteins.df <- proteins.temp
for (i in 2:length(pathways.files)) {
  proteins.temp <- read.table(file = paste(prot.dir.path, pathways.files[i], sep = ""), sep = "\t", header = TRUE)
  proteins.df <- rbind(proteins.df, proteins.temp)
}
proteins.df.uniq <- unique(proteins.df)

proteins.df.MoleculeName <- str_split_fixed(proteins.df.uniq$MoleculeName, " ", 2)
proteins.df.MoleculeName <- as.data.frame(proteins.df.MoleculeName)
colnames(proteins.df.MoleculeName) <- c("uniprot", "symbol")
repair.genes <- unique(proteins.df.MoleculeName["Symbol"])
write.csv(repair.genes, "saved_objects/repair.genes.csv")

rm(list=c("prot.dir.path", "pathways.files", "proteins.temp", "proteins.df", "proteins.df.uniq", "proteins.df.MoleculeName"))
