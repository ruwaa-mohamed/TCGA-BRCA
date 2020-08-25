################################################################################
## Merging DEGs and DEMs by dplyr
library(dplyr)

hsa_MTI <- read.xlsx("raw_data/miRNA-DBs/hsa_MTI.xlsx")
substring(hsa_MTI$miRNA, 7) = 'r'

mir.gene.pairs <- unique(hsa_MTI[, c("miRNA", "Target.Gene")])
mir.gene.pairs.rep <- mir.gene.pairs[mir.gene.pairs$Target.Gene %in% repair.genes$symbol,]
################################################################################
## IDC DEGs
idc.degs.rep <- as.data.frame(res.idc.degs.rep)
idc.degs.rep <- cbind(Target.Gene=rownames(idc.degs.rep), idc.degs.rep)
idc.degs.rep <- mutate(idc.degs.rep, degs = log2FoldChange > 0)
idc.degs.rep$degs <- as.factor(idc.degs.rep$degs)
levels(idc.degs.rep$degs) <- c("up")
table(idc.degs.rep$degs)
idc.degs.rep <- idc.degs.rep[,c("Target.Gene", "degs")]
################################################################################
## IDC DEMs
idc.dems <- as.data.frame(mirna.res.idc.dems) 
idc.dems <- cbind(miRNA=rownames(idc.dems), idc.dems)
idc.dems <- mutate(idc.dems, dems = log2FoldChange > 0)
idc.dems$dems <- as.factor(idc.dems$dems)
levels(idc.dems$dems) <- c("down", "up")
table(idc.dems$dems)
idc.dems <- idc.dems[,c("miRNA", "dems")]
################################################################################
## Linking IDC DEGs and DEMs
idc.degs.dems <- mir.gene.pairs.rep
idc.degs.dems <- full_join(idc.degs.dems, idc.dems, by="miRNA")
idc.degs.dems <- full_join(idc.degs.dems, idc.degs.rep, by="Target.Gene")
idc.degs.dems <- idc.degs.dems[complete.cases(idc.degs.dems),]
write.table(idc.degs.dems, "saved_objects/idc.degs.dems.csv", sep=",", quote=FALSE, row.names=FALSE)
################################################################################
## There's no IDCLC DEGs Repair!