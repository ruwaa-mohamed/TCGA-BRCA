plotting_gene <- function(gene){
  boxplot(t(assays(dds.run.plot[gene])[["counts"]])~dds.run.plot$Sample.Type, 
          range=0, las=1, log='y', boxwex=.4,
          at=c(0.0, 0.5), cex.lab = 1.5, cex.axis = 1.5,
          xlab=NULL, ylab=NULL,
          col=c("darksalmon", "darkred"))
  stripchart(t(assays(dds.run.plot[gene])[["counts"]])~dds.run.plot$Sample.Type, 
             vertical=TRUE, method='jitter', add=TRUE, pch=23, col="black", cex=1, at=c(0, 0.5)) 
}