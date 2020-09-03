plotting_gene <- function(dds.run, gene, l='y', atr="counts"){
  boxplot(t(assays(dds.run[gene])[[atr]])~dds.run$Sample.Type, 
          range=0, las=1, log=l, boxwex=.4,
          at=c(0.0, 0.5), #cex.lab = 1.5, cex.axis = 1.5,
          xlab=NULL, ylab=NULL,
          main=gene,
          col=c("darksalmon", "darkred"))
  stripchart(t(assays(dds.run[gene])[[atr]])~dds.run$Sample.Type, 
             vertical=TRUE, method='jitter', add=TRUE, 
             pch=20, col="black", cex=0.5, at=c(0, 0.5)) 
}