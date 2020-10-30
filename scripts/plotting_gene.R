plotting_gene <- function(dds.run, gene, colorslist = c("darksalmon", "darkred"),l='y', atr="counts"){
  boxplot(t(assays(dds.run[gene])[[atr]])~dds.run$group, 
          range=0, las=1, log=l, boxwex=.4,
          at=c(0.0, 0.5, 1.0), #cex.lab = 1.25, cex.axis = 1.25,
          xlab=NULL, ylab=NULL, #ylim = c(2.7, 18.3),
          main=gene,
          col=colorslist)
  stripchart(t(assays(dds.run[gene])[[atr]])~dds.run$group, 
             vertical=TRUE, method='jitter', add=TRUE, 
             pch=20, col="black", cex=0.5, at=c(0.0, 0.5, 1.0)) 
}