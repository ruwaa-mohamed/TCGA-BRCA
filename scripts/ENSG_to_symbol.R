library(org.Hs.eg.db)
ENSG_to_symbol <- function(df.in){
  df.ensemble.id <- sapply(rownames(df.in), function(x) strsplit(as.character(x),"\\.")[[1]][1])
  df.in <- cbind(df.ensemble.id, df.in)
  mapper <- mapIds(org.Hs.eg.db, keys=df.ensemble.id, keytype="ENSEMBL", column="SYMBOL", multiVals="first")
  mapper.df <- as.data.frame(mapper)
  mapper.df <- cbind(rownames(mapper.df), mapper.df)
  names(mapper.df) <- c("df.ensemble.id", "symbol")
  df.out <- merge(df.in, mapper.df ,by="df.ensemble.id", all.x=TRUE)
  df.out <- df.out[-1]
  df.out <- df.out[ ! is.na(df.out$symbol),]
  # df.out <- aggregate_rows(df.out, agg.var=df.out$symbol)
  return(df.out)
}