aggregate_rows <- function(df.in, agg.var){
  df.in.data <- df.in[-dim(df.in)[2]]
  df.in.data <- apply(df.in.data, 2, as.numeric)
  df.in.agg <- aggregate(df.in.data, list(agg.var), FUN=mean)
  rownames(df.in.agg) <- df.in.agg$Group.1
  df.in.agg <- df.in.agg[-1]   ### This is our final DF
  return(df.in.agg)
}