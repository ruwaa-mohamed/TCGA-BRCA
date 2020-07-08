get_df_from_gdc_mirna <- function(mirna.path, file_ext){
  files <- list.files(path=mirna.path, recursive=TRUE, pattern=file_ext)
  
  # Read the first file 
  file <- files[1]
  file.id <- strsplit(file,"/")[[1]][1]
  temp <- read.csv(paste(mirna.path, file, sep="/"), header=TRUE, sep="\t")[, c("miRNA_ID", "read_count")]
  
  # Put it in DF *exp.df*
  exp.df <- temp
  rownames(exp.df) <- exp.df$miRNA_ID
  exp.df <- exp.df[-1]
  colnames(exp.df) <- file.id
  
  # Add all other files iteratively
  for (i in 2:length(files)) {
    file <- files[i]
    file.id <- strsplit(file,"/")[[1]][1]
    temp <- read.csv(paste(mirna.path, file, sep="/"), header=TRUE, sep="\t")[, c("miRNA_ID", "read_count")]
    temp <- temp[-1] # remove the first column.. we already have it in the DF!
    colnames(temp) <- c(file.id)
    exp.df <- cbind(exp.df, temp)
  }
  return(exp.df)
}