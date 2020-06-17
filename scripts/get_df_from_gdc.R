get_df_from_gdc <- function(data.path, file_ext){
  files <- list.files(path=data.path, recursive=TRUE, pattern=file_ext)
  
  # Read the first file 
  file <- files[1]
  file.id <- strsplit(file,"/")[[1]][1]
  gz.con <- gzfile(file.path(data.path, files[1]))
  temp <- read.table(gz.con, header=FALSE)
  
  # Put it in DF *exp.df*
  exp.df <- temp
  rownames(exp.df) <- exp.df[,1]
  exp.df <- exp.df[-1]
  colnames(exp.df) <- c(file.id)
  
  # Add all other files iteratively
  for (i in 2:length(files)) {
    file <- files[i]
    file.id <- strsplit(file,"/")[[1]][1]
    gz.con <- gzfile(file.path(data.path, files[i]))
    temp <- read.table(gz.con, header=FALSE)
    temp <- temp[-1] # remove the first column.. we already have it in the DF!
    colnames(temp) <- c(file.id)
    exp.df <- cbind(exp.df, temp)
  }
  # rm(list = file.id, gz.con, temp, i)
  # rm(list = ls(file))
  return(exp.df)
}