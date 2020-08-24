to_tiff <- function(f, n){
  figures_dir = "final_figures"
  tiff(paste(figures_dir, n, sep="/"),
       width = 18, height = 21, units = 'cm', res = 300)
  f
  dev.off()
}
  