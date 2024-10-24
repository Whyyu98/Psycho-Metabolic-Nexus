library(garfield)
library(data.table)

files <- list.files(path = "PD-AD", pattern = "*.txt", full.names = TRUE)

output_dir <- "garfield-data/pval"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

for (file_path in files) {
  data <- fread(file_path)
  data_list <- split(data, data$CHR_PD)
  
  file_prefix <- tools::file_path_sans_ext(basename(file_path))
  
  subfolder <- file.path(output_dir, file_prefix)
  if (!dir.exists(subfolder)) {
    dir.create(subfolder)
  }
  
  for (chr in names(data_list)) {
    chr_data <- data_list[[chr]]
    
    chr_data <- chr_data[, .(BP_PD, P_PD_AD_shared)]
    chr_data <- chr_data[order(BP_PD)]
    
    file_name <- file.path(subfolder, paste0("chr", chr))
    fwrite(chr_data, file_name, col.names = FALSE, sep = "\t")
  }
}
