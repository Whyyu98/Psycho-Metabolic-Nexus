library(data.table)
library(dplyr)

sample_size_data <- fread("sample.txt", header = FALSE, col.names = c("Group", "SampleSize"))

file_list <- list.files("DATA/", pattern = "*.txt", full.names = TRUE)

for (file_path in file_list) {
  file_name <- basename(file_path)
  file_name_no_ext <- tools::file_path_sans_ext(file_name)

  parts <- unlist(strsplit(file_name_no_ext, "-"))
  group1 <- parts[1]
  group2 <- parts[2]
  
  sample_size1 <- sample_size_data[sample_size_data$Group == group1,]$SampleSize
  sample_size2 <- sample_size_data[sample_size_data$Group == group2,]$SampleSize
  
  if (length(sample_size1) > 0 & length(sample_size2) > 0) {
    SampleSize <- c(sample_size1, sample_size2)
    
    X <- read.table(file_path, header = TRUE, sep = "\t", check.names = FALSE)
    X <- X %>% distinct(SNP, .keep_all = TRUE)  
    rownames(X) <- X[,1]  
    X <- X[,-1]  
    
    CorrMatrix <- cor(X)
    
    source("CPASSOC.R")
    x <- SHet(X = X, SampleSize = SampleSize, CorrMatrix = CorrMatrix, correct = 1, isAllpossible = TRUE)
    SHet <- as.matrix(x)
    
    p.shet <- pchisq(SHet, df = 1, ncp = 0, lower.tail = FALSE)
    
    result <- data.frame(SNP = rownames(p.shet), P.SHet = p.shet)
    colnames(result) <- c("SNP", "P.SHet")
    
    output_file <- paste0("results/", file_name_no_ext, "-SHet.txt")

    write.table(result, file = output_file, row.names = FALSE, quote = FALSE, sep = "\t")
  } else {
    cat("done:", group1, group2, "\n")
  }
}
