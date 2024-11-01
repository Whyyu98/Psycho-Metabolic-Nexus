library(data.table)
library(dplyr)

file_list <- list.files("E:/GXN/MTAG/MTAG/", pattern = "*.txt", full.names = TRUE)

for (file_path in file_list) {

  file_name <- basename(file_path)
  file_name_no_ext <- tools::file_path_sans_ext(file_name)
  
  data1 <- fread(file_path)
  data1 <- data1[, c("SNP", "mtag_pval")]
  colnames(data1) <- c("SNP", "P.MTAG")
  
  data2_file_path <- paste0("E:/GXN/CPASSOC/results/", file_name_no_ext, "-SHet.txt")
  data2 <- fread(data2_file_path)
  
  data <- merge(data1, data2, by = "SNP")
  
  parts <- unlist(strsplit(file_name_no_ext, "-"))
  group1 <- parts[1]
  group2 <- parts[2]
  
  PD_file_path <- paste0("E:/GXN/datasets/PD/", group1, ".tsv.gz")
  PD <- fread(PD_file_path)
  PD <- PD[, c("SNP", "P")]
  colnames(PD) <- c("SNP", "P.PD")
  
  data <- merge(data, PD, by = "SNP")
  
  MeTs_file_path <- paste0("E:/GXN/datasets/data/", group2, ".tsv.gz")
  MeTs <- fread(MeTs_file_path)
  MeTs <- MeTs[, c("SNP", "P")]
  colnames(MeTs) <- c("SNP", "P.MeTs")
  
  data <- merge(data, MeTs, by = "SNP")
  
  f <- data %>%
    filter(P.MTAG < 5e-8 ) %>%
    filter(P.SHet < 5e-8 ) %>%
    filter(P.MTAG < P.PD) %>%
    filter(P.MTAG < P.MeTs)%>%
    filter(P.PD < 1e-4) %>%
    filter(P.MeTs < 1e-4)
  
  output_file <- paste0("MTAG_CPASSOC/", file_name_no_ext, ".txt")
  
  fwrite(f, output_file, sep = "\t")
}
