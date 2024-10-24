library(data.table)
library(dplyr)
library(coloc)

snp_data <- fread("SNP.txt")

unique_groups <- unique(snp_data$Group)

for (group_val in unique_groups) {
  parts <- unlist(strsplit(group_val, "-"))
  group1 <- parts[1]
  group2 <- parts[2]
  
  data1 <- fread(paste0("E:/GXN/datasets/PD/", group1, ".tsv.gz"))
  data2 <- fread(paste0("E:/GXN/datasets/data/", group2, ".tsv.gz"))
  data1$VAR <- data1$SE^2
  data2$VAR <- data2$SE^2
  
  group_snp_data <- snp_data %>% filter(Group == group_val)
  
  result_list <- list()
  result_res_list <- list()
  
  for (i in 1:nrow(group_snp_data)) {
    chr_val <- group_snp_data$CHR[i]
    start_val <- group_snp_data$BP[i]
    
    data3 <- data1 %>% filter(CHR == chr_val, BP >= (start_val - 100000), BP <= (start_val + 100000))
    data4 <- data2 %>% filter(CHR == chr_val, BP >= (start_val - 100000), BP <= (start_val + 100000))
    
    if (nrow(data3) == 0 || nrow(data4) == 0) {
      cat(paste("No overlapping SNPs for CHR", chr_val, "Group", group_val, "\n"))
      result_list[[i]] <- data.frame(CHR = chr_val, BP = NA, SNP = NA, coloc.post.prob = NA)
      result_res_list[[i]] <- data.frame(SNP = NA, SNP.PP.H4 = NA)
      next  
    }
    
    data <- merge(data3, data4, by = "SNP")
    data <- data[!duplicated(data$SNP), ]
    data <- data %>% mutate(BETA.y = ifelse(A1.x == A1.y, BETA.y, -BETA.y))
    
    data3 <- data[, c("BETA.x", "VAR.x", "SNP")]
    data4 <- data[, c("BETA.y", "VAR.y", "SNP")]
    colnames(data3) <- c("beta", "varbeta", "snp")
    colnames(data4) <- c("beta", "varbeta", "snp")
    data3 <- as.list(data3)
    data4 <- as.list(data4)
    data3$type <- "cc"
    data4$type <- "cc"
    
    res <- coloc.abf(data3, data4, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    
    result_list[[i]] <- res$summary
    
    max_h4_rows <- res$results %>% filter(SNP.PP.H4 == max(SNP.PP.H4, na.rm = TRUE))
    max_h4_row <- max_h4_rows[sample(nrow(max_h4_rows), 1), ]
    result_res_list[[i]] <- max_h4_row
  }
  
  result_df <- bind_rows(result_list)
  result_res_df <- bind_rows(result_res_list)
  
  print(result_df)
  print(result_res_df)
  
  write.csv(result_df, paste0("result/", group_val, "_summary.csv"), row.names = FALSE)
  write.csv(result_res_df, paste0("SNP/", group_val, "_max_SNP.PP.H4.csv"), row.names = FALSE)
}
