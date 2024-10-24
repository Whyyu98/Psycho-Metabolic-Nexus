setwd("E:/GXN/COLOC")

library(data.table)
library(dplyr)
library(coloc)
library(ieugwasr)

snp_data <- fread("SNP.txt")

a <- fread("COLOC/susie/EUR2m.bim")
a <- a[, c("V2", "V5", "V6")]
colnames(a) <- c("SNP", "A1_ref", "A2_ref")

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
      result_list[[i]] <- data.frame(CHR = chr_val, BP = NA, SNP = NA, coloc.post.prob = NA, group = paste(chr_val, start_val, sep = ":"))
      result_res_list[[i]] <- data.frame(SNP = NA, SNP.PP.H4 = NA)
      next  
    }
    
    data <- merge(data3, data4, by = "SNP")
    data <- merge(a, data, by = "SNP")
    data <- data[!duplicated(data$SNP), ]
    data <- data %>% filter(!(VAR.x >= 1e-06 & VAR.x < 1e-05) & !(VAR.y >= 1e-06 & VAR.y < 1e-05))
    for (i in 1:nrow(data)) {
      if ((data$A1.x[i] == data$A2_ref[i] && data$A2.x[i] == data$A1_ref[i])) {
        data$BETA.x[i] <- -data$BETA.x[i]
        temp <- data$A1.x[i]
        data$A1.x[i] <- data$A2.x[i]
        data$A2.x[i] <- temp
      } else if (!(data$A1.x[i] == data$A1_ref[i] && data$A2.x[i] == data$A2_ref[i])) {
        cat("Warning: SNP", data$SNP[i], "in A1.x and A2.x does not match reference alleles\n")
      }
    }
    
    for (i in 1:nrow(data)) {
      if ((data$A1.y[i] == data$A2_ref[i] && data$A2.y[i] == data$A1_ref[i])) {
        data$BETA.y[i] <- -data$BETA.y[i]
        temp <- data$A1.y[i]
        data$A1.y[i] <- data$A2.y[i]
        data$A2.y[i] <- temp
      } else if (!(data$A1.y[i] == data$A1_ref[i] && data$A2.y[i] == data$A2_ref[i])) {
        cat("Warning: SNP", data$SNP[i], "in A1.y and A2.y does not match reference alleles\n")
      }
    }
    
    data3 <- data[, c("BETA.x", "VAR.x", "SNP")]
    data4 <- data[, c("BETA.y", "VAR.y", "SNP")]
    
    eqtlld_result <- ld_matrix(
      data$SNP,
      plink_bin = "./plink.exe",
      keep-allele-order,
      bfile = "./COLOC/susie/EUR2m",
      with_alleles = FALSE
    )
    
    snp_order_ld <- rownames(eqtlld_result)
    
    data3 <- data3[match(snp_order_ld, data3$SNP), ]
    data4 <- data4[match(snp_order_ld, data4$SNP), ]
    
    colnames(data3) <- c("beta", "varbeta", "snp")
    colnames(data4) <- c("beta", "varbeta", "snp")
    data3 <- as.list(data3)
    data4 <- as.list(data4)
    data3$type <- "cc"
    data4$type <- "cc"
    data3$LD <- eqtlld_result
    data4$LD <- eqtlld_result
    
    D3 <- runsusie(data3)
    D4 <- runsusie(data4)
    
    res <- coloc.susie(dataset1 = D3, dataset2 = D4)
    
    res_summary <- res$summary
    res_summary$group <- paste(chr_val, start_val, sep = ":")
    result_list[[i]] <- res_summary
  }
  
  result_df <- bind_rows(result_list)
  
  print(result_df)
  
  write.csv(result_df, paste0("susie_/", group_val, "_summary.csv"), row.names = FALSE)
}
