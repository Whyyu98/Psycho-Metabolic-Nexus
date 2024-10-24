library(MRlap)
library(data.table)

groups <- fread("Group.txt", header = FALSE)$V1

samples <- fread("sample.txt", header = TRUE, col.names = c("NAMES", "N"))

pd_datasets <- c("MDD", "BD", "SCZ", "ADHD", "ANX", "ASD", "OCD", "PTSD")

for (group in groups) {
  parts <- unlist(strsplit(group, "-"))
  data1_name <- parts[1]
  data2_name <- parts[2]
  
  if (data1_name %in% pd_datasets) {
    data1_path <- paste0("E:/GXN/datasets/PD/", data1_name, ".tsv.gz")
  } else {
    data1_path <- paste0("E:/GXN/datasets/data/", data1_name, ".tsv.gz")
  }
  
  data1 <- fread(data1_path)
  data1 <- data1[, c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P")]
  colnames(data1) <- c("snp", "chr", "pos", "a1", "a2", "beta", "se", "p")
  data1$N <- samples[samples$NAMES == data1_name, N]
  
  if (data2_name %in% pd_datasets) {
    data2_path <- paste0("E:/GXN/datasets/PD/", data2_name, ".tsv.gz")
  } else {
    data2_path <- paste0("E:/GXN/datasets/data/", data2_name, ".tsv.gz")
  }
  
  data2 <- fread(data2_path)
  data2 <- data2[, c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P")]
  data2$P <- as.numeric(data2$P)
  colnames(data2) <- c("snp", "chr", "pos", "a1", "a2", "beta", "se", "p")
  data2$N <- samples[samples$NAMES == data2_name, N]
  
  A <- MRlap(exposure = data1,
             exposure_name = data1_name,
             outcome = data2,
             outcome_name = data2_name,
             MR_pruning_dist = 1000,
             MR_pruning_LD = 0.001,
             ld = "E:/GXN/LDSC/eur_w_ld_chr",
             hm3 = "E:/GXN/LDSC/w_hm3.snplist")
  
  corrected_effect <- A$MRcorrection$corrected_effect
  corrected_effect_se <- A$MRcorrection$corrected_effect_se
  corrected_effect_p <- A$MRcorrection$corrected_effect_p
  df <- data.frame(
    corrected_effect = corrected_effect,
    corrected_effect_se = corrected_effect_se,
    corrected_effect_p = corrected_effect_p
  )
  
  output_file <- paste0("result/", group, ".csv")
  write.csv(df, file = output_file, row.names = FALSE)
}
