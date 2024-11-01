library(data.table)
library(dplyr)

group_data <- fread("Group.txt")

pd_folder <- "E:/GXN/datasets/PD"
data_folder <- "E:/GXN/datasets/data"

result_data <- data.frame()

for (group in unique(group_data$Group)) {
  group_parts <- strsplit(group, "-")[[1]]
  group_pd <- group_parts[1]
  group_mets <- group_parts[2]
  
  pd_file <- file.path(pd_folder, paste0(group_pd, ".tsv.gz"))
  pd_data <- fread(pd_file)
  
  mets_file <- file.path(data_folder, paste0(group_mets, ".tsv.gz"))
  mets_data <- fread(mets_file)
  
  snps <- group_data %>% filter(Group == group) %>% pull(SNP)
  
  pd_subset <- pd_data %>%
    filter(SNP %in% snps) %>%
    select(SNP, BETA, P,SE) %>%
    rename(BETA.PD = BETA, P.PD = P,SE.PD=SE)
  
  mets_subset <- mets_data %>%
    filter(SNP %in% snps) %>%
    select(SNP, BETA, P,SE) %>%
    rename(BETA.MeTs = BETA, P.MeTs = P,SE.MeTs=SE)
  
  combined_data <- pd_subset %>%
    inner_join(mets_subset, by = "SNP")
  
  combined_data$VAR.PD=combined_data$SE.PD^2
  combined_data$VAR.MeTs=combined_data$SE.MeTs^2
  combined_data$Group <- group
  
  result_data <- rbind(result_data, combined_data)
}
print(result_data)

fwrite(result_data, "extracted_SNPs.csv", row.names = FALSE)
