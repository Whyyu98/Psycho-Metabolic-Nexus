library(data.table)
library(TwoSampleMR)
library(plinkbinr)
library(ieugwasr)

pd_files <- list.files("E:/GXN/datasets/data", pattern = "*.tsv.gz", full.names = TRUE)

output_dir <- "MeTs-exp"
dir.create(output_dir, showWarnings = FALSE)

for (pd_file in pd_files) {
  a <- fread(pd_file)
  a$P <- as.numeric(a$P)
  
  pd_file_name <- tools::file_path_sans_ext(basename(pd_file))
  
  p_thresh <- 5e-08

  b <- subset(a, P < p_thresh)
  
  if (nrow(b) < 3) {
    p_thresh <- 5e-06
    b <- subset(a, P < p_thresh)
  }
  
  exposure_file <- paste0(output_dir, "/exposure_", pd_file_name, ".csv")
  write.csv(b, file = exposure_file, row.names = FALSE)
  
  exp_dat <- read_exposure_data(filename = exposure_file, sep = ",", 
                                snp_col = "SNP", beta_col = "BETA", se_col = "SE", 
                                effect_allele_col = "A1", other_allele_col = "A2", 
                                eaf_col = "FRQ", pval_col = "P")
  
  exp_dat_clumped <- ld_clump(clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99, 
                              pop = "EUR", 
                              dplyr::tibble(rsid = exp_dat$SNP, pval = exp_dat$pval.exposure, id = exp_dat$id.exposure), 
                              plink_bin = "E:/SJ/plink_win64_20230116/plink.exe", 
                              bfile = "E:/GXN/SMR/g1000_eur/g1000_eur")
  
  if (nrow(exp_dat_clumped) < 3) {
    p_thresh <- 5e-06
    b <- subset(a, P < p_thresh)
    write.csv(b, file = exposure_file, row.names = FALSE)
    exp_dat <- read_exposure_data(filename = exposure_file, sep = ",", 
                                  snp_col = "SNP", beta_col = "BETA", se_col = "SE", 
                                  effect_allele_col = "A1", other_allele_col = "A2", 
                                  eaf_col = "FRQ", pval_col = "P")
    exp_dat_clumped <- ld_clump(clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99, 
                                pop = "EUR", 
                                dplyr::tibble(rsid = exp_dat$SNP, pval = exp_dat$pval.exposure, id = exp_dat$id.exposure), 
                                plink_bin = "E:/SJ/plink_win64_20230116/plink.exe", 
                                bfile = "E:/GXN/SMR/g1000_eur/g1000_eur")
  }
  
  exp_dat_clumped <- subset(exp_dat, SNP %in% exp_dat_clumped$rsid)
  
  file_names <- list.files("E:/GXN/datasets/PD", pattern = "*.tsv.gz", full.names = TRUE)
  
  for (file in file_names) {
    c <- fread(file)
    
    file_name <- tools::file_path_sans_ext(basename(file))
    
    d <- merge(exp_dat_clumped, c, by.x = "SNP", by.y = "SNP")
    
    if (nrow(d) < 3) {
      p_thresh <- 5e-06
      b <- subset(a, P < p_thresh)
      write.csv(b, file = exposure_file, row.names = FALSE)
      exp_dat <- read_exposure_data(filename = exposure_file, sep = ",", 
                                    snp_col = "SNP", beta_col = "BETA", se_col = "SE", 
                                    effect_allele_col = "A1", other_allele_col = "A2", 
                                    eaf_col = "FRQ", pval_col = "P")
      exp_dat_clumped <- ld_clump(clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99, 
                                  pop = "EUR", 
                                  dplyr::tibble(rsid = exp_dat$SNP, pval = exp_dat$pval.exposure, id = exp_dat$id.exposure), 
                                  plink_bin = "E:/SJ/plink_win64_20230116/plink.exe", 
                                  bfile = "E:/GXN/SMR/g1000_eur/g1000_eur")
      exp_dat_clumped <- subset(exp_dat, SNP %in% exp_dat_clumped$rsid)
      d <- merge(exp_dat_clumped, c, by.x = "SNP", by.y = "SNP")
    }
    
    d <- d[is.finite(d$BETA) & is.finite(d$SE), ]
    
    write.csv(d, file = paste0(output_dir, "/outcome_", pd_file_name, "_", file_name, ".csv"), row.names = FALSE)
    
    outcome_file <- paste0(output_dir, "/outcome_", pd_file_name, "_", file_name, ".csv")
    outcome_dat <- read_outcome_data(snps = exp_dat_clumped$SNP, filename = outcome_file,
                                     sep = ",", snp_col = "SNP", beta_col = "BETA",
                                     se_col = "SE", effect_allele_col = "A1",
                                     other_allele_col = "A2", pval_col = "P")
    
    dat <- harmonise_data(exposure_dat = exp_dat_clumped, outcome_dat = outcome_dat)
    
    mr_res <- mr(dat)
    e <- generate_odds_ratios(mr_res = mr_res)
    fdr <- p.adjust(e$pval, method = "BH")
    e <- data.frame(e, FDR = fdr)
    
    output_file <- paste0(output_dir, "/", pd_file_name, "_", file_name, ".csv")
    write.csv(e, file = output_file, row.names = FALSE)
  }
}
