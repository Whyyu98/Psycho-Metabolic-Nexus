library(dplyr)
library(data.table)

tissue_groups <- list(
  "Adrenal gland" = c("fetal_adrenal_gland"),
  "Testis" = c("testis", "fetal_testes"),
  "Brain" = c("brain_hippocampus", "brain", "cerebellar", "fetal_brain"),
  "Nervous" = c("nervous", "spinal_cord", "fetal_spinal_cord"),
  "Heart" = c("heart", "fetal_heart"),
  "Stem cell" = c("es_cell", "ips_cell", "blastula"),
  "Lung" = c("lung", "fetal_lung", "embryonic_lung"),
  "Stomach" = c("stomach", "fetal_stomach"),
  "Kidney" = c("kidney", "fetal_kidney"),
  "Urothelium" = c("urothelium"),
  "Epithelium" = c("epithelium"),
  "Eye" = c("eye"),
  "Intestine" = c("colon", "fetal_intestine,_large", "fetal_intestine,_small"),
  "Pancreas/spleen" = c("pancreas", "pancreatic_duct", "fetal_spleen"),
  "Skin" = c("skin", "foreskin", "fetal_skin"),
  "Prostate" = c("prostate"),
  "Vascular endothelial" = c("blood_vessel"),
  "Thymus" = c("fetal_thymus"),
  "Liver" = c("liver"),
  "Connective" = c("connective", "fibroblast"),
  "Muscle" = c("muscle", "fetal_muscle", "fetal_muscle,_lower_limb", "fetal_muscle,_trunk", "fetal_muscle,_upper_trunk"),
  "Breast" = c("breast"),
  "Cervix" = c("cervix"),
  "Embryo/uterus" = c("uterus", "myometrium", "fetal_membrane", "fetal_placenta"),
  "Blood/immune" = c("blood", "multi-tissue"),
  "Others" = c("gingival")
)

output_dir <- "garfield-data/output"
subdirs <- list.dirs(output_dir, full.names = TRUE, recursive = FALSE)

for (subdir in subdirs) {
  subdir_name <- basename(subdir)
  
  file_path <- file.path(subdir, paste0("garfield.test.", subdir_name, ".out"))
  
  data <- fread(file_path)
  
  data$Group <- sapply(data$Tissue, function(tissue) {
    group <- NA
    for (key in names(tissue_groups)) {
      if (tissue %in% tissue_groups[[key]]) {
        group <- key
        break
      }
    }
    return(group)
  })
  
  data <- data %>%
    group_by(Group) %>%
    mutate(FDR = p.adjust(Pvalue, method = "bonferroni")) %>%
    ungroup()
  
  filtered_data <- data %>% filter(FDR < 0.05)
  
  significant_proportion <- data %>%
    group_by(Group) %>%
    summarise(
      total_annotations = n(),
      significant_annotations = sum(FDR < 0.05),
      proportion_significant = significant_annotations / total_annotations
    )
  
  significant_proportion <- na.omit(significant_proportion)
  
  # 保存结果到文�?
  fwrite(filtered_data, file.path("out", paste0(subdir_name, ".csv")), row.names = FALSE)
  fwrite(significant_proportion, file.path("out_", paste0(subdir_name, ".csv")), row.names = FALSE)
}
