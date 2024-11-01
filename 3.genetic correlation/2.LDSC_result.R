library(stringr)
library(dplyr)

extract_genetic_correlation <- function(file_path) {
  lines <- readLines(file_path)
  gc_indexes <- grep("Genetic Correlation", lines)
  gc_results <- list()
  
  if (length(gc_indexes) > 0) {
    for (gc_index in gc_indexes) {
      gc_lines <- lines[gc_index:(gc_index+2)]
      gc_results[[length(gc_results) + 1]] <- gc_lines
    }
  } else {
    gc_results[[1]] <- list(NA, NA, NA)
  }
  
  return(gc_results)
}

parse_data_to_df <- function(data, group) {
  genetic_correlation <- z_score <- p_value <- NA
  if (!is.na(data[[1]])) {
    genetic_correlation <- str_extract(data[1], "Genetic Correlation: [-0-9.]+ \\([-0-9.]+\\)")
    z_score <- as.numeric(str_extract(data[2], "(?<=Z-score: )[-0-9.]+"))
    p_value <- as.numeric(str_extract(data[3], "(?<=P: )[-0-9.e]+"))
  }
  
  df <- data.frame(Group = group, Genetic_Correlation = genetic_correlation,
                   Z_score = z_score, P_value = p_value)
  return(df)
}

file_paths <- list.files(pattern = "*.log", full.names = TRUE)

final_df <- data.frame(Group = character(), Genetic_Correlation = character(),
                       Z_score = numeric(), P_value = numeric())

for (file_path in file_paths) {
  gc_data_list <- extract_genetic_correlation(file_path)
  
  group <- basename(file_path)
  group <- str_match(group, "^(.+)\\.sumstats-(.+)\\.sumstats\\.results\\.log$")[, c(2,3)]
  group <- paste(group, collapse = "-")
  
  for (gc_data in gc_data_list) {
    df <- parse_data_to_df(gc_data, group)
    final_df <- rbind(final_df, df)
  }
}

final_df <- final_df %>%
  group_by(Group) %>%
  arrange(desc(Genetic_Correlation), desc(Z_score), desc(P_value)) %>%
  filter(if (all(is.na(Genetic_Correlation))) TRUE else !is.na(Genetic_Correlation)) %>%
  distinct(Group, .keep_all = TRUE) %>%
  ungroup()

final_df$Genetic_Correlation <- gsub("Genetic Correlation: ", "", final_df$Genetic_Correlation)

fdr <- p.adjust(final_df$P_value, method = "bonferroni")

final_df<- data.frame(final_df, FDR = fdr)

print(final_df)

write.csv(final_df, "LDSC.csv", row.names = FALSE)
