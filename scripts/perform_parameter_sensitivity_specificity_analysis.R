library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("--input_pc_metadata"), type="character",
              help="Path to positive control metadata TSV file."),
  make_option(c("--positive_control"), type="character",
              help="Common name of positive control."),
  make_option(c("--output_full"), type="character",
              help="Path to output TSV file containing all results."),
  make_option(c("--output_best"), type="character",
              help="Path to output TSV file containing the parameters with highest youden index.")
)

args <- parse_args(OptionParser(option_list=option_list))

# functions for sensitivity and specificity analysis ----------------------

# The foldseek results are big. Write a function that reads them in but only
# keep the minimum number of columns so that we save on RAM 
read_foldseek <- function(foldseek_path){
  foldseek <- read_tsv(foldseek_path, show_col_types = FALSE) %>%
    select(query, target, host_gene_names_primary, alnlen, alntmscore, qtmscore,
           ttmscore, qlen, tlen, alnlen, qcov, tcov, evalue)
}

read_gtalign <- function(gtalign_path){
  gtalign <- read_tsv(gtalign_path, show_col_types = FALSE) %>%
    select(query, target, host_gene_names_primary, alnlen, alntmscore, qtmscore,
           ttmscore, qlen, tlen, alnlen, qcov, tcov, q2tmscore, t2tmscore)
}

# We'll use dynamic thresholding to assign positive and negative labels that are
# then used to calculate true positive, false positives, true negatives, and
# false negatives. Since we have positive controls, we use a threshold to bound
# the measurement space within which we assess performance. Note that direction
# indication the direction in which to decide whether something is classified as
# a match or not at a given threshold. For TM-scores, the direction should by >,
# the default. For E-values, the direction should be < (less than). 
label_classification_outcomes_by_threshold <- function(df, threshold, score_column, direction = ">") {
  if(direction == ">"){
    # for TM-score, where bigger values indication better matches
    df <- df %>%
      mutate(predicted_label = ifelse(.data[[score_column]] >= threshold, "Positive", "Negative"))
  } else {
    # For evalue, where smaller values indicate better matches
    df <- df %>%
      mutate(predicted_label = ifelse(.data[[score_column]] <= threshold, "Positive", "Negative"))
  }
  
  
  df %>%
    group_by(comparison) %>%
    summarize(tp = sum(positive_control == "Positive" & predicted_label == "Positive", na.rm = TRUE),
              fp = sum(positive_control == "Negative" & predicted_label == "Positive", na.rm = TRUE),
              fn = sum(positive_control == "Positive" & predicted_label == "Negative", na.rm = TRUE),
              # Foldseek isn't always exhaustive. Set this to the maximum number
              # of comparisons by counting comparisons not returned as true 
              # negatives.
              tn = 201920 - tp - fp - fn) %>%
    mutate(sensitivity  = tp / (tp + fn),
           specificity  = tn / (tn + fp),
           youden_index = sensitivity + specificity - 1,
           threshold    = threshold,
           score_column = score_column)
}

# Create a function to always round down, to the number of digits specified.
# This behaves like the floor() function, but allows the user to specify how
# many decimal places to keep.
floor_dec <- function(x, level=1) {
  round(x - 5*10^(-level-1), level)
}

# Determine performance of different sets of structural comparison parameters
# at different thresholds. The goal is to determine the parameter, measurement,
# and threshold combination that leads to the highest sensitivity and
# specificity (youden index) for each positive control structure (ex. IL10).
run_sensitivity_specificity_analysis <- function(df, score_column) {
  min_positive_control_score <- df %>%
    filter(positive_control == "Positive") %>%
    summarize(min_score = min(.data[[score_column]], na.rm = TRUE)) %>%
    pull(min_score) %>%
    # take the floor of the decimal point, so this value is always below the
    # value of the score in the positive controls
    floor_dec(level = 2)
  
  # Calculate across the range from the minimum positive control score
  all_thresholds <- seq(from = min_positive_control_score, to = 1, by = 0.01)
  
  df_metrics <- map_dfr(
    all_thresholds, 
    ~ label_classification_outcomes_by_threshold(df, .x, score_column)
  ) %>%
    arrange(desc(youden_index))
  
  return(df_metrics)
}

run_sensitivity_specificity_analysis_evalue <- function(df, score_column) {
  max_positive_control_evalue <- df %>%
    filter(positive_control == "Positive") %>%
    summarize(max_evalue = max(.data[[score_column]], na.rm = TRUE)) %>%
    pull(max_evalue)
  
  # Generate a range of thresholds from 0 to the maximum positive control evalue
  # plus a buffer
  all_thresholds <- seq(from = 0, to = max_positive_control_evalue * 10, length.out = 50)
  
  df_metrics <- map_dfr(
    all_thresholds, 
    ~ label_classification_outcomes_by_threshold(df, .x, score_column, direction = "<")
  ) %>%
    arrange(desc(youden_index))
  
  return(df_metrics)
}

# read in and format metadata and results ---------------------------------

# Read in metadata about structures
structure_metadata <- read_tsv(args$input_pc_metadata, show_col_types=FALSE) %>%
  mutate(query = str_remove(string = structure_file, pattern = "\\.pdb"), .after = structure_file) %>%
  select(-structure_file)

# Identify the type of mimic being analyzed
mimic_type <- structure_metadata %>%
  select(target_gene, mimic_type) %>%
  mutate(target_gene = tolower(target_gene)) %>%
  distinct() %>%
  filter(target_gene == args$positive_control) %>%
  pull(mimic_type)

# single mimic analysis ---------------------------------------------------

if(mimic_type == "single_mimic"){
  # Constmimic_type# Construct path strings from positive control input value
  foldseek_glob <- paste0("outputs/human/foldseek/", args$positive_control, "/processed/*tsv")
  gtalign_glob <- paste0("outputs/human/gtalign/", args$positive_control, "/processed/*tsv")
  
  # Read in structure comparison results
  foldseek_results <- Sys.glob(foldseek_glob) %>%
    set_names() %>%
    map_dfr(read_foldseek, .id = "comparison") %>%
    mutate(comparison = gsub(".tsv", "", basename(comparison)))
  
  gtalign_results <- Sys.glob(gtalign_glob) %>%
    set_names() %>%
    map_dfr(read_gtalign, .id = "comparison") %>%
    mutate(comparison = gsub(".tsv", "", basename(comparison)))
  
  # Combine results across method types
  all_results <- bind_rows(foldseek_results, gtalign_results)
  
  # Join results to metadata and label whether a match is a true hit or not.
  all_results <- left_join(all_results, structure_metadata, by = "query") %>%
    mutate(positive_control = ifelse(target_uniprot == target, "Positive", "Negative"))
  
  # run analysis on different measurements  ---------------------------------
  
  # assess alntmscore, qtmscore, and ttmscore for both foldseek and gtalign results
  score_columns <- c("alntmscore", "qtmscore", "ttmscore")
  tmscore_results <- map(score_columns, ~ run_sensitivity_specificity_analysis(all_results, .x))
  names(tmscore_results) <- score_columns
  
  # Filter to foldseek alignment type 2 (the only type that has evalues) and
  # evaluate accuracy at different evalue thresholds.
  alignmenttype2_results <- all_results %>% 
    filter(str_detect(string = comparison, pattern = "alignmenttype2"))
  
  score_columns <- c("evalue")
  evalue_results <- map(score_columns, ~ run_sensitivity_specificity_analysis_evalue(alignmenttype2_results, .x))
  names(evalue_results) <- score_columns
  
  # evaluate q2tmscore and t2tmscore for gtalign results
  gtalign_results <- all_results %>%
    filter(str_detect(string = comparison, pattern = "gtalign"))
  
  score_columns <- c("q2tmscore", "t2tmscore")
  twotmscore_results <- map(score_columns, ~ run_sensitivity_specificity_analysis(gtalign_results, .x))
  names(twotmscore_results) <- score_columns
  
  # pull out results and select the best ones -------------------------------
  
  all_results <- bind_rows(tmscore_results[["alntmscore"]],
                           tmscore_results[["qtmscore"]],
                           tmscore_results[["ttmscore"]],
                           evalue_results[["evalue"]],
                           twotmscore_results[['q2tmscore']],
                           twotmscore_results[['t2tmscore']]) %>%
    arrange(desc(youden_index))
  
  best_result <- all_results %>%
    filter(youden_index == max(youden_index, na.rm = TRUE))
  
  write_tsv(all_results, args$output_full)
  write_tsv(best_result, args$output_best)
} else {
  file.create(args$output_full)
  file.create(args$output_best)
}
