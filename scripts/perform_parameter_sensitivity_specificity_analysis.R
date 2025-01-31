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
              help="Path to output TSV file containing the parameters with highest youden index."),
  make_option(c("--output_full_viral"), type="character",
              help="Path to output TSV file containing results for viruses only."),
  make_option(c("--output_best_viral"), type="character",
              help="Path to output TSV file containing the parameters with highest youden index for viruses only.")
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
# df <- all_results
# tp_metadata <- single_mimic_tp_count_all
# threshold <- 0.7
# score_column <- "alntmscore"

label_classification_outcomes_by_threshold <- function(df, tp_metadata, threshold, score_column, direction = ">") {
  
  df <- df %>%
    # join to metadata that reports the number of true positives we actually have
    left_join(tp_metadata, by = "target_gene")
    
  if(direction == ">"){
    # for *TM-score, where bigger values indication better matches
    df <- df %>%
      mutate(predicted_label = ifelse(.data[[score_column]] >= threshold, "Positive", "Negative")) 
  } else {
    # For evalue, where smaller values indicate better matches
    df <- df %>%
      mutate(predicted_label = ifelse(.data[[score_column]] <= threshold, "Positive", "Negative"))
  }
  
  df_summary <- df %>%
    ungroup() %>%
    group_by(comparison) %>%
    summarize(tp = sum(positive_control == "Positive" & predicted_label == "Positive", na.rm = TRUE),
              fp = sum(positive_control == "Negative" & predicted_label == "Positive", na.rm = TRUE),
              # Hack to keep the number of positive controls column. 
              # In this case, should always return only one number.
              num_positive_controls = unique(num_positive_controls)) %>%
    mutate(
      # This was dropping some comparisons if they weren't returned by foldseek.
      # Therefore, we report the number of positive controls we started with
      # in a metadata table and use this number, minus the number of positive
      # controls we detect, to calculate the number of false negatives
      fn = num_positive_controls - tp,
      # Foldseek isn't always exhaustive. Set this to the maximum number
      # of comparisons by counting comparisons not returned as true 
      # negatives.
      tn = 201920 - tp - fp - fn,
      sensitivity  = tp / (tp + fn),
      specificity  = tn / (tn + fp),
      youden_index = sensitivity + specificity - 1,
      threshold    = threshold,
      score_column = score_column) %>%
    select(-num_positive_controls)
  return(df_summary)
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
run_sensitivity_specificity_analysis <- function(df, tp_metadata, score_column) {
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
    ~ label_classification_outcomes_by_threshold(df = df,
                                                 tp_metadata = tp_metadata,
                                                 threshold = .x,
                                                 score_column = score_column)
  ) %>%
    arrange(desc(youden_index))
  
  return(df_metrics)
}

run_sensitivity_specificity_analysis_evalue <- function(df, tp_metadata, score_column) {
  max_positive_control_evalue <- df %>%
    filter(positive_control == "Positive") %>%
    summarize(max_evalue = max(.data[[score_column]], na.rm = TRUE)) %>%
    pull(max_evalue)
  
  # Generate a range of thresholds from 0 to the maximum positive control evalue
  # plus a buffer
  all_thresholds <- seq(from = 0, to = max_positive_control_evalue * 10, length.out = 50)
  
  df_metrics <- map_dfr(
    all_thresholds, 
    ~ label_classification_outcomes_by_threshold(df = df, 
                                                 tp_metadata = tp_metadata,
                                                 threshold = .x,
                                                 score_column = score_column,
                                                 direction = "<")
  ) %>%
    arrange(desc(youden_index))
  
  return(df_metrics)
}

run_full_analysis <- function(df, tp_metadata) {
  # Assess alntmscore, qtmscore, and ttmscore for both foldseek and gtalign results
  score_columns <- c("alntmscore", "qtmscore", "ttmscore")
  tmscore_results <- map(score_columns, ~run_sensitivity_specificity_analysis(df = df, 
                                                                              tp_metadata = tp_metadata,
                                                                              score_column = .x))
  names(tmscore_results) <- score_columns
  
  # Filter to foldseek alignment type 2 and evaluate accuracy at different evalue thresholds
  alignmenttype2_results <- df %>%
    filter(str_detect(string = comparison, pattern = "alignmenttype2"))
  
  score_columns <- c("evalue")
  evalue_results <- map(score_columns, ~run_sensitivity_specificity_analysis_evalue(df = alignmenttype2_results, 
                                                                                    tp_metadata = tp_metadata,
                                                                                    score_column = .x))
  names(evalue_results) <- score_columns
  
  # Evaluate q2tmscore and t2tmscore for gtalign results
  gtalign_results <- df %>%
    filter(str_detect(string = comparison, pattern = "gtalign"))
  
  score_columns <- c("q2tmscore", "t2tmscore")
  twotmscore_results <- map(score_columns, ~run_sensitivity_specificity_analysis(df = gtalign_results,
                                                                                 tp_metadata = tp_metadata,
                                                                                 score_column = .x))
  names(twotmscore_results) <- score_columns
  
  # Pull out results and select the best ones
  all_results <- bind_rows(tmscore_results[["alntmscore"]],
                           tmscore_results[["qtmscore"]],
                           tmscore_results[["ttmscore"]],
                           evalue_results[["evalue"]],
                           twotmscore_results[['q2tmscore']],
                           twotmscore_results[['t2tmscore']]) %>%
    arrange(desc(youden_index))
  
  best_result <- all_results %>%
    filter(youden_index == max(youden_index, na.rm = TRUE))
  
  return(list(all_results = all_results, best_result = best_result))
}

# read in and format metadata and results ---------------------------------

# Read in metadata about structures
structure_metadata <- read_tsv(args$input_pc_metadata, show_col_types=FALSE) %>%
  mutate(query = str_remove(string = structure_file, pattern = "\\.pdb"), .after = structure_file) %>%
  select(-structure_file)

# Identify the type of mimic being analyzed
mimic_type <- structure_metadata %>%
  select(target_gene, mimic_type) %>%
  distinct() %>%
  mutate(target_gene = tolower(target_gene)) %>%
  filter(target_gene == args$positive_control) %>%
  pull(mimic_type)

single_mimic_tp_count_all <- structure_metadata %>%
  filter(mimic_type == "single_mimic") %>%
  group_by(target_gene) %>%
  tally() %>%
  select(target_gene, num_positive_controls = n)

single_mimic_tp_count_viral <- structure_metadata %>%
  filter(!structure_species %in% c("chimp", "human", "mouse", "macaque")) %>%
  filter(mimic_type == "single_mimic") %>%
  group_by(target_gene) %>%
  tally() %>%
  select(target_gene, num_positive_controls = n)

# single mimic analysis ---------------------------------------------------

if(mimic_type == "single_mimic"){
  # Construct path strings from positive control input value
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
  
  # Run analysis on viral and euk results
  viral_and_euk_results <- run_full_analysis(df = all_results,
                                             tp_metadata = single_mimic_tp_count_all)
  write_tsv(viral_and_euk_results$all_results, args$output_full)
  write_tsv(viral_and_euk_results$best_result, args$output_best)
  
  # Run the analysis on only the viral results
  viral_results <- run_full_analysis(df = all_results %>% filter(!structure_species %in% c("chimp", "human", "mouse", "macaque")),
                                     tp_metadata = single_mimic_tp_count_viral)
  write_tsv(viral_results$all_results, args$output_full_viral)
  write_tsv(viral_results$best_result, args$output_best_viral)

} else {
  file.create(args$output_full)
  file.create(args$output_best)
  
  file.create(args$output_full_viral)
  file.create(args$output_best_viral)
}
