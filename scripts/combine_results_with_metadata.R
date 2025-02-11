library(tidyverse)
library(janitor)
library(optparse)

option_list <- list(
  make_option(c("--input_results"), type="character",
              help="Path to TSV file from foldseek or gtalign (parsed)."),
  make_option(c("--target"), type="character",
              help="Database used for search. Either human or viral."),
  make_option(c("--input_human_metadata"), type="character",
              help="Path to human metadata CSV file."),
  make_option(c("--input_host_lddt"), type="character",
              help="Path to host structure quality measurement TSV file."),
  make_option(c("--input_viral_metadata"), type="character",
              help="Path to viral metadata TSV file."),
  make_option(c("--output"), type="character",
              help="Path to output TSV file.")
)

args <- parse_args(OptionParser(option_list=option_list))

results <- read_tsv(args$input_results, show_col_types = FALSE)

if(nrow(results)  > 0){
  results <- results %>%
    mutate(query = str_remove(string = query, pattern = "\\.pdb"),
           query = str_remove(string = query, pattern = "\\-F1-model_v4"),
           query = str_remove(string = query, pattern = "^AF-"),
           target = str_remove(string = target, pattern = "\\.pdb")) %>%
    # Calculate the alntmscore (alignment TM-score). 
    # Foldseek outputs the incorrect alntmscore in some fraction of results, and
    # gtalign does not output the alntmscore
    mutate(tmraw1 = qtmscore * qlen,
           tmraw2 = ttmscore * tlen,
           tmraw = (tmraw1 + tmraw2) / 2,
           alntmscore = tmraw/alnlen) %>%
    select(-tmraw1, -tmraw2, -tmraw)
  
  host_pdb_lddt <- read_tsv(args$input_host_lddt, show_col_types = FALSE) %>%
    select(protid, pdb_plddt = pdb_confidence) %>%
    distinct()

  viral_metadata <- read_tsv(args$input_viral_metadata, show_col_types = FALSE) %>%
    clean_names() %>%
    # names might not match, so edit to original file names
    mutate(protid = str_remove(string = structure_file, pattern = "\\.pdb$")) %>%
    rename_with(.cols = everything(), function(x){paste0("viral_", x)}) %>%
    distinct() 
 
  host_metadata <- read_csv(args$input_human_metadata, show_col_types = FALSE) %>% 
    rename(protid = uniprot) %>%
    left_join(host_pdb_lddt, by = c("protid")) %>%
    rename(host_protid = protid, host_pdb_plddt = pdb_plddt) %>%
    rename_with(.cols = starts_with("uniprot"), function(x){gsub("uniprot_", "host_", x)}) %>%
    group_by(host_protid) %>%
    slice_head(n = 1) %>%
    distinct()
  
  if(args$target == "human"){
    results <- results %>%
      left_join(host_metadata, by = c("target" = "host_protid")) %>%
      left_join(viral_metadata, by = c("query" = "viral_protid"))
  } else { 
    results <- results %>%
      left_join(host_metadata, by = c("query" = "host_protid")) %>%
      left_join(viral_metadata, by = c("target" = "viral_protid"))
  }
  results <- results %>%
      dplyr::relocate(viral_virus_name) %>%
      dplyr::relocate(all_of(c("host_gene_names_primary", "host_function_cc",
                               "host_tissue_specificity", "host_subcellular_location_cc")),
                      .after = ttmscore) %>%
      arrange(desc(alntmscore)) %>%
      distinct()

  write_tsv(results, args$output)
} else {
  file.create(args$output)
}
