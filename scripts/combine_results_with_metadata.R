library(tidyverse)
library(janitor)
library(optparse)

option_list <- list(
  make_option(c("--input_results"), type="character",
              help="Path to TSV file from foldseek or gtalign (parsed)."),
  make_option(c("--input_human_metadata"), type="character",
              help="Path to human metadata CSV file."),
  make_option(c("--input_host_lddt"), type="character",
              help="Path to host structure quality measurement TSV file."),
  make_option(c("--input_query_metadata"), type="character",
              help="Path to query metadata TSV file."),
  make_option(c("--input_query_lddt"), type="character",
              help="Path to query structure quality measurement TSV file."),
  make_option(c("--output"), type="character",
              help="Path to output TSV file.")
)

args <- parse_args(OptionParser(option_list=option_list))

results <- read_tsv(args$input_foldseek_results, show_col_types = FALSE) %>%
#results <- read_tsv("~/Downloads/tmp.tsv", show_col_types = FALSE) %>%
  mutate(query = str_remove(string = query, pattern = "\\.pdb"),
         target = str_remove(string = target, pattern = "\\.pdb")) %>%
  # Calculate the alntmscore (alignment TM-score). 
  # Foldseek outputs the incorrect alntmscore in some fraction of results, and
  # gtalign does not output the alntmscore
  mutate(tmraw1 = qtmscore * qlen,
         tmraw2 = ttmscore * tlen,
         tmraw = (tmraw1 + tmraw2) / 2,
         alntmscore = tmraw/alnlen) %>%
  select(-tmraw1, -tmraw2, -tmraw)

#host_pdb_lddt <- read_tsv("inputs/human_proteome_pdb_structure_quality.tsv") %>%
host_pdb_lddt <- read_tsv(args$input_host_lddt, show_col_types = FALSE) %>%
  select(protid, pdb_plddt = pdb_confidence) %>%
  distinct()

#query_pdb_lddt <- read_tsv("benchmarking_data/positive_controls/positive_controls_plddt.tsv") %>%
query_pdb_lddt <- read_tsv(args$input_query_lddt, show_col_types = FALSE) %>%
  select(protid, pdb_plddt = pdb_confidence) %>%
  distinct()

#query_metadata <- read_tsv("inputs/viral/viral_structure_metadata.tsv") %>%
query_metadata <- read_tsv(args$input_query_metadata, show_col_types = FALSE) %>%
  clean_names() %>%
  left_join(query_pdb_lddt, by = c("nomburg_protein_name" = "protid")) %>%
  rename_with(.cols = everything(), function(x){paste0("query_", x)}) %>%
  select(-query_structure_filepaths) %>%
  distinct()
  
#host_metadata <- read_csv("inputs/human_metadata_combined.csv.gz") %>%
host_metadata <- read_csv(args$input_human_metadata, show_col_types = FALSE) %>% 
  rename(protid = uniprot) %>%
  left_join(host_pdb_lddt, by = c("protid")) %>%
  rename(host_protid = protid, host_pdb_plddt = pdb_plddt) %>%
  rename_with(.cols = starts_with("uniprot"), function(x){gsub("uniprot_", "host_", x)}) %>%
  group_by(host_protid) %>%
  slice_head(n = 1) %>%
  distinct()

results <- results %>%
  left_join(host_metadata, by = c("target" = "host_protid")) %>%
  left_join(query_metadata, by = c("query" = "query_nomburg_protein_name")) %>%
  dplyr::relocate(query_species) %>%
  dplyr::relocate(all_of(c("host_gene_names_primary", "host_function_cc",
                           "host_tissue_specificity", "host_subcellular_location_cc")),
                  .after = ttmscore) %>%
  arrange(desc(alntmscore)) %>%
  distinct()


write_tsv(foldseek_results, args$output)
