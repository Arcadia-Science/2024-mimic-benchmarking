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

results <- read_tsv(args$input_results, show_col_types = FALSE)

if(nrow(results)  > 0){
results <- results %>%
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

host_pdb_lddt <- read_tsv(args$input_host_lddt, show_col_types = FALSE) %>%
  select(protid, pdb_plddt = pdb_confidence) %>%
  distinct()

query_pdb_lddt <- read_tsv(args$input_query_lddt, show_col_types = FALSE) %>%
  select(protid, viro3d_filename, pdb_plddt = pdb_confidence) %>%
  distinct()

query_metadata <- read_tsv(args$input_query_metadata, show_col_types = FALSE) %>%
  clean_names() %>%
  # names might not match, so edit to original file names
  mutate(viro3d_filename = str_remove(string = structure_file, pattern = ".*(?=[CE]F)"),
         viro3d_filename = str_remove(string = viro3d_filename, pattern = "(?<=_relaxed).*")) %>%
  left_join(query_pdb_lddt, by = "viro3d_filename") %>%
  rename_with(.cols = everything(), function(x){paste0("query_", x)}) %>%
  select(-query_viro3d_filename) %>%
  distinct()
  
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
  left_join(query_metadata, by = c("query" = "query_protid")) %>%
  dplyr::relocate(query_virus_name) %>%
  dplyr::relocate(all_of(c("host_gene_names_primary", "host_function_cc",
                           "host_tissue_specificity", "host_subcellular_location_cc")),
                  .after = ttmscore) %>%
  arrange(desc(alntmscore)) %>%
  distinct()


write_tsv(results, args$output)
} else {
file.create(args$output)
}
