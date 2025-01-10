library(tidyverse)
library(janitor)
library(openxlsx)
library(optparse)

option_list <- list(
  make_option(c("--host_organism"), type="character",
              help="Host organism to select structures for."),
  make_option(c("--host_metadata"), type="character",
              help="Path to host organism metadata CSV."),
  make_option(c("--nomburg_metadata"), type="character",
              help="Path to Nomburg et al. supplementary table 1 XLSX."),
  make_option(c("--virushostdb"), type="character",
              help="Path to KEGG virushostdb TSV."),
  make_option(c("--output_tsv"), type="character",
              help="Path to output TSV file."),
  make_option(c("--output_txt"), type="character",
              help="Path to output TXT file.")
)

args <- parse_args(OptionParser(option_list=option_list))

HOST_TAX_ID <- read_csv(args$host_metadata) %>%
  filter(organism == args$host_organism) %>%
  pull(taxon_id)

#virushostdb <- read_tsv("outputs/viral/viral_lineages.tsv") %>%
virushostdb <- read_tsv(args$virushostdb) %>%
  clean_names() %>%
  rename("ncbi_virus_lineage"= x15, "ncbi_virus_tax_id_lineage"=x16) %>%
  filter(host_tax_id == HOST_TAX_ID) %>%
  # change the taxon_id so that it ends with a ";" to prevent partial matching of taxon_ids when grepping
  mutate(ncbi_virus_tax_id_lineage = paste0(ncbi_virus_tax_id_lineage, ";"))

# Read in supplementary table from Nomburg, J., Doherty, E.E., Price, N. et al.
# Birth of protein folds and functions in the virome. Nature 633, 710–717 (2024).
# https://doi.org/10.1038/s41586-024-07809-y
# Table includes names and taxonimic identifiers for viral proteins in the database.
# supp_table1 <- read.xlsx("inputs/41586_2024_7809_MOESM4_ESM.xlsx", sheet = 1)
supp_table1 <- read.xlsx(args$nomburg_metadata, sheet = 1)

# Parse viral protein metadata information from protein names.
supp_table1_parsed <- supp_table1 %>%
  # Fix row that causes parsing error.
  mutate(cluster_member = gsub("putative_membrane__secreted_protein", "putative_membrane_secreted_protein", cluster_member)) %>%
  separate(col = cluster_member, into = c("gene_name", "ncbi_id", "virus", "taxon_id"), sep = "__", remove = FALSE)

# Some nomburg viruses seem to have legacy or wrong taxonomy IDs. Correct these.
# We can correct these no matter what host we're working with
supp_table1_parsed <- supp_table1_parsed %>%
  mutate(taxon_id = ifelse(taxon_id == 11269, "3052505", taxon_id),     # Marburg_marburgvirus
         taxon_id = ifelse(taxon_id == 11623, "3052303", taxon_id),     # Lymphocytic_choriomeningitis_mammarenavirus
         taxon_id = ifelse(taxon_id == 11628, "3052317", taxon_id)) %>% # Machupo_mammarenavirus
  # change the taxon_id so that it ends with a ";" to prevent partial matching of taxon_ids when grepping
  mutate(taxon_id = paste0(";", taxon_id, ";"))

# when host is human, do additional correction that we encoutered through various
# curation efforts.
if(HOST_TAX_ID == 9606){
  # Some viruses that definitely infect humans are in nomburg but not in KEGG.
  # Create a vector containing these viruses so we can include them too.
  other_human_viruses <- c(";11569;", # Thogotovirus_thogotoense
                           ";11274;", # Piry_virus
                           ";11280;", # Vesicular_stomatitis_New_Jersey_virus
                           ";11984;") # Southampton_virus)
  human_infecting_viral_lineages <- c(virushostdb$ncbi_virus_tax_id_lineage, other_human_viruses)
  
  # filter the nomburg viral structures to those that belong to human lineages
  supp_table1_filtered <- supp_table1_parsed %>%
    filter(sapply(taxon_id, function(x) any(grepl(pattern = x,
                                                  x = human_infecting_viral_lineages))))
} else {
  # otherwise, if the host is not human, just work with the annotations in the KEGG table
  # filter the nomburg viral structures to those that belong to host lineages
  supp_table1_filtered <- supp_table1_parsed %>%
    filter(sapply(taxon_id, function(x) any(grepl(pattern = x,
                                                  x = virushostdb$ncbi_virus_tax_id_lineage))))
}

supp_table1_filtered <- supp_table1_filtered %>%
  # remove TaxonID as it is redundant, and remove subcluster_rep bc it isn't
  # clearly described in Nomburg et al.
  select(-taxonID, -subcluster_rep) %>%
  # Rename the cluster columns to something more intuitive, since it won't be
  # clear what these are outside of the context of the supplementary table.
  rename(nomburg_cluster_id = cluster_ID, nomburg_cluster_representative = cluster_rep,
         nomburg_cluster_count = cluster_count, nomburg_protein_name = cluster_member) %>%
  mutate(family = ifelse(is.na(family), "undefined_family", family)) %>%
  # generate the filepaths for the structures in the zip archive
  mutate(structure_filepaths = paste0("viral_structures", "/", family, "/", nomburg_protein_name, ".pdb"))

write_tsv(supp_table1_filtered %>% distinct(), args$output_tsv)
write_tsv(supp_table1_filtered %>% select(structure_filepaths) %>% distinct(), args$output_txt, col_names = FALSE)
