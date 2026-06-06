# initialise
library(tidyverse)
setwd("~/GitHub/MBDP_Metagenomics_2026/06_ANVIO/")

# read sample metadata
metadata <- read_tsv("metadata.tsv", show_col_types = FALSE)

# read anvio summary
summary <- tibble()
for (file in list.files(pattern = "*_bins_summary.txt")) {
  cols <- c("MAG", "total_length", "num_contigs", "N50", "GC_content", 
            "percent_completion", "percent_redundancy", "anvio_taxonomy")
  
  df <- read_tsv(file, show_col_types = FALSE)
  df <- mutate(df, bins = str_c(str_extract(file, "ERR500034[2,3]"), bins, sep = "_"))
  df <- mutate(df, anvio_taxonomy = str_c(t_domain, t_phylum, t_class, t_order, t_family, t_genus, t_species, sep = ";"))
  df <- rename(df, MAG = bins)
  df <- select(df, all_of(cols))
  
  summary <- bind_rows(summary, df)
}
summary

# read gtdb taxonomy
gtdb <- tibble()
for (file in list.files(pattern = "gtdbtk.*.summary.tsv")) {
  cols <- c("MAG", "gtdb_taxonomy")
  
  df <- read_tsv(file, show_col_types = FALSE)
  df <- rename(df, MAG = user_genome)
  df <- rename(df, gtdb_taxonomy = classification)
  df <- select(df, all_of(cols))
  
  gtdb <- bind_rows(gtdb, df)
}
gtdb

# read MAG coverage
coverage <- tibble()
for (file in list.files(pattern = "*_mean_coverage.txt")) {
  df <- read_tsv(file, show_col_types = FALSE)
  df <- mutate(df, bins = str_c(str_extract(file, "ERR500034[2,3]"), bins, sep = "_"))
  df <- rename(df, MAG = bins)
  
  coverage <- bind_rows(coverage, df)
}
coverage

# transform to relative abundance
relabund <- mutate(coverage, across(where(is.numeric), function(x) x/sum(x)))

# read anvio metabolism
metabolism <- tibble()
for (file in list.files(pattern = "*_modules.txt")) {
  cols <- c("module", "MAG", "module_name", "module_class", "module_category", "module_subcategory", 
            "module_completeness", "module_is_complete", "enzyme_hits_in_module")
            
  df <- read_tsv(file, show_col_types = FALSE)
  df <- mutate(df, bin_name = str_c(str_extract(file, "ERR500034[2,3]"), bin_name, sep = "_"))
  df <- rename(df, MAG = bin_name)
  df <- rename(df, module_completeness = pathwise_module_completeness)
  df <- rename(df, module_is_complete = pathwise_module_is_complete)
  df <- select(df, all_of(cols))
  
  metabolism <- bind_rows(metabolism, df)
}
metabolism

# transform to presence/absence
metabolism_pa <- select(metabolism, MAG, module_name, module_is_complete)
metabolism_pa <- spread(metabolism_pa, module_name, module_is_complete, fill = FALSE)
metabolism_pa <- mutate(metabolism_pa, across(where(is.logical), as.numeric))
metabolism_pa

# plot MAG taxonomy
full_join(gtdb, relabund, by = "MAG") %>% 
  select(MAG, gtdb_taxonomy, starts_with("ERR")) %>% 
  gather(accession, relabund, -MAG, -gtdb_taxonomy) %>% 
  group_by(gtdb_taxonomy, accession) %>% 
  summarise(relabund = sum(relabund)) %>% 
  ungroup() %>% 
  left_join(metadata, by = "accession") %>% 
  mutate(relabund = na_if(relabund, 0)) %>% 
  filter(str_detect(gtdb_taxonomy, "Unclassified", negate = TRUE)) %>% 
  ggplot(aes(x = accession, 
             y = fct_rev(gtdb_taxonomy))) +
  geom_point(aes(size = relabund,
                 fill = relabund),
             shape = 21) +
  ggforce::facet_row(vars(vegetation), 
                     scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5),
        axis.text.y = element_text(hjust = 0),
        axis.title = element_blank()) +
  scale_fill_viridis_c(labels = scales::percent) +
  scale_size_continuous(guide = "none")

# see what Methanoflorens is doing
full_join(gtdb, metabolism_pa, by = "MAG") %>% 
  filter(str_detect(gtdb_taxonomy, "Methanoflorens")) %>% 
  select(-gtdb_taxonomy) %>% 
  gather(module_name, module_is_complete, -MAG) %>% 
  filter(module_is_complete == 1) %>% 
  print(n = Inf)

# check methane genes
filter(metabolism, MAG == "ERR5000343_Bin_31") %>% 
  filter(str_detect(module_name, "Methanogen")) %>% 
  filter(module_is_complete == TRUE) %>% 
  select(module_name, enzyme_hits_in_module) %>% 
  pull(enzyme_hits_in_module)
