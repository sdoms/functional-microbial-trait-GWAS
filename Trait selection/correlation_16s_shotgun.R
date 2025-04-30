setwd("Documents/Research/Experiments/shotgun_hybrids/")
library(tidyverse)
library(PurrfectPlots)
library(ggpubr)

load("results_shotgun_mapping/abund_ps.Rdata")
# abundances are already normalized
name_trans <- readxl::read_excel("name_translation_file.xlsx")
metadata <- readxl::read_excel("metadata_extractions.xlsx")


all_meta <- name_trans %>%
  full_join(metadata, by = (c("subject_name" = "LabID"))) %>%
  mutate(Mouse_name = gsub(x = Mouse_name, pattern = "\\/", replacement = "."))


sample_names <- sample_names(abund_ps)
new_names <- c()
for (name in sample_names) {
  new_name <- all_meta %>%
    filter(DS_name == name) %>%
    pull(Mouse_name)
  new_names <- c(new_names, new_name)
}
sample_names(abund_ps) <- new_names

ps_bac <- subset_taxa(abund_ps, Kingdom == "Bacteria")
ps_bac_F2 <- subset_samples(ps_bac, sample_names(ps_bac) %in% F2)
ps_bac_gen <- tax_glom(ps_bac_F2, taxrank = "Genus")

tops <- names(sort(taxa_sums(ps_bac_gen), TRUE)[1:20])
ps_top_gen <- prune_taxa(tops, ps_bac_gen)


plot_bar(ps_top_gen, x = "Mouse_name", fill = "Genus") +
  geom_bar(stat = "identity", aes(color = Genus)) +
  scale_color_purrs("disney") + scale_fill_purrs("disney") + theme_cat() +
  theme(axis.text.x = element_blank()) + facet_wrap(~Strain, scales = "free_x") + theme(legend.position = "top")

load("../Final_QTL_mapping/Phenotyping_27.02.2020/psdada2_genus.RData")
taxa_names(psdada2.gen.glom.rar)

# problem is that the naming of taxa is not the same => need to harmonize using external database

# Install taxize if not already installed
if (!requireNamespace("taxize", quietly = TRUE)) {
  install.packages("taxize")
}

# Load the library
library(taxize)
# Example genera from RDP training set 16
rdp_genera <- gsub("G_", "", taxa_names(psdada2.gen.glom.rar))
rdp_genera <- rdp_genera[-starts_with("unclassified_", vars = rdp_genera)]
# Example genera from Kraken2
kraken2_genera <- tax_table(ps_bac_gen)[, "Genus"]

# Function to retrieve NCBI Taxonomy IDs
get_tax_ids <- function(genera) {
  sapply(genera, function(genus) {
    tryCatch(
      taxize::get_uid(genus),
      error = function(e) NA # Handle errors gracefully
    )
  })
}

# Get NCBI IDs for each dataset
rdp_ncbi_ids <- get_tax_ids(rdp_genera)
kraken2_ncbi_ids <- get_tax_ids(kraken2_genera)

# Create a data frame for mapping
mappingRDP <- data.frame(
  RDP_Genus = rdp_genera,
  RDP_TaxID = rdp_ncbi_ids,
  stringsAsFactors = FALSE
)

mappingKraken2 <- data.frame(
  Kraken2_Genus = kraken2_genera,
  Kraken2_TaxID = kraken2_ncbi_ids,
  stringsAsFactors = FALSE
)
mapping <- inner_join(mappingRDP, mappingKraken2, by = c("RDP_TaxID" = "Kraken2_TaxID"), na_matches = "never")


psdada2.gen.glom.rar.overlap <- subset_taxa(psdada2.gen.glom.rar, gsub("G_", "", taxa_names(psdada2.gen.glom.rar)) %in% mapping$RDP_Genus)
# Rename RDP genera to match Kraken2 genera using the mapping
# make a column name with Genus species
taxdat <- as.data.frame(tax_table(psdada2.gen.glom.rar.overlap))
taxdat$Kraken2_Genus <- mapping$Genus[match(gsub("G_", "", taxa_names(psdada2.gen.glom.rar.overlap)), mapping$RDP_Genus)]
tax_table(psdada2.gen.glom.rar.overlap) <- tax_table(as.matrix(taxdat))

psdada2_16S_Kraken2_tax <- tax_glom(psdada2.gen.glom.rar.overlap, taxrank = "Kraken2_Genus")
ps_bac_gen_subset <- subset_taxa(ps_bac_gen, Genus %in% mapping$Genus)

# only DNA from 16S
psdada2_16S_Kraken2_tax_DNA <- subset_samples(psdada2_16S_Kraken2_tax, DNA_RNA == "DNA")
otu_16S <- psmelt(psdada2_16S_Kraken2_tax_DNA)
otu_Kraken2 <- psmelt(ps_bac_gen_subset)

otu_16S |> select(Name, Kraken2_Genus, Abundance)
otu_Kraken2 |> select(Sample, Genus, Abundance)
all_otu <- inner_join(otu_16S |> select(Name, Kraken2_Genus, Abundance),
  otu_Kraken2 |> select(Sample, Genus, Abundance),
  by = c("Name" = "Sample", "Kraken2_Genus" = "Genus"), suffix = c("_16S", "_shotgun")
)
cor(all_otu$Abundance_16S, all_otu$Abundance_shotgun, method = "spearman")
cor.test(all_otu$Abundance_16S, all_otu$Abundance_shotgun, method = "spearman")
cor_by_genus <- all_otu |>
  group_by(Kraken2_Genus) |>
  summarise(
    cor = cor(Abundance_16S, Abundance_shotgun, method = "spearman"),
    cor.p = cor.test(Abundance_16S, Abundance_shotgun, method = "spearman")
  )

# core
core_genera <- read.table("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Phenotyping_27.02.2020/tables_core/dna_core_genus.txt")$taxa
core_genera <- gsub("G_", "", core_genera[1:17])
cor_by_genus |> filter(Kraken2_Genus %in% all_of(core_genera))
