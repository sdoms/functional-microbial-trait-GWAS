# Set working directory to results folder
setwd("~/Documents/Research/Experiments/shotgun_hybrids/mapping_results/Results_genes/")

# Load necessary library
library(tidyverse)

# Load pre-saved GWAS results
load(file = "../all_results_GWAS_functional_traits_SW.RData")

# Create output directory for extracted gene results
output_dir <- "KEGG_parsing_genes_interval_per_trait/genes"
dir.create(output_dir, showWarnings = F)

# Get list of unique traits from the results
traits <- unique(EC_BRITE_GBM_GMM_PCoA$trait)

# Loop through all traits
for (tt in traits) {
  # Check if KEGG parsing results already exist for this trait
  if (file.exists(paste0("KEGG_parsing_genes_interval_per_trait", "/", tt, "_KEGG_results_genes.csv"))) {
    # Load trait-specific KEGG results
    final_results <- read_csv(paste0("KEGG_parsing_genes_interval_per_trait", "/", tt, "_KEGG_results_genes.csv"))

    # ----------------- BRITE and EC parsing -----------------

    # Split BRITE column into separate rows per gene
    brite_split <- final_results %>%
      filter(!is.na(BRITE), BRITE != "NA") %>%
      dplyr::select(external_gene_name, BRITE) |>
      separate_rows(BRITE, sep = ";") %>%
      mutate(BRITE = trimws(BRITE)) |> # remove leading/trailing spaces
      distinct() # remove duplicates

    # Split EC column into separate rows per gene
    ec_split <- final_results %>%
      filter(!is.na(EC)) %>%
      dplyr::select(external_gene_name, EC) |>
      separate_rows(EC, sep = ";") %>%
      mutate(EC = trimws(EC)) |> # remove leading/trailing spaces
      distinct()

    # ----------------- Gene matching logic -----------------

    # Check if trait name corresponds to an EC or BRITE term
    if (startsWith(tt, "EC_level")) {
      # Remove prefix to get the actual EC term
      tt_name <- gsub("EC_level_", "", tt)

      # Try matching the trait to the cleaned EC categories
      if (tt_name %in% gsub(pattern = "[^[:alnum:] ]", "_", ec_split$EC)) {
        ec_names <- gsub(pattern = "[^[:alnum:] ]", "_", ec_split$EC)
        gene <- ec_split[which(ec_names == tt_name), 1]
      } else if (tt_name %in% gsub(pattern = " ", "", ec_split$EC)) {
        ec_names <- gsub(pattern = " ", "", ec_split$EC)
        gene <- ec_split[which(ec_names == tt_name), 1]
      } else {
        gene <- NULL
      }
    } else if ((startsWith(tt, "BRITE"))) {
      # Remove prefix to get the actual BRITE term
      tt_name <- gsub("BRITE_", "", tt)

      # Try matching the trait to the cleaned BRITE categories
      if (tt_name %in% gsub(pattern = "[^[:alnum:] ]", "_", brite_split$BRITE)) {
        brite_names <- gsub(pattern = "[^[:alnum:] ]", "_", brite_split$BRITE)
        gene <- brite_split[which(brite_names == tt_name, 1)]
      } else if (tt_name %in% gsub(pattern = " ", "", brite_split$BRITE)) {
        brite_names <- gsub(pattern = " ", "", brite_split$BRITE)
        gene <- brite_split[which(brite_names == tt_name), 1]
      } else {
        gene <- NULL
      }
    } else {
      gene <- NULL
    }

    # Save matching genes if found
    if (!is.null(gene)) write_csv(as.data.frame(gene), paste0(output_dir, "/genes_with_same_KO_", tt, ".csv"))
  }
}

# ----------------- Plotting section -----------------

# Load full KEGG results table for all genes
kegg_results <- read_csv("kegg_cache.csv")

# Split EC column by gene
ec_split <- kegg_results %>%
  filter(!is.na(EC)) %>%
  select(EntrezID, EC) |>
  separate_rows(EC, sep = ";") %>%
  mutate(EC = trimws(EC)) %>%
  distinct()

# Split BRITE column by gene
brite_split <- kegg_results %>%
  filter(!is.na(BRITE), BRITE != "NA") %>%
  select(EntrezID, BRITE) |>
  separate_rows(BRITE, sep = ";") %>%
  mutate(BRITE = trimws(BRITE)) %>%
  distinct()

# Split REACTION column by gene
reaction_split <- kegg_results %>%
  filter(!is.na(REACTION), REACTION != "NA") %>%
  select(EntrezID, REACTION) |>
  separate_rows(REACTION, sep = ";") %>%
  mutate(REACTION = trimws(REACTION)) %>%
  distinct()

# Split MODULE column by gene
module_split <- kegg_results %>%
  filter(!is.na(MODULE), MODULE != "NA") %>%
  select(EntrezID, MODULE) |>
  separate_rows(MODULE, sep = ";") %>%
  mutate(MODULE = trimws(MODULE)) %>%
  distinct()

# ----------------- EC frequency plot -----------------
if (nrow(ec_split) > 0) {
  ec_count <- ec_split %>%
    group_by(EC) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  ec_count <- ec_count[-which(ec_count$EC == "NA"), ]
  ec_plot <- ec_count |>
    filter(count > 5) |> # only show EC terms with more than 5 genes
    ggplot(aes(x = reorder(EC, count), y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(x = "EC Categories", y = "Frequency") +
    theme_minimal()
  print(ec_plot)
}
ggsave("EC_count_more5.pdf", width = 12, height = 6, plot = ec_plot)

# ----------------- BRITE frequency plot -----------------
if (nrow(brite_split) > 0) {
  brite_count <- brite_split %>%
    group_by(BRITE) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  brite_plot <- brite_count |>
    mutate(BRITE = reorder(BRITE, count)) |>
    slice(-1) |> # remove top BRITE (likely uninformative)
    slice_head(n = 50) |> # top 50 BRITE categories
    ggplot(aes(x = reorder(BRITE, count), y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(x = "BRITE Categories", y = "Frequency") +
    theme_minimal()
  print(brite_plot)
}
ggsave("BRITE_count_top50.pdf", width = 12, height = 6, plot = brite_plot)

# Alternative BRITE plot excluding categories starting with "09" (e.g., diseases)
brite_plot2 <- brite_count |>
  filter(!startsWith(BRITE, "09")) |>
  slice_head(n = 50) |>
  ggplot(aes(x = reorder(BRITE, count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "BRITE Categories", y = "Frequency") +
  theme_minimal()
brite_plot2
ggsave("BRITE_count_top50_without09.pdf", width = 12, height = 6, plot = brite_plot2)

# ----------------- GPCR genes (BRITE category 04030) -----------------
# Load gene annotations from Ensembl
library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Filter Entrez IDs of GPCRs
gpcrs <- brite_split |>
  filter(BRITE == "04030 G protein-coupled receptors [BR:mmu04030]") |>
  pull(EntrezID)

# Retrieve gene names for GPCR Entrez IDs
gpcr_names <- getBM(
  attributes = c("external_gene_name", "entrezgene_id"), filters = "entrezgene_id",
  values = gpcrs, mart = mart
)

# ----------------- REACTION frequency plot -----------------
if (nrow(reaction_split) > 0) {
  reaction_count <- reaction_split %>%
    group_by(REACTION) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  reaction_plot <- ggplot(reaction_count, aes(x = reorder(REACTION, count), y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(x = "Reaction Categories", y = "Frequency", title = tt) +
    theme_minimal()
  print(reaction_plot)
}

# ----------------- MODULE frequency plot -----------------
if (nrow(module_split) > 0) {
  module_count <- module_split %>%
    group_by(MODULE) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  module_plot <- ggplot(module_count, aes(x = reorder(MODULE, count), y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(x = "Module Categories", y = "Frequency") +
    theme_minimal()
  print(module_plot)
}
ggsave("module_count_.pdf", width = 12, height = 10, plot = module_plot)
