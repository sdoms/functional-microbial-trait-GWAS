setwd("~/Documents/Research/Experiments/shotgun_hybrids/mapping_results/Results_genes/")
library(KEGGREST)
library(biomaRt)
library(tidyverse)
library(future.apply)
# library(purrr)
library(furrr)

# Set up parallel processing
plan(multisession, workers = 1) # Adjust the number of workers as needed
# cannot use multisession because the Kegg api is rate-limiting

# Load data
load(file = "../all_results_GWAS_functional_traits_SW.RData")

# Set up biomaRt connection for mouse genes
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

output_dir <- "KEGG_parsing_genes_interval_per_trait"
dir.create(output_dir, showWarnings = FALSE)

# Caching function for KEGG queries
cache_file <- paste0(output_dir, "/kegg_cache.csv")


# Define a function that fetches KEGG data with rate-limiting
get_kegg_entry <- function(entrez_id) {
  # Add a delay to maintain 3 requests per second rate
  Sys.sleep(0.35) # Sleep for 350 ms

  kegg_entry <- tryCatch(keggGet(paste0("mmu:", entrez_id)), error = function(e) NULL)

  return(kegg_entry)
}

get_kegg_info <- function(entrez_id) {
  # Check if data is cached
  if (entrez_id %in% kegg_cache$EntrezID) {
    return(kegg_cache[kegg_cache$EntrezID == entrez_id, ])
  }

  # Query KEGG and check response
  # Check if KEGG returns anything
  kegg_entry <- get_kegg_entry(entrez_id)

  # Print the result

  if (is.null(kegg_entry) || length(kegg_entry) == 0) {
    cat("No data found for Entrez ID:", entrez_id, "\n")
    return(data.frame(EntrezID = entrez_id, all_info = NA, BRITE = NA, EC = NA, MODULE = NA, REACTION = NA))
  }

  # Process KEGG data with extra checks
  all_info <- if (!is.null(kegg_entry[[1]]$BRITE)) paste(unique(kegg_entry[[1]]$BRITE), collapse = "; ") else NA
  brite_info <- if (!is.na(all_info)) {
    matches <- unique(str_extract_all(all_info, "\\b\\d{5} [^;]*")[[1]])
    paste(matches, collapse = "; ")
  } else {
    NA
  }

  # EC info retrieval with enhanced handling
  ec_info <- if (!is.null(kegg_entry[[1]]$ENZYME)) {
    paste(kegg_entry[[1]]$ENZYME, collapse = "; ")
  } else {
    ec_section <- if (!is.na(all_info)) str_extract(all_info, "(?<=Enzymes \\[BR:mmu01000\\];).*") else NULL
    ec_info <- if (!is.null(ec_section) && nzchar(ec_section)) {
      paste(trimws(str_split(ec_section, ";", simplify = TRUE)[1:4]), collapse = "; ")
    } else {
      NA
    }
  }

  module_info <- if (!is.null(kegg_entry[[1]]$MODULE)) paste(kegg_entry[[1]]$MODULE, collapse = "; ") else NA
  reaction_info <- if (!is.null(kegg_entry[[1]]$REACTION)) paste(kegg_entry[[1]]$REACTION, collapse = "; ") else NA

  # Save results to cache if data is not empty

  result <- data.frame(EntrezID = entrez_id, all_info = all_info, BRITE = brite_info, EC = ec_info, MODULE = module_info, REACTION = reaction_info)

  return(result)
}
# Required libraries
library(KEGGREST)
library(dplyr)
library(purrr)
library(biomaRt)
library(tidyr)
library(ggplot2)
library(readr)

# Define cache file and output directory
cache_file <- "kegg_cache.csv"

dir.create(output_dir, showWarnings = FALSE)

# Define traits
traits <- unique(EC_BRITE_GBM_GMM_PCoA$trait)[1:10]

# Loop over each trait
for (tt in traits) {
  # Load cache if exists, otherwise create a new data frame
  kegg_cache <- if (file.exists(cache_file)) read_csv(cache_file) else data.frame(EntrezID = NA, all_info = NA, BRITE = NA, EC = NA, MODULE = NA, REACTION = NA)

  # Prepare data for the current trait
  trait_info <- EC_BRITE_GBM_GMM_PCoA |> filter(trait == tt)
  all_genes_trait <- na.omit(unique(trimws(c(str_split(trait_info$list_total_genes, ",", simplify = TRUE)))))

  # Retrieve KEGG gene IDs for the gene symbols
  gene_mapping <- tryCatch(
    {
      getBM(
        attributes = c("external_gene_name", "entrezgene_id"), filters = "external_gene_name",
        values = all_genes_trait, mart = mart
      )
    },
    error = function(e) {
      message("No genes found for trait: ", tt)
      return(NULL)
    }
  )

  # If gene_mapping is NULL or has zero rows, skip to the next trait
  if (is.null(gene_mapping) || nrow(gene_mapping) == 0) {
    next
  }

  # Process KEGG info in parallel with rate limit
  kegg_results <- future_map_dfr(gene_mapping$entrezgene_id, get_kegg_info)
  kegg_cache <- rbind(kegg_cache, kegg_results) |> distinct()

  # Save the cache for future sessions
  write_csv(kegg_cache, cache_file)

  # Merge gene symbols and KEGG results
  final_results <- merge(gene_mapping, kegg_results |> distinct(), by.x = "entrezgene_id", by.y = "EntrezID") |> distinct()
  write_csv(final_results, paste0(output_dir, "/", tt, "_KEGG_results_genes.csv"))

  # Data splits for EC, BRITE, Reaction, and Module visualizations
  ec_split <- final_results %>%
    filter(!is.na(EC)) %>%
    select(external_gene_name, EC) |>
    separate_rows(EC, sep = ";") %>%
    mutate(EC = trimws(EC)) %>%
    distinct()
  brite_split <- final_results %>%
    filter(!is.na(BRITE), BRITE != "NA") %>%
    select(external_gene_name, BRITE) |>
    separate_rows(BRITE, sep = ";") %>%
    mutate(BRITE = trimws(BRITE)) %>%
    distinct()
  reaction_split <- final_results %>%
    filter(!is.na(REACTION), REACTION != "NA") %>%
    select(external_gene_name, REACTION) |>
    separate_rows(REACTION, sep = ";") %>%
    mutate(REACTION = trimws(REACTION)) %>%
    distinct()
  module_split <- final_results %>%
    filter(!is.na(MODULE), MODULE != "NA") %>%
    select(external_gene_name, MODULE) |>
    separate_rows(MODULE, sep = ";") %>%
    mutate(MODULE = trimws(MODULE)) %>%
    distinct()

  # Generate combined PDF with plots
  pdf(file = paste0("plots/summary_plots_", tt, ".pdf"), height = 15, onefile = TRUE)

  # EC plot
  if (nrow(ec_split) > 0) {
    ec_count <- ec_split %>%
      group_by(EC) %>%
      summarise(count = n()) %>%
      arrange(desc(count))
    ec_plot <- ggplot(ec_count, aes(x = reorder(EC, count), y = count)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(x = "EC Categories", y = "Frequency", title = tt) +
      theme_minimal()
    print(ec_plot)
  }

  # BRITE plot
  if (nrow(brite_split) > 0) {
    brite_count <- brite_split %>%
      group_by(BRITE) %>%
      summarise(count = n()) %>%
      arrange(desc(count))
    brite_plot <- ggplot(brite_count, aes(x = reorder(BRITE, count), y = count)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(x = "BRITE Categories", y = "Frequency", title = tt) +
      theme_minimal()
    print(brite_plot)
  }

  # Reaction plot
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

  # Module plot
  if (nrow(module_split) > 0) {
    module_count <- module_split %>%
      group_by(MODULE) %>%
      summarise(count = n()) %>%
      arrange(desc(count))
    module_plot <- ggplot(module_count, aes(x = reorder(MODULE, count), y = count)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(x = "Module Categories", y = "Frequency", title = tt) +
      theme_minimal()
    print(module_plot)
  }

  dev.off()
}
