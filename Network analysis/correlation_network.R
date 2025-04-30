# In this script I want to make a correlation network
# 10/10/24: updated script to include the taxonomic 16S traits and results
# 16/10/24: study-wide Meff threshold, changed gene names so that it is closest and removed abs from correlation threshold
# 22/11/24: make combined edge weight by overlaps between traits * correlation
# load corralation network made in prep_for_matSpD.R script
setwd("~/Documents/Research/Experiments/shotgun_hybrids/mapping_results/correlation_overlap/correlation_network/")
# # cor_mat <- read.csv("correlation.matrix_functional_traits_rownames.csv", row.names = 1) # without the taxa
# cor_mat <- read.csv("correlation.matrix_functional_traits_and_taxa_16S_rownames.csv", row.names = 1)

cor_edges <- read.csv("../overlaps_with_correlation.csv")
# 1. Trait network construction
# Load necessary libraries
library(igraph)
library(ggraph)
library(tidygraph)
library(tidyverse)

# Assuming `cor_matrix` is your correlation matrix
# Set threshold to include only correlations above a certain value
# threshold <- 0.8
# cor_edges <- cor_mat %>%
#   as.data.frame() %>%
#   rownames_to_column("trait1") %>%
#   gather(key = "trait2", value = "correlation", -trait1) %>%
#   filter(correlation> threshold) |>
#   filter(trait1 != trait2) |>
#   mutate(type = "correlation")

cor_edges$edge_weight_perc <- (cor_edges$edge_weight / max(cor_edges$edge_weight))
# Create a graph object
threshold <- 0.4
cor_edges <- cor_edges |> filter(edge_weight_perc > threshold)
g <- graph_from_data_frame(cor_edges, directed = FALSE)

# Plot the network with correlation as edge width
ggraph(g, layout = "fr") + # "fr" = Fruchterman-Reingold layout
  geom_edge_link(aes(width = edge_weight_perc), color = "blue", alpha = 0.7) +
  geom_node_point(size = 5, color = "red") +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void() # Cleaner plot
ggsave("network_edge_weight_perc_cor_overlap_0.4.pdf", width = 10, height = 10)
# 2. Incorporating GWAS results
# Example GWAS results dataframe
# gwas_results <- data.frame(trait = c(...), SNP = c(...), p_value = c(...), z_score = c(...))
load("../../all_results_GWAS_functional_traits_AND_taxonomic_16S_SW.RData")

#

# Add SNPs as additional nodes and link them to traits
gwas_edges <- KO_taxa %>%
  select(trait, peak.snp, peak.Pval, closest_gene, closest_gene_type) %>%
  mutate(weight = -log10(peak.Pval), type = "GWAS") |> # Convert p-value to more interpretable scale
  mutate(weight01 = weight / max(weight))

snps_gene_annotation <- read_delim("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/snps_gene_annotation.csv", delim = ";")
snps_gene_annotation <- snps_gene_annotation |>
  separate(gene, into = c("gene1", "gene2")) |>
  filter(gene1 != "" | gene2 != "") |>
  select(marker, gene1, gene2) |>
  mutate(
    gene1 = na_if(gene1, ""), # Convert "" to NA in gene1
    gene2 = na_if(gene2, "")
  ) %>% # Convert "" to NA in gene2
  mutate(gene = coalesce(gene1, gene2))

gwas_edges_new <- left_join(gwas_edges, snps_gene_annotation, by = c("peak.snp" = "marker")) |>
  distinct_all() |>
  mutate(gene = coalesce(closest_gene, gene))

# Add SNPs to the existing network
nodes <- unique(c(cor_edges$trait1, cor_edges$trait2, gwas_edges_new$trait, gwas_edges_new$gene))
edges <- bind_rows(
  cor_edges %>% rename(weight = edge_weight_perc), # Trait correlations
  gwas_edges_new %>% rename(trait1 = trait, trait2 = gene) # GWAS results
)

g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

# Visualize with SNPs
ggraph(g, layout = "fr") +
  geom_edge_link(aes(edge_alpha = abs(weight), color = weight > 0)) +
  geom_node_point(aes(
    size = ifelse(name %in% gwas_edges$SNP, 5, 3),
    color = ifelse(name %in% gwas_edges$SNP, "blue", "black")
  )) +
  geom_node_text(aes(label = name), vjust = 1.5, size = 3) +
  scale_edge_color_manual(values = c("lightgray", "red")) + # Separate colors for positive/negative weights
  theme_void()

# 3. Gene interactions
gene_list_string <- data.frame(unique(gwas_edges_new$gene))
write.table(gene_list_string, file = "../../STRING_interactions_taxa_KO/gene_list_20241122.txt", quote = F, row.names = F, col.names = F)
gene_interactions <- read_table("../../STRING_interactions_taxa_KO/string_interactions_short_medium_confidence_20241122.tsv")
# Convert gene interactions to edges and combine with trait-gene network
gene_edges <- gene_interactions %>%
  rename(trait1 = `#node1`, trait2 = node2, weight = combined_score) |>
  mutate(type = "gene2gene")

combined_all_edges <- bind_rows(edges, gene_edges)

# Create final network with traits, SNPs, genes, and interactions
final_network <- graph_from_data_frame(d = combined_all_edges, directed = FALSE)

# Visualize the final integrated network
ggraph(final_network, layout = "fr") +
  geom_edge_link(aes(width = ifelse(is.na(weight), 0.5, abs(weight))), color = "gray") +
  geom_node_point(aes(color = ifelse(name %in% gwas_edges$closest_gene, "blue",
    ifelse(name %in% gene_edges$trait1 | name %in% gene_edges$trait2, "green", "red")
  )), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void()



# Install RCy3 package if not already installed

# BiocManager::install("RCy3")

# Load the package
library(RCy3)

# Convert igraph network to data frames for nodes and edges
nodes_igraph <- data.frame(id = V(final_network)$name) # Create node data frame
edges_igraph <- igraph::as_data_frame(final_network, what = "edges") # Create edge data frame
# Initialize connection to Cytoscape
cytoscapePing()
# Create a network in Cytoscape from the node and edge data frames

ig <- graph_from_data_frame(combined_all_edges, directed = FALSE, vertices = nodes_igraph)

createNetworkFromIgraph(ig, "Taxa_KO")
# Export to GraphML
write_graph(ig, file = "network_GWAS_functions_taxa_SW_20241122.graphml", format = "graphml")

# Export to SIF (simple interaction format)
write_graph(ig, file = "network_GWAS_functions_taxa_SW_20241122.sif", format = "ncol")

#
# Convert igraph network to data frames for nodes and edges
nodes_igraph <- data.frame(id = V(ig)$name) # Create node data frame
nodes_igraph <- nodes_igraph |> mutate("trait_cat" = ifelse(startsWith(id, "DNA_"), "DNA_taxon",
  ifelse(startsWith(id, "RNA_"), "RNA_taxon",
    ifelse(startsWith(id, "BRITE_"), "BRITE",
      ifelse(startsWith(id, "EC_level"), "EC_level",
        ifelse(startsWith(id, "MF"), "GMM", ifelse(startsWith(id, "MGB"), "GBM", "gene"))
      )
    )
  )
))
edges_igraph <- igraph::as_data_frame(ig, what = "edges") # Create edge data frame
# Export node table to CSV
write.csv(nodes_igraph, file = "../../STRING_interactions_taxa_KO/nodes_graph_SW_20241122.csv", row.names = FALSE)

#
# Optional: Add additional edge attributes if they exist
edges_igraph$weight <- E(ig)$weight
edges_igraph$correlation <- E(ig)$correlation

# Export edge table to CSV
write.csv(edges_igraph, file = "../../STRING_interactions_taxa_KO/edges_graph_SW_20241122.csv", row.names = FALSE)

