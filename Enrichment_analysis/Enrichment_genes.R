# Uncomment the following lines if 'BiocManager' is not already installed.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("clusterProfiler")

# Load required libraries
library(clusterProfiler) # For gene enrichment analysis
library(tidyverse) # For data manipulation and visualization
library("AnnotationDbi") # For working with annotation databases
library(org.Mm.eg.db) # Mouse gene annotation database
library(ReactomePA) # For Reactome pathway analysis

# Set the working directory for file operations
setwd("~/Documents/Research/Experiments/shotgun_hybrids/mapping_results/Results_genes/enrichment/")

# Read the GWAS results CSV file containing functional traits data
results_file <- read_csv("../../All_functional_traits_GWAS_results_SW.csv")

# Group genes by trait and save the results to a CSV file
genes_by_trait <- results_file %>%
  group_by(trait) %>%
  summarise(all_genes = paste0(unique(na.omit(list_total_genes)), collapse = ","))
write_csv(genes_by_trait, "../Genes_by_trait.csv")

# Group protein-coding genes by trait and save the results
coding_genes_by_trait <- results_file %>%
  group_by(trait) %>%
  summarise(all_genes = paste0(unique(na.omit(list_protein_coding_genes)), collapse = ","))
write_csv(coding_genes_by_trait, "../coding_genes_by_trait.csv")

# Loop through each trait and perform gene enrichment analyses
for (i in 1:nrow(genes_by_trait)) {
  kk <- NULL
  mkk <- NULL
  x <- NULL
  ego_mf <- NULL
  ego_bp <- NULL

  # Extract the trait and associated genes
  trait <- as.character(genes_by_trait[i, 1])
  geneset <- trimws(as.character(str_split(genes_by_trait[i, "all_genes"], ",", simplify = T)), "both")

  # Convert gene symbols to Entrez IDs
  genes_entrez <- as.numeric(na.omit(mapIds(org.Mm.eg.db, keys = geneset, keytype = "SYMBOL", column = "ENTREZID")))

  # KEGG pathway enrichment analysis
  kk <- enrichKEGG(
    gene = genes_entrez,
    organism = "mmu", # Mouse
    pvalueCutoff = 1, # Set significance threshold for p-value
    qvalueCutoff = 1 # Set significance threshold for q-value
  )
  if (!is.null(kk)) {
    kk <- setReadable(kk, org.Mm.eg.db, keyType = "ENTREZID")
    write_csv(as.data.frame(kk), paste0(trait, "_enrichmentKEGG_pathways_results.csv"))
  }

  # KEGG modules enrichment analysis
  mkk <- enrichMKEGG(
    gene = genes_entrez,
    organism = "mmu", # Mouse
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  if (!is.null(mkk)) {
    mkk <- setReadable(mkk, org.Mm.eg.db, keyType = "ENTREZID")
    write_csv(as.data.frame(mkk), paste0(trait, "_enrichmentKEGG_modules_results.csv"))
  }

  # Reactome pathway enrichment analysis
  x <- enrichPathway(gene = genes_entrez, pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE, organism = "mouse")
  if (!is.null(x)) {
    x <- setReadable(x, org.Mm.eg.db, keyType = "ENTREZID")
    write_csv(as.data.frame(x), paste0(trait, "_enrichment_Reactome_pathways_results.csv"))
  }

  # GO term enrichment for Biological Process (BP)
  ggo_bp <- groupGO(
    gene = as.character(genes_entrez),
    OrgDb = org.Mm.eg.db,
    ont = "BP", # Biological Process
    level = 4, # Restrict to level 4 in the GO hierarchy
    readable = TRUE
  )
  ggo_bp <- as.tibble(ggo_bp) %>% filter(Count > 0)
  if (!is.null(ggo_bp)) {
    write_csv(as.data.frame(ggo_bp), paste0(trait, "_all_GO_BP_results.csv"))
  }

  # GO term enrichment for Molecular Function (MF)
  ggo_mf <- groupGO(
    gene = as.character(genes_entrez),
    OrgDb = org.Mm.eg.db,
    ont = "MF", # Molecular Function
    level = 4,
    readable = TRUE
  )
  ggp_mf <- as.tibble(ggo_mf) %>% filter(Count > 0)
  if (!is.null(ggo_mf)) {
    write_csv(as.data.frame(ggo_mf), paste0(trait, "_all_GO_MF_results.csv"))
  }

  # Enrichment analysis for GO Molecular Function (MF)
  ego_mf <- enrichGO(
    gene = genes_entrez,
    OrgDb = org.Mm.eg.db,
    ont = "MF", # Molecular Function
    pAdjustMethod = "BH", # Benjamini-Hochberg correction
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
  if (!is.null(ego_mf)) {
    write_csv(as.data.frame(ego_mf), paste0(trait, "_enrichment_GO_MF_results.csv"))
  }

  # Enrichment analysis for GO Biological Process (BP)
  ego_bp <- enrichGO(
    gene = genes_entrez,
    OrgDb = org.Mm.eg.db,
    ont = "BP", # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
  if (!is.null(ego_bp)) {
    write_csv(as.data.frame(ego_bp), paste0(trait, "_enrichment_GO_BP_results.csv"))
  }
}

# Additional enrichment analysis for marker genes
marker_genes <- read_csv("../../markers_with_genes_functional_traits_SW.csv")
all_genes <- na.omit(unique(trimws(unlist(str_split(marker_genes$gene, pattern = " ")), which = "both")))

# Save all marker genes to a CSV file
write_csv(as.data.frame(all_genes[-4]), file = "../marker_genes.csv", col_names = F)

# Convert marker gene symbols to Entrez IDs
genes_entrez <- as.numeric(na.omit(mapIds(org.Mm.eg.db, keys = all_genes[-4], keytype = "SYMBOL", column = "ENTREZID")))

# Perform KEGG pathway enrichment for marker genes
kk <- enrichKEGG(
  gene = genes_entrez,
  organism = "mmu", # Mouse
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
kk <- setReadable(kk, org.Mm.eg.db, keyType = "ENTREZID")
write_csv(as.data.frame(kk@result), "closest_genes_enrichmentKEGG_pathways_results.csv")

# Perform KEGG module enrichment for marker genes
mkk <- enrichMKEGG(
  gene = genes_entrez,
  organism = "mmu", # Mouse
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
if (!is.null(mkk)) {
  mkk <- setReadable(mkk, org.Mm.eg.db, keyType = "ENTREZID")
  write_csv(as.data.frame(mkk@result), "closest_genes_enrichmentKEGG_modules_results.csv")
}

# Perform Reactome pathway enrichment for marker genes
x <- enrichPathway(gene = as.character(genes_entrez), pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE, organism = "mouse")
if (!is.null(x)) {
  x <- setReadable(x, org.Mm.eg.db, keyType = "ENTREZID")
  write_csv(as.data.frame(x@result), "closest_genes_enrichment_Reactome_pathways_results.csv")
}

# Perform GO term enrichment for Molecular Function (MF) for marker genes
ego_mf <- enrichGO(
  gene = genes_entrez,
  OrgDb = org.Mm.eg.db,
  ont = "MF", # Molecular Function
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)
if (!is.null(ego_mf)) {
  egoSimGO_mf <- clusterProfiler::simplify(ego_mf,
    cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang",
    semData = NULL
  )
  write_csv(as.data.frame(egoSimGO_mf@result), "closest_genes_enrichment_GO_MF_results.csv")
}

# Perform GO term enrichment for Biological Process (BP) for marker genes
ego_bp <- enrichGO(
  gene = genes_entrez,
  OrgDb = org.Mm.eg.db,
  ont = "BP", # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)
if (!is.null(ego_bp)) {
  egoSimGO_bp <- clusterProfiler::simplify(ego_bp,
    cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang",
    semData = NULL
  )
  write_csv(as.data.frame(egoSimGO_bp), "closest_genes_enrichment_GO_BP_results.csv")
}

# Perform GO term enrichment for Cellular Component (CC) for marker genes
ego_cc <- enrichGO(
  gene = genes_entrez,
  OrgDb = org.Mm.eg.db,
  ont = "CC", # Cellular Component
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)
if (!is.null(ego_cc)) {
  egoSimGO_cc <- clusterProfiler::simplify(ego_cc,
    cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang",
    semData = NULL
  )
  write_csv(as.data.frame(egoSimGO_cc), "closest_genes_enrichment_GO_CC_results.csv")
}

# Visualize the enrichment results for Cellular Component (CC), Molecular Function (MF), and Biological Process (BP)
p_cc <- barplot(egoSimGO_cc, showCategory = 20) +
  ggtitle("GO-CC") +
  xlab("Enriched terms") + ylab("Count") + scale_fill_viridis_c(trans = "log10")
ggsave("GO-CC.pdf", height = 10, width = 12, plot = p_cc)

p_mf <- barplot(egoSimGO_mf, showCategory = 20) +
  ggtitle("GO-MF") +
  xlab("Enriched terms") + ylab("Count") + scale_fill_viridis_c(trans = "log10")
ggsave("GO-MF.pdf", height = 10, width = 12, plot = p_mf)

p_bp <- barplot(egoSimGO_bp, showCategory = 20) +
  ggtitle("GO-BP") +
  xlab("Enriched terms") + ylab("Count") + scale_fill_viridis_c(trans = "log10")
ggsave("GO-BP.pdf", height = 10, width = 12, plot = p_bp)
