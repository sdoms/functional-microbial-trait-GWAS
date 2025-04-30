# This script will find overlap between sig regions from different traits
# and overlap with previous studies
setwd("~/Documents/Research/Experiments/shotgun_hybrids/mapping_results/overlap_other_studies")
library(readxl)
library(tidyverse)

loadRData <- function(fileName) {
  # loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

Meff.threshold <- 0.05 / 4603

study.wide.meff <- Meff.threshold / 140.9546
# Format our results into bed files
result_table <- read_csv("../All_functional_traits_GWAS_results_SW.csv")[, -1]

results_table_SW <- result_table |> filter(peak.Pval < study.wide.meff)

results.bed <- results_table_SW[, c("chr", "start.pos", "stop.pos")]
results.bed$chr <- paste0("chr", results.bed$chr)
write.table(results.bed, file = "all_combined_traits.bed", row.names = F, col.names = F, sep = "\t", quote = F, options(scipen = 10))
# Used bedtools to sort and merge overlapping regions to get combined sig regions all_combined_traits_merged.bed

# remove duplicate regions to get file to use for intersect
results.unique.bed <- unique(results.bed)
write.table(results.unique.bed, file = "all_combined_traits_unique_SW.bed", row.names = F, col.names = F, sep = "\t", quote = F, options(scipen = 10))

# Make bed file for regions from previous studies
previous_table <- loadRData("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/all_DNA_RNA_22022022_SW_pval_corr.Rdata")
previous_table <- previous_table[, -c(1, 12:23)]
previous_table <- previous_table %>%
  group_by(trait, chr, start.LD.pos, stop.LD.pos) %>%
  mutate(P.type = paste(P.type, collapse = ", "), trait = paste0(dna.rna, "_", tax_level, "_", trait)) %>%
  distinct()

previous.bed <- previous_table[, c("chr", "start.LD.pos", "stop.LD.pos")]
previous.bed$chr <- paste0("chr", previous.bed$chr)
previous.unique.bed <- unique(previous.bed)
write.table(previous.unique.bed, file = "sQTLprevStudies_SW.bed", row.names = F, col.names = F, sep = "\t", quote = F, options(scipen = 10))
# Use bedtools to get intersecting region pairs
command <- "/opt/homebrew/bin/bedtools intersect -a all_combined_traits_unique_SW.bed -b sQTLprevStudies_SW.bed -wo > results.previous.matchList.20241024.SW.bed"
cat(command, "\n")
try(system(command))

# test if there is significant overlap

system("python3 /Users/doms/Documents/Research/Software/poverlap/poverlap.py poverlap --a sQTLprevStudies_SW.bed --b all_combined_traits_unique_SW.bed -g ~/Documents/Research/Experiments/shotgun_hybrids/mouse.mm10.genome --n 9999")
# {"'wc -l'": {"observed": 5328.0, "shuffle_cmd": "bedtools intersect -wo -a sQTLprevStudies_SW.bed -b <(bedtools shuffle -allowBeyondChromEnd   -i all_combined_traits_unique_SW.bed -g /Users/doms/Documents/Research/Experiments/shotgun_hybrids/mouse.mm10.genome )", "metric": "'wc -l'", "simulated mean metric": 3658.5286528652864, "simulated_p": 0.0001}}
match_list <- read.table("results.previous.matchList.20241024.SW.bed")
colnames(match_list) <- c("Chr.res", "Start.res", "Stop.res", "Chr.prev", "Start.prev", "Stop.prev", "overlap")
match_list$interval.res <- paste(match_list$Chr.res, match_list$Start.res, match_list$Stop.res, sep = "_")
match_list$interval.prev <- paste(match_list$Chr.prev, match_list$Start.prev, match_list$Stop.prev, sep = "_")
match_list$interval.res <- gsub("chr", "", match_list$interval.res)
match_list$interval.prev <- gsub("chr", "", match_list$interval.prev)


previous_table$interval <- paste(previous_table$chr, previous_table$start.LD.pos, previous_table$stop.LD.pos, sep = "_")

output_table <- tibble()
for (row in 1:nrow(result_table)) {
  interval <- paste(result_table[row, "chr"], result_table[row, "start.pos"], result_table[row, "stop.pos"], sep = "_")
  matching.intervals <- match_list[match_list$interval.res == interval, "interval.prev"]
  overlap <- previous_table %>% filter(interval %in% matching.intervals)
  output_table[row, 1:5] <- result_table[row, c(1, 2, 3, 5, 6)]
  if (!is_empty(overlap$trait)) {
    output_table[row, 6] <- paste(unique(overlap$trait), collapse = ", ")
    output_table[row, 7] <- length(unique(overlap$trait))
  } else {
    output_table[row, 6] <- NA
    output_table[row, 7] <- 0
  }
}
colnames(output_table) <- c("trait", "chr_bin", "chr", "start.pos", "stop.pos", "overlap.phenotype.taxa", "number.of.overlaps.taxa")
write_csv(output_table, "overlap_with_sig_regions_16S_GWAS_20241121_SW.csv")

### regions within our study

# Use bedtools to get intersecting region pairs
command <- "/opt/homebrew/bin/bedtools intersect -a all_combined_traits_unique_SW.bed -b all_combined_traits_unique_SW.bed -wo > results.self.matchList.20241024.SW.bed"
cat(command, "\n")
try(system(command))

match_list_results <- read.table("results.self.matchList.20241024.SW.bed")
colnames(match_list_results) <- c("Chr.res1", "Start.res1", "Stop.res1", "Chr.res2", "Start.res2", "Stop.res2", "overlap")
match_list_results$interval.res1 <- paste(match_list_results$Chr.res1, match_list_results$Start.res1, match_list_results$Stop.res1, sep = "_")
match_list_results$interval.res2 <- paste(match_list_results$Chr.res2, match_list_results$Start.res2, match_list_results$Stop.res2, sep = "_")
match_list_results$interval.res1 <- gsub("chr", "", match_list_results$interval.res1)
match_list_results$interval.res2 <- gsub("chr", "", match_list_results$interval.res2)

result_table$interval <- paste(result_table$chr, result_table$start.pos, result_table$stop.pos, sep = "_")

homoreg_table <- tibble()
for (row in 1:nrow(result_table)) {
  interval <- result_table$interval[row]
  pheno <- result_table$trait[row]
  matching.intervals <- match_list_results[match_list_results$interval.res1 == interval, "interval.res2"]

  overlap <- result_table %>% filter(interval %in% matching.intervals & trait != pheno)
  homoreg_table[row, 1:5] <- result_table[row, c(1, 2, 3, 5, 6)]
  if (!is_empty(overlap$trait)) {
    homoreg_table[row, 6] <- paste(overlap$trait, overlap$chr_bin, collapse = ", ")
    homoreg_table[row, 7] <- length(overlap$trait)
  } else {
    homoreg_table[row, 6] <- NA
    homoreg_table[row, 7] <- 0
  }
}
colnames(homoreg_table) <- c("trait", "chr_bin", "chr", "start.pos", "stop.pos", "overlap.phenotype.GMMs", "number.of.overlaps.GMMs")

all_overlaps <- cbind(output_table, homoreg_table[, 6:7])
all_overlaps_info <- cbind(all_overlaps, result_table)

write_csv(all_overlaps, "overlaps_other_studies_and_own_study20241024_SW.csv")
write_csv(all_overlaps_info, "overlaps_other_studies_and_own_study_with_info20241024_SW.csv")

#### mills genes ####
library(biomaRt)
mills_genes <- read_excel("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/overexpressed_protein_GF_con.xlsx", sheet = "ALL_DF", col_names = "genes")
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "www") # can take long



attributes.1 <- c("chromosome_name", "start_position", "end_position", "strand", "mgi_symbol")
results.1 <- getBM(attributes = attributes.1, filters = c("mgi_symbol"), values = mills_genes$genes, mart = gene.ensembl)

bed_results <- results.1 %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  dplyr::select(chr, start_position, end_position)
write.table(bed_results, "mills_genes_bed.txt", row.names = F, col.names = F, sep = "\t", quote = F, options(scipen = 10))

system("python3 /Users/doms/Documents/Research/Software/poverlap/poverlap.py poverlap --a mills_genes_bed.txt --b all_combined_traits_unique_SW.bed -g ~/Documents/Research/Experiments/shotgun_hybrids/mouse.mm10.genome --n 9999")
# output:
# {"'wc -l'": {"observed": 312.0, "shuffle_cmd": "bedtools intersect -wo -a mills_genes_bed.txt -b <(bedtools shuffle -allowBeyondChromEnd   -i all_combined_traits_unique_SW.bed -g /Users/doms/Documents/Research/Experiments/shotgun_hybrids/mouse.mm10.genome )", "metric": "'wc -l'", "simulated mean metric": 244.87988798879888, "simulated_p": 0.0168}}
# Use bedtools to get intersecting region pairs
command <- "/opt/homebrew/bin/bedtools intersect -a all_combined_traits_unique_SW.bed -b mills_genes_bed.txt -wo > results.mills.matchList.SW.20241024.bed"
cat(command, "\n")
try(system(command))


match_list_mills <- read.table("results.mills.matchList.SW.20241024.bed")
colnames(match_list_mills) <- c("Chr.res", "Start.res", "Stop.res", "Chr.prev", "Start.prev", "Stop.prev", "overlap")
match_list_mills$interval.res <- paste(match_list_mills$Chr.res, match_list_mills$Start.res, match_list_mills$Stop.res, sep = "_")
match_list_mills$interval.prev <- paste(match_list_mills$Chr.prev, match_list_mills$Start.prev, match_list_mills$Stop.prev, sep = "_")
match_list_mills$interval.res <- gsub("chr", "", match_list_mills$interval.res)
match_list_mills$interval.prev <- gsub("chr", "", match_list_mills$interval.prev)


results.1$interval <- paste(results.1$chr, results.1$start_position, results.1$end_position, sep = "_")

output_table_mills <- tibble()
for (row in 1:nrow(result_table)) {
  interval <- paste(result_table[row, "chr"], result_table[row, "start.pos"], result_table[row, "stop.pos"], sep = "_")
  matching.intervals <- match_list_mills[match_list_mills$interval.res == interval, "interval.prev"]
  overlap <- results.1 %>% filter(interval %in% matching.intervals)
  output_table_mills[row, 1:5] <- result_table[row, c(1, 2, 3, 5, 6)]
  if (!is_empty(overlap$mgi_symbol)) {
    output_table_mills[row, 6] <- paste(unique(overlap$mgi_symbol), collapse = ", ")
    output_table_mills[row, 7] <- length(unique(overlap$mgi_symbol))
  } else {
    output_table_mills[row, 6] <- NA
    output_table_mills[row, 7] <- 0
  }
}
colnames(output_table_mills) <- c("trait", "chr_bin", "chr", "start.pos", "stop.pos", "overlap.genes.mills", "number.of.overlaps.mills")
output_total <- cbind(output_table, output_table_mills)
write_csv(output_total, "overlap_with_sig_regions_16S_GWAS_mills_20241024_SW.csv")

write.table(as.data.frame(na.omit(unique(unlist(str_split(output_table_mills$overlap.genes.mills, ", "))))), "genes_in_mill_SW.txt",
  col.names = F, row.names = F, quote = F
)
