setwd("~/Documents/Research/Experiments/shotgun_hybrids/mapping_results/GMMs/")


loadRData <- function(fileName) {
  # loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

library(tidyverse)
library(stringr)
library(patchwork)
library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # for annotation
library(org.Mm.eg.db)
library(PurrfectPlots)

threshold <- 0.05 / 30401

musdom <- read.csv("~/Documents/Research/Experiments/Final_QTL_mapping/Results/allele_freq_consensus_mus_dom.csv", sep = ",")
musdom$dd_frew <- NULL
load("~/Documents/Research/Experiments/Final_QTL_mapping/Cleaning_snps/consensus_F0_all.Rdata")
parents <- as.data.frame(reduced_geno)
parents$X <- rownames(parents)
musdom <- merge(musdom, parents, by = "X", )

outputdir <- "all_genes_snps/"
dir.create(outputdir, showWarnings = F)

#####################################################################
##                        Summary table 1:                         ##
##  a table that lists all significant SNPs                        ##
#####################################################################

trait_list <- read.table("../../GMM_abundance.tsv")[, 1:2]

tot_files <- tibble()
for (i in 1:nrow(trait_list)) {
  gws <- trait_list$Module[i]

  output.name <- paste0(trait_list$Module[i], "_", gsub(" ", "_", gsub(pattern = "[^[:alnum:] ]", "_", trait_list$Description[i])))

  if (file.exists(paste0("sig.summaries/", gws, "/", output.name, "_gwscan_sigSNPs.Rdata"))) {
    infile <- loadRData(paste0("sig.summaries/", gws, "/", output.name, "_gwscan_sigSNPs.Rdata"))
    if (nrow(infile) > 0) {
      infile$trait <- gws
    }
    tot_files <- rbind(tot_files, infile)
  }
}

write_csv(tot_files, file = paste0(outputdir, "all_sig_snps_GW.csv"))
save(tot_files, file = paste0(outputdir, "all_sig_snps_GW.Rdata"))



################################################################################################
##                                      Summary table 2:                                      ##
##  a table that lists all uniaue significant SNPs for all P values for all taxonomic levels  ##
##              with a count of the times it was significant in different traits              ##
################################################################################################


all_sig_snps <- tot_files

ex <- all_sig_snps %>%
  left_join(musdom, by = c("marker" = "X")) %>%
  dplyr::group_by(marker) %>%
  add_count(marker, name = "count") %>%
  dplyr::mutate(all_traits = paste(trait, collapse = " | ")) %>%
  ungroup() %>%
  mutate(chr = as.numeric(as.character(chr)), pos = as.numeric(as.character(pos))) %>%
  dplyr::select(marker, chr, pos,
    A = A2, B = A1, mus, dom, AA, AB, BB, count, all_traits,
    FSP3, FSP5, HAP1, HAP2, HOP3, HOP6, TUP3, TUP4
  ) %>%
  distinct() %>%
  arrange(chr, pos) %>%
  dplyr::mutate(chr = replace_na(as.character(chr), "X"))
write_csv(ex, paste0(outputdir, "summary_table_markers_count_GW.csv"))



################################################################################################
##                                      Summary table 3:                                      ##
##  a table that lists all unique significant SNPs for all P values for all taxonomic levels  ##
##              with a count of the times it was significant in different traits              ##
##                            with the closest genes to the marker                            ##
################################################################################################

# Add closest genes to al snps
library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # for annotation
library(org.Mm.eg.db)

input <- ex %>%
  drop_na(chr) %>%
  dplyr::mutate(chr = paste0("chr", chr)) %>%
  dplyr::select(rsid = marker, chr, pos)
# input <- input[1:3,]
input
#         rsid  chr       pos
# 1  rs3753344 chr1   1142150
# 2 rs12191877 chr6  31252925
# 3   rs881375 chr9 123652898

target <- with(
  input,
  GRanges(
    seqnames = Rle(chr),
    ranges = IRanges(pos, end = pos, names = rsid),
    strand = Rle(strand("*"))
  )
)
target







loc <- locateVariants(target, TxDb.Mmusculus.UCSC.mm10.knownGene, AllVariants())
loc
names(loc) <- NULL

p_ids <- unlist(loc$PRECEDEID, use.names = FALSE)
exp_ranges <- rep(loc, elementNROWS(loc$PRECEDEID))
p_dist <- GenomicRanges::distance(exp_ranges, TxDb.Mmusculus.UCSC.mm10.knownGene, id = p_ids, type = "gene")
head(p_dist)
exp_ranges$PRECEDE_DIST <- p_dist
exp_ranges

## Collapsed view of ranges, gene id and distance:
loc$PRECEDE_DIST <- relist(p_dist, loc$PRECEDEID)
loc


f_ids <- unlist(loc$FOLLOWID, use.names = FALSE)
exp_ranges_f <- rep(loc, elementNROWS(loc$FOLLOWID))

f_dist <- GenomicRanges::distance(exp_ranges_f, TxDb.Mmusculus.UCSC.mm10.knownGene, id = f_ids, type = "gene")
head(f_dist)
exp_ranges_f$FOLLOW_DIST <- f_dist
exp_ranges_f

## Collapsed view of ranges, gene id and distance:
loc$FOLLOW_DIST <- relist(f_dist, loc$FOLLOWID)
loc

# find closest
# find closest
loc$PRECEDEID_place <- lapply(loc$PRECEDE_DIST, function(x) which.min(x))
loc$PRECEDEID_place <- as.numeric(as.character(loc$PRECEDEID_place))

# loc$PRECEDEID_closest <- lapply(loc$PRECEDEID, function(x) x[unlist(loc$PRECEDEID_place)])

loc$FOLLOWID_place <- lapply(loc$FOLLOW_DIST, function(x) which.min(x))
loc$FOLLOWID_place <- as.numeric(as.character(loc$FOLLOWID_place))

# out$FOLLOWID_closest <- apply(out, 1, function(x) dplyr::nth(x[14],unlist(x[20]), default=NULL))


out <- as.data.frame(loc)


out$names <- names(target)[out$QUERYID]

precede_ID <- enframe(as.list(out$PRECEDEID))[, 2]
colnames(precede_ID) <- "ID"
precede_ID$place <- out$PRECEDEID_place
for (x in 1:nrow(precede_ID)) {
  n <- as.numeric(precede_ID[x, "place"])
  precede_ID[x, "closest"] <- unlist(precede_ID[x, 1])[n]
}

follow_ID <- enframe(as.list(out$FOLLOWID))[, 2]
colnames(follow_ID) <- "ID"
follow_ID$place <- out$FOLLOWID_place
for (x in 1:nrow(follow_ID)) {
  n <- as.numeric(follow_ID[x, "place"])
  follow_ID[x, "closest"] <- unlist(follow_ID[x, 1])[n]
}


out <- out[, c("names", "seqnames", "start", "end", "LOCATION", "GENEID")]
out <- cbind(out, precede_ID$closest)
out <- cbind(out, follow_ID$closest)

out <- unique(out)
out

colnames(out) <- c("names", "seqnames", "start", "end", "LOCATION", "GENEID", "PRECEDEID", "FOLLOWID")
Symbol2id <- as.list(org.Mm.egSYMBOL2EG)
id2Symbol <- rep(names(Symbol2id), sapply(Symbol2id, length))
names(id2Symbol) <- unlist(Symbol2id)

x <- unique(with(out, c(levels(as.factor(GENEID)), levels(as.factor(PRECEDEID)), levels(as.factor(FOLLOWID)))))
table(x %in% names(id2Symbol)) # good, all found

out$GENESYMBOL <- id2Symbol[as.character(out$GENEID)]
out$PRECEDESYMBOL <- id2Symbol[as.character(out$PRECEDEID)]
out$FOLLOWSYMBOL <- id2Symbol[as.character(out$FOLLOWID)]
out
final_out <- merge(ex, out, by.x = "marker", by.y = "names", all.x = T)

ex_dist <- final_out %>%
  group_by(marker) %>%
  dplyr::mutate(locations = paste(LOCATION, collapse = " | "), all_genes = paste(GENESYMBOL, collapse = " | "), all_preceding_genes = paste(PRECEDESYMBOL, collapse = " | "), all_following_genes = paste(FOLLOWSYMBOL, collapse = " | ")) %>%
  distinct(marker, locations, all_genes, all_preceding_genes, all_following_genes) %>%
  left_join(final_out[, c(1:20)], by = "marker") %>%
  distinct() %>%
  dplyr::mutate(chr = as.numeric(as.character(chr)), pos = as.numeric(as.character(pos))) %>%
  arrange(chr, pos) %>%
  dplyr::mutate(chr = replace_na(as.character(chr), "X"))

write_csv(ex_dist, paste0(outputdir, "markers_with_genes_GW.csv"))




all_final <- read_csv(paste0(outputdir, "markers_with_genes_GW.csv"))




################################################################################################
##                                      Summary table 4:                                      ##
##  a table that calculates the amount of significant regions per trait                       ##
################################################################################################
library(ggsci)
load("/Users/doms/Documents/Research/Experiments/shotgun_hybrids/mapping_results/GMMs/sig.summaries/GMMs_allresults_GW.Rdata")


count_gmm <- all_files %>%
  mutate(P.type = sapply(P.type, function(x) {
    paste(sort(strsplit(x, ", ")[[1]]), collapse = ", ")
  })) |>
  group_by(trait, P.type) %>%
  dplyr::summarise(n = n()) %>%
  ungroup()


ggplot(count_gmm, aes(x = reorder(trait, n, sum), y = n, fill = P.type)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  theme_minimal() +
  scale_fill_purrs("disney") +
  labs(x = "", y = "Number of siginificant loci", fill = "P value")
ggsave(paste0(outputdir, "significant_loci_per_trait_GW.pdf"))
ggsave(paste0(outputdir, "significant_loci_per_trait_GW.png"))

### SW ####
trait_list <- read.table("~/Documents/Research/Experiments/shotgun_hybrids/GMM_abundance_F2_core.tsv")
KO <- "GMM"
tot_files <- tibble()
for (i in 1:nrow(trait_list)) {
  gws <- trait_list$Module[i]

  output.name <- paste0(trait_list$Module[i], "_", gsub(" ", "_", gsub(pattern = "[^[:alnum:] ]", "_", trait_list$Description[i])))

  if (file.exists(paste0("sig.summaries/", gws, "/", output.name, "_gwscan_sigSNPs.Rdata"))) {
    infile <- loadRData(paste0("sig.summaries/", gws, "/", output.name, "_gwscan_sigSNPs.Rdata"))
    if (nrow(infile) > 0) {
      infile$trait <- gws
    }
    tot_files <- rbind(tot_files, infile)
  }
}

tot_files_SW <- tot_files |> filter(P.sig.Meff_SW == "Y" | add.P.sig.Meff_SW == "Y" | dom.P.sig.Meff_SW == "Y")

write_csv(tot_files_SW, file = paste0(outputdir, "all_sig_snps_", KO, "_SW.csv"))
save(tot_files_SW, file = paste0(outputdir, "all_sig_snps_", KO, "_SW.Rdata"))


################################################################################################
##                                      Summary table 2:                                      ##
##  a table that lists all uniaue significant SNPs for all P values for all taxonomic levels  ##
##              with a count of the times it was significant in different traits              ##
################################################################################################


all_sig_snps <- tot_files_SW

ex <- all_sig_snps %>%
  left_join(musdom, by = c("marker" = "X")) %>%
  dplyr::group_by(marker) %>%
  add_count(marker, name = "count") %>%
  dplyr::mutate(all_traits = paste(trait, collapse = " | ")) %>%
  ungroup() %>%
  mutate(chr = as.numeric(as.character(chr)), pos = as.numeric(as.character(pos))) %>%
  dplyr::select(marker, chr, pos,
    A = A2, B = A1, mus, dom, AA, AB, BB, count, all_traits,
    FSP3, FSP5, HAP1, HAP2, HOP3, HOP6, TUP3, TUP4
  ) %>%
  distinct() %>%
  arrange(chr, pos) %>%
  dplyr::mutate(chr = replace_na(as.character(chr), "X"))
write_csv(ex, paste0(outputdir, "summary_table_markers_count_", KO, "_SW.csv"))



################################################################################################
##                                      Summary table 3:                                      ##
##  a table that lists all unique significant SNPs for all P values for all taxonomic levels  ##
##              with a count of the times it was significant in different traits              ##
##                            with the closest genes to the marker                            ##
################################################################################################
snps_gene_annotation <- read_delim("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/snps_gene_annotation.csv", delim = ";")

final_out <- left_join(ex, snps_gene_annotation |> dplyr::select(-chr, -pos), by = "marker")

ex_dist <- final_out %>%
  dplyr::mutate(chr = as.numeric(as.character(chr)), pos = as.numeric(as.character(pos)), gene = trimws(gene)) %>%
  arrange(chr, pos) %>%
  dplyr::mutate(chr = replace_na(as.character(chr), "X"))

write_csv(ex_dist, paste0(outputdir, "markers_with_genes_", KO, "_SW.csv"))





all_final <- read_csv(paste0(outputdir, "markers_with_genes_", KO, "_SW.csv"))




################################################################################################
##                                      Summary table 4:                                      ##
##  a table that calculates the amount of significant regions per trait                       ##
################################################################################################
library(ggsci)
load("/Users/doms/Documents/Research/Experiments/shotgun_hybrids/mapping_results/GMMs/sig.summaries/GMMs_allresults_SW.Rdata")


count_gmm <- all_files %>%
  mutate(P.type = sapply(P.type, function(x) {
    paste(sort(strsplit(x, ", ")[[1]]), collapse = ", ")
  })) |>
  group_by(trait, P.type) %>%
  dplyr::summarise(n = n()) %>%
  ungroup()


ggplot(count_gmm, aes(x = reorder(trait, n, sum), y = n, fill = P.type)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  theme_minimal() +
  scale_fill_purrs("disney") +
  labs(x = "", y = "Number of siginificant loci", fill = "P value")
ggsave(paste0(outputdir, "significant_loci_per_trait_GMM_SW.pdf"), width = 12, height = 10)
ggsave(paste0(outputdir, "significant_loci_per_trait_GMM_SW.png"), width = 12, height = 10)
