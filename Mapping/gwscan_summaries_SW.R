require(biomaRt)
require(stringr)
require(forcats)
# require(RColorBrewer)
require(viridis)
require(plyr)
require(patchwork)
require(tidyverse)
library(VariantAnnotation)

KO <- "EC_level" # change this to EC_level, ...
plink <- "/Users/doms/Documents/Research/Software/plink/plink"

# load an Rdata file an rename it
loadRData <- function(fileName) {
  # loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# Biomarts - set for mm10
snp.ensembl <- useEnsembl(biomart = "snp", dataset = "mmusculus_snp", version = 102)
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", version = 102)

# load all snp info for expanding intervals by including snps in LD - see below
load("~/Documents/Research/Experiments/Final_QTL_mapping/Scripts/Data/genotypes/preLDpruneSNPs.Rdata")
preLDpruneSNPs$pos <- preLDpruneSNPs$pos * 1e6

# INPUTS##
# list of traits mapped - "Trait" matches folder name,"TraitFileName" matches .rds prefixes (differ for covars)

# trait_list <- read.csv("~/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/EC_level_1_core.csv")
# trait_list <- read.csv("~/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/EC_level_2_core.csv")
trait_list <- read.csv("~/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/EC_level_3_core.csv")
# trait_list <- read.csv("~/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/BRITE_categories_all_core_median_50.csv")

# SW threshold
# correlation matrix
# cor_df <- trait_list[, -c(1:2)]
# rownames(cor_df)<- trait_list$Module
# correlation.matrix <- cor(t(cor_df), method = "spearman", use = "complete.obs")
# write.table(correlation.matrix, file = "correlation.matrix.GBM.txt")
# source("../../scripts/mapping/matSpDlite.R")

# X chromosome results?
Xtrue <- F # nolint: object_name_linter.
# folder that has a results folder for each trait
# results.folder <- "~/Documents/Research/Experiments/shotgun_hybrids/mapping_results/EC_level_1_core_out/"
# results.folder <- "~/Documents/Research/Experiments/shotgun_hybrids/mapping_results/EC_level_2_core_out/"
results.folder <- "~/Documents/Research/Experiments/shotgun_hybrids/mapping_results/EC_level_3_core_out/"
# results.folder <- "~/Documents/Research/Experiments/shotgun_hybrids/mapping_results/BRITE/"

# make output folder for sig summaries
output.folder <- paste0(results.folder, "sig.summaries")
dir.create(output.folder, showWarnings = FALSE)
setwd(output.folder)

for (i in 1:nrow(trait_list)) {
  gws <- trait_list$Identifier[i]
  output.name <- paste0(KO, "_", gsub(" ", "", gsub("[[:punct:]]", "_", gws)))


  gwscan <- data.frame()

  # combine all chr results into 1 file
  if (file.exists(paste0(results.folder, "raw_results/RDS_files/", gws, "_with_add_dom.rds"))) {
    gwscan <- readRDS(paste0(results.folder, "raw_results/RDS_files/", gws, "_with_add_dom.rds"))
    save(gwscan, file = paste0(results.folder, "raw_results/Rdata_files/", output.name, "_all_chr.RData"))
  } else if (file.exists(paste0(results.folder, "raw_results/RDS_files/", gws, "_chr_1with_add_dom.rds"))) {
    for (chr in 1:19) {
      tx <- readRDS(paste0(results.folder, gws, "_chr_", chr, "with_add_dom.rds"))
      gwscan <- rbind(gwscan, tx)
    }
    save(gwscan, file = paste0(results.folder, "raw_results/Rdata_files/", output.name, "_all_chr.RData"))
  } else if (file.exists(paste0(results.folder, "raw_results/RDS_files/", gsub("[[:punct:]]", "", gws), "_with_add_dom.rds"))) {
    gwscan <- readRDS(paste0(results.folder, gsub("[[:punct:]]", "", gws), "_with_add_dom.rds"))
    save(gwscan, file = paste0(results.folder, "raw_results/Rdata_files/", output.name, "_all_chr.RData"))
  } else if (file.exists(paste0(results.folder, "raw_results/RDS_files/", gsub("\\/", "", gws), "_with_add_dom.rds"))) {
    gwscan <- readRDS(paste0(results.folder, gsub("\\/", "", gws), "_with_add_dom.rds"))
    save(gwscan, file = paste0(results.folder, "raw_results/Rdata_files/", output.name, "_all_chr.RData"))
  } else {
    print(paste(gws, " File does not exist!"))
    trait_list <- trait_list[-i, ]
  }
}



for (i in 1:nrow(trait_list)) {
  gws <- trait_list$Identifier[i]
  output.name <- paste0(KO, "_", gsub(" ", "", gsub("[[:punct:]]", "_", gws)))



  chrom_range <- c(1:19)

  gwscan <- loadRData(paste0(results.folder, "raw_results/Rdata_files/", output.name, "_all_chr.RData"))

  gwscan <- gwscan[, 1:20] # remove index column

  # round to 3 sig figs
  round.colnames <- c("add.Beta", "add.StdErr", "add.T", "dom.Beta", "dom.StdErr", "dom.T", "P", "add.P", "dom.P")
  gwscan[, round.colnames] <- signif(gwscan[, round.colnames], 3)
  colnames(gwscan)[which(colnames(gwscan) == "tax")] <- "trait"
  # change position from Mb to bp
  gwscan$pos <- gwscan$pos * 1e6

  # make file with all sig SNPs
  bon.threshold <- 0.05 / 30401
  Meff.threshold <- 0.05 / 4603

  study.wide.bon <- bon.threshold / 140.9546


  study.wide.meff <- Meff.threshold / 140.9546



  sig.rows.P.Meff <- which(gwscan$P < Meff.threshold)
  sig.rows.add.P.Meff <- which(gwscan$add.P < Meff.threshold)
  sig.rows.dom.P.Meff <- which(gwscan$dom.P < Meff.threshold)
  sig.rows.all <- unique(c(sig.rows.P.Meff, sig.rows.add.P.Meff, sig.rows.dom.P.Meff))
  if (length(sig.rows.all) > 0) {
    gwscan_sig <- gwscan[sig.rows.all, ]
    # make a numeric chr column for convenience
    gwscan_sig$chr.num <- gsub("X", "20", gwscan_sig$chr)
    gwscan_sig$chr.num <- as.numeric(gwscan_sig$chr.num)
    gwscan_sig <- gwscan_sig[order(gwscan_sig$chr.num, gwscan_sig$pos), ]
    gwscan_sig$P.sig.Meff <- "N"
    gwscan_sig$add.P.sig.Meff <- "N"
    gwscan_sig$dom.P.sig.Meff <- "N"
    if (sum(gwscan_sig$P < Meff.threshold, na.rm = T) >= 1) gwscan_sig[is.na(gwscan_sig$P) == F & gwscan_sig$P < Meff.threshold, ]$P.sig.Meff <- "Y"
    if (sum(gwscan_sig$add.P < Meff.threshold, na.rm = T) >= 1) gwscan_sig[is.na(gwscan_sig$add.P) == F & gwscan_sig$add.P < Meff.threshold, ]$add.P.sig.Meff <- "Y"
    if (sum(gwscan_sig$dom.P < Meff.threshold, na.rm = T) >= 1) gwscan_sig[is.na(gwscan_sig$dom.P) == F & gwscan_sig$dom.P < Meff.threshold, ]$dom.P.sig.Meff <- "Y"

    gwscan_sig$P.sig.bon <- "N"
    gwscan_sig$add.P.sig.bon <- "N"
    gwscan_sig$dom.P.sig.bon <- "N"
    if (sum(gwscan_sig$P < bon.threshold, na.rm = T) >= 1) gwscan_sig[is.na(gwscan_sig$P) == F & gwscan_sig$P < bon.threshold, ]$P.sig.bon <- "Y"
    if (sum(gwscan_sig$add.P < bon.threshold, na.rm = T) >= 1) gwscan_sig[is.na(gwscan_sig$add.P) == F & gwscan_sig$add.P < bon.threshold, ]$add.P.sig.bon <- "Y"
    if (sum(gwscan_sig$dom.P < bon.threshold, na.rm = T) >= 1) gwscan_sig[is.na(gwscan_sig$dom.P) == F & gwscan_sig$dom.P < bon.threshold, ]$dom.P.sig.bon <- "Y"
    # study-wide significant
    gwscan_sig$P.sig.Meff_SW <- "N"
    gwscan_sig$add.P.sig.Meff_SW <- "N"
    gwscan_sig$dom.P.sig.Meff_SW <- "N"
    if (sum(gwscan_sig$P < study.wide.meff, na.rm = T) >= 1) gwscan_sig[is.na(gwscan_sig$P) == F & gwscan_sig$P < study.wide.meff, ]$P.sig.Meff_SW <- "Y"
    if (sum(gwscan_sig$add.P < study.wide.meff, na.rm = T) >= 1) gwscan_sig[is.na(gwscan_sig$add.P) == F & gwscan_sig$add.P < study.wide.meff, ]$add.P.sig.Meff_SW <- "Y"
    if (sum(gwscan_sig$dom.P < study.wide.meff, na.rm = T) >= 1) gwscan_sig[is.na(gwscan_sig$dom.P) == F & gwscan_sig$dom.P < study.wide.meff, ]$dom.P.sig.Meff_SW <- "Y"

    gwscan_sig$P.sig.bon_SW <- "N"
    gwscan_sig$add.P.sig.bon_SW <- "N"
    gwscan_sig$dom.P.sig.bon_SW <- "N"
    if (sum(gwscan_sig$P < study.wide.bon, na.rm = T) >= 1) gwscan_sig[is.na(gwscan_sig$P) == F & gwscan_sig$P < study.wide.bon, ]$P.sig.bon_SW <- "Y"
    if (sum(gwscan_sig$add.P < study.wide.bon, na.rm = T) >= 1) gwscan_sig[is.na(gwscan_sig$add.P) == F & gwscan_sig$add.P < study.wide.bon, ]$add.P.sig.bon_SW <- "Y"
    if (sum(gwscan_sig$dom.P < study.wide.bon, na.rm = T) >= 1) gwscan_sig[is.na(gwscan_sig$dom.P) == F & gwscan_sig$dom.P < study.wide.bon, ]$dom.P.sig.bon_SW <- "Y"
  }
  assign(paste(output.name, "_gwscan_sigSNPs", sep = ""), gwscan_sig)
  save(list = paste(output.name, "_gwscan_sigSNPs", sep = ""), file = paste(output.name, "_gwscan_sigSNPs.Rdata", sep = ""))
  rm(list = paste(output.name, "_gwscan_sigSNPs", sep = ""))
}


# Get significant intervals & annotate ####
for (i in 1:nrow(trait_list)) {
  gws <- trait_list$Identifier[i]
  output.name <- paste0(KO, "_", gsub(" ", "", gsub("[[:punct:]]", "_", gws)))

  load(paste0(output.name, "_gwscan_sigSNPs.Rdata"))
  gwscan_sig <- eval(parse(text = paste0(output.name, "_gwscan_sigSNPs")))
  rm(list = paste(output.name, "_gwscan_sigSNPs", sep = ""))
  sig.cats <- "n"
  for (sig.type in c("Meff", "bon")) {
    for (P.type in c("P", "add.P", "dom.P")) {
      if (length(which(gwscan_sig[, paste0(P.type, ".sig.", sig.type)] == "Y")) > 0) {
        sig.cats <- c(sig.cats, paste0(P.type, ".", sig.type))
        gwscan_sig_sub <- gwscan_sig[gwscan_sig[, paste0(P.type, ".sig.", sig.type)] == "Y", ]
        sig.chrs <- unique(gwscan_sig_sub$chr.num)
        for (chr in sig.chrs) {
          gwscan_sig_chr <- gwscan_sig_sub[gwscan_sig_sub$chr.num == chr, ]
          gwscan_sig_chr$bin <- 1
          # combine sig snps <10 Mb apart into bins
          if (nrow(gwscan_sig_chr) > 1) {
            for (j in 2:nrow(gwscan_sig_chr)) {
              gwscan_sig_chr$bin[j] <- gwscan_sig_chr$bin[j - 1]
              if (gwscan_sig_chr$pos[j] - gwscan_sig_chr$pos[j - 1] > 10 * 1e6) gwscan_sig_chr$bin[j] <- gwscan_sig_chr$bin[j] + 1
            }
          }
          if (chr == sig.chrs[1]) gwscan_sigsnp_bins <- gwscan_sig_chr
          if (chr != sig.chrs[1]) gwscan_sigsnp_bins <- rbind(gwscan_sigsnp_bins, gwscan_sig_chr)
        }
        gwscan_sigsnp_bins$chr_bin <- paste(gwscan_sigsnp_bins$chr.num, gwscan_sigsnp_bins$bin, sep = "_")
        sig.bins <- data.frame(chr_bin = unique(gwscan_sigsnp_bins$chr_bin), chr.num = 0, start.snp = "a", start.pos = 0, peak.snp = "a", peak.Pval = 0, peak.pos = 0, stop.snp = "a", stop.pos = 0, sig.snps = 0)
        for (k in 1:nrow(sig.bins)) {
          snps.sub <- gwscan_sigsnp_bins[gwscan_sigsnp_bins$chr_bin == sig.bins$chr_bin[k], ]
          sig.bins$chr.num[k] <- snps.sub$chr.num[1]
          sig.bins$start.pos[k] <- min(snps.sub$pos)
          sig.bins$stop.pos[k] <- max(snps.sub$pos)
          sig.bins$start.snp[k] <- snps.sub[which.min(snps.sub$pos), ]$marker
          sig.bins$stop.snp[k] <- snps.sub[which.max(snps.sub$pos), ]$marker
          sig.bins$sig.snps[k] <- nrow(snps.sub)
          sig.bins$peak.snp[k] <- snps.sub[which.min(snps.sub[, P.type]), ]$marker
          sig.bins$peak.pos[k] <- snps.sub[which.min(snps.sub[, P.type]), ]$pos
          sig.bins$peak.Pval[k] <- snps.sub[which.min(snps.sub[, P.type]), P.type]
        }
        # go back to pre-LD-filtered snps (hybrid_f2_5 from cleaning) to extend each bin to include SNPs in LD>0.9 with start/stop of interval #####Note Shauni's script used hybrid_f2
        sig.bins$start.LD.pos <- sig.bins$start.pos
        sig.bins$stop.LD.pos <- sig.bins$stop.pos
        for (l in 1:nrow(sig.bins)) {
          snp <- sig.bins$start.snp[l]
          system(paste0(plink, " --bfile ~/Documents/Research/Experiments/Final_QTL_mapping/Scripts/Data/genotypes/hybrid_f2_5 -r2 --ld-snp ", snp, " --ld-window-kb 10000 --ld-window-r2 0.9 --ld-window 99999 --allow-extra-chr"))
          ld_set <- read.table("plink.ld", header = T)
          rownames(ld_set) <- ld_set$SNP_B
          ld_set$marker <- ld_set$SNP_B
          ld_set <- match_df(preLDpruneSNPs, ld_set, on = "marker")
          sig.bins$start.LD.pos[l] <- min(c(ld_set$pos, sig.bins$start.pos[l]))
          if (sig.bins$sig.snps[l] == 1) sig.bins$stop.LD.pos[l] <- max(c(ld_set$pos, sig.bins$stop.pos[l]))
          if (sig.bins$sig.snps[l] > 1) {
            snp <- sig.bins$stop.snp[l]
            system(paste0(plink, " --bfile ~/Documents/Research/Experiments/Final_QTL_mapping/Scripts/Data/genotypes/hybrid_f2_5 -r2 --ld-snp ", snp, " --ld-window-kb 10000 --ld-window-r2 0.9 --ld-window 99999 --allow-extra-chr"))
            ld_set <- read.table("plink.ld", header = T)
            rownames(ld_set) <- ld_set$SNP_B
            ld_set$marker <- ld_set$SNP_B
            ld_set <- match_df(preLDpruneSNPs, ld_set, on = "marker")
            sig.bins$stop.LD.pos[l] <- max(c(ld_set$pos, sig.bins$stop.pos[l]))
          }
        }
        # check to make sure all intervals are still >10Mb apart after LD expansion
        sig.bin.dist <- rep(100 * 1e6, nrow(sig.bins))
        if (nrow(sig.bins) == 1) sig.bins.all <- sig.bins
        if (nrow(sig.bins) > 1) {
          for (l in 2:nrow(sig.bins)) {
            if (sig.bins$chr.num[l] == sig.bins$chr.num[l - 1]) sig.bin.dist[l] <- sig.bins$start.LD.pos[l] - sig.bins$stop.LD.pos[l - 1]
          }
          redo.bin <- which(sig.bin.dist > 0 & sig.bin.dist < 10 * 1e6)
          if (length(redo.bin) == 0) sig.bins.all <- sig.bins
          if (length(redo.bin) > 0) {
            chr.redo <- unique(sig.bins[redo.bin, ]$chr.num)
            for (chr in sig.chrs) {
              if (chr %in% chr.redo == F) sig.bins.chr <- sig.bins[sig.bins$chr.num == chr, ]
              if (chr %in% chr.redo == T) {
                sig.bins.chr.old <- sig.bins[sig.bins$chr.num == chr, ]
                sig.bins.chr.old$bin.new <- 1
                for (m in 2:nrow(sig.bins.chr.old)) {
                  sig.bins.chr.old$bin.new[m] <- sig.bins.chr.old$bin.new[m - 1]
                  if (sig.bins.chr.old$start.LD.pos[m] - sig.bins.chr.old$stop.LD.pos[m - 1] > 10 * 1e6) sig.bins.chr.old$bin.new[m] <- sig.bins.chr.old$bin.new[m] + 1
                  sig.bins.chr.old$chr_bin.new <- paste(sig.bins.chr.old$chr.num, sig.bins.chr.old$bin.new, sep = "_")
                }
                sig.bins.chr <- data.frame(chr_bin = unique(sig.bins.chr.old$chr_bin.new), chr.num = chr, start.snp = "a", start.pos = 0, peak.snp = "a", peak.Pval = 0, peak.pos = 0, stop.snp = "a", stop.pos = 0, sig.snps = 0, start.LD.pos = 0, stop.LD.pos = 0)
                for (k in 1:nrow(sig.bins.chr)) {
                  sig.bins.chr$start.pos[k] <- min(sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ]$start.pos)
                  sig.bins.chr$stop.pos[k] <- max(sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ]$stop.pos)
                  sig.bins.chr$start.snp[k] <- sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ]$start.snp[1]
                  sig.bins.chr$stop.snp[k] <- sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ]$start.snp[nrow(sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ])]
                  sig.bins.chr$start.LD.pos[k] <- min(sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ]$start.LD.pos)
                  sig.bins.chr$stop.LD.pos[k] <- max(sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ]$stop.LD.pos)
                  sig.bins.chr$sig.snps[k] <- nrow(gwscan_sig_sub[gwscan_sig_sub$chr.num == chr & gwscan_sig_sub$pos >= sig.bins.chr$start.LD.pos[k] & gwscan_sig_sub$pos <= sig.bins.chr$stop.LD.pos[k], ])
                  peak.row <- which.min(sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ]$peak.Pval)
                  sig.bins.chr$peak.snp[k] <- sig.bins.chr.old$peak.snp[peak.row]
                  sig.bins.chr$peak.Pval[k] <- sig.bins.chr.old$peak.Pval[peak.row]
                  sig.bins.chr$peak.pos[k] <- sig.bins.chr.old$peak.pos[peak.row]
                }
              }
              if (chr == sig.chrs[1]) sig.bins.all <- sig.bins.chr
              if (chr != sig.chrs[1]) sig.bins.all <- rbind(sig.bins.all, sig.bins.chr)
            }
          }
        }
        sig.bins.all$P.type <- P.type
        sig.bins.all$sig.type <- sig.type
        sig.bins.all$trait <- gws
        assign(paste(output.name, "_gwscan_sigIntervals.", P.type, ".", sig.type, sep = ""), sig.bins.all)
        save(list = paste(output.name, "_gwscan_sigIntervals.", P.type, ".", sig.type, sep = ""), file = paste(output.name, "_gwscan_sigIntervals.", P.type, ".", sig.type, ".Rdata", sep = ""))
      }
    }
  }
  # Annotation of significant intervals
  if (length(sig.cats) > 1) {
    sig.cats <- sig.cats[-1]
    sig.intervals.comb <- get(paste0(output.name, "_gwscan_sigIntervals.", sig.cats[1]))
    if (length(sig.cats) > 1) {
      for (k in 2:length(sig.cats)) {
        sig.intervals.comb <- rbind(sig.intervals.comb, get(paste0(output.name, "_gwscan_sigIntervals.", sig.cats[k])))
      }
    }
    sig.intervals.comb$chr <- sig.intervals.comb$chr.num
    if (20 %in% sig.intervals.comb$chr.num == T) sig.intervals.comb$chr <- gsub("20", "X", as.character(sig.intervals.comb$chr.num))
    rownames(sig.intervals.comb) <- 1:nrow(sig.intervals.comb)
    for (j in 1:nrow(sig.intervals.comb)) {
      P.type <- sig.intervals.comb$P.type[j]
      sig.type <- sig.intervals.comb$sig.type[j]
      max_snp <- sig.intervals.comb$peak.snp[j]
      pos_max_snp <- sig.intervals.comb$peak.pos[j]
      gwscan_sig_sub <- gwscan_sig[gwscan_sig[, paste0(P.type, ".sig.", sig.type)] == "Y", ]
      # get genes in region using bioMart
      out.bm.genes.region <- getBM(
        attributes = c("start_position", "end_position", "ensembl_gene_id", "external_gene_name", "gene_biotype", "description"),
        filters = c("chromosome_name", "start", "end"),
        values = list(sig.intervals.comb$chr[j], sig.intervals.comb$start.LD.pos[j], sig.intervals.comb$stop.LD.pos[j]),
        mart = gene.ensembl
      )
      interval.annot <- data.frame(
        chr_bin = sig.intervals.comb$chr_bin[j], P.type = sig.intervals.comb$P.type[j], peak_snp = max_snp,
        rsId = preLDpruneSNPs[max_snp, "rsID"], add.Beta_peak_snp = gwscan_sig_sub[max_snp, "add.Beta"], add.StdErr_peak_snp = gwscan_sig_sub[max_snp, "add.StdErr"], add.Zscore_peak_snp = gwscan_sig_sub[max_snp, "add.T"], add.P_peak_snp = gwscan_sig_sub[max_snp, "add.P"], dom.Beta_peak_snp = gwscan_sig_sub[max_snp, "dom.Beta"], dom.StdErr_peak_snp = gwscan_sig_sub[max_snp, "dom.StdErr"], dom.Zscore_peak_snp = gwscan_sig_sub[max_snp, "dom.T"], dom.P_peak_snp = gwscan_sig_sub[max_snp, "dom.P"],
        markers = paste(gwscan_sig_sub[gwscan_sig_sub$pos >= sig.intervals.comb$start.LD.pos[j] & gwscan_sig_sub$pos <= sig.intervals.comb$stop.LD.pos[j], ]$marker, collapse = ", "),
        closest_gene = "", closest_gene_type = "",
        closest_gene_description = "",
        total_genes = nrow(out.bm.genes.region),
        list_total_genes = "",
        total_protein_coding_genes = 0,
        list_protein_coding_genes = ""
      )
      if ((dim(out.bm.genes.region)[1] > 0)) {
        absol_diff <- matrix(c(
          abs(out.bm.genes.region$start_position - pos_max_snp),
          abs(out.bm.genes.region$end_position - pos_max_snp)
        ), ncol = 2, byrow = F)
        closest_gene <- out.bm.genes.region[which(absol_diff == min(absol_diff), arr.ind = T)[1], ]
        interval.annot$closest_gene[1] <- paste(closest_gene$external_gene_name, collapse = ", ")
        interval.annot$closest_gene_type[1] <- paste(closest_gene$gene_biotype, collapse = ", ")
        interval.annot$closest_gene_description[1] <- paste(closest_gene$description, collapse = ", ")
        interval.annot$list_total_genes[1] <- paste(out.bm.genes.region$external_gene_name, collapse = ", ")
        interval.annot$total_protein_coding_genes[1] <- nrow(out.bm.genes.region[which(out.bm.genes.region$gene_biotype == "protein_coding"), ])
        interval.annot$list_protein_coding_genes[1] <- paste(out.bm.genes.region[which(out.bm.genes.region$gene_biotype == "protein_coding"), "external_gene_name"], collapse = ", ")
      }
      # For intervals not containing any genes - get closest gene within 0.5 Mb on either side of peak SNP
      if (dim(out.bm.genes.region)[1] == 0) {
        nearby.genes <- getBM(
          attributes = c("start_position", "end_position", "ensembl_gene_id", "external_gene_name", "gene_biotype", "description"),
          filters = c("chromosome_name", "start", "end"),
          values = list(sig.intervals.comb$chr[j], sig.intervals.comb$start.LD.pos[j] - 5e5, sig.intervals.comb$stop.LD.pos[j] + 5e5),
          mart = gene.ensembl
        )
        if (nrow(nearby.genes) > 0) {
          absol_diff <- matrix(c(
            abs(nearby.genes$start_position - pos_max_snp),
            abs(nearby.genes$end_position - pos_max_snp)
          ), ncol = 2, byrow = F)
          closest_gene <- nearby.genes[which(absol_diff == min(absol_diff), arr.ind = T)[1], ]
          interval.annot$closest_gene[1] <- paste(closest_gene$external_gene_name, collapse = ", ")
          interval.annot$closest_gene_type[1] <- paste(closest_gene$gene_biotype, collapse = ", ")
          interval.annot$closest_gene_description[1] <- paste(closest_gene$description, collapse = ", ")
        }
      }
      if (j == 1) sig.intervals.annot <- interval.annot
      if (j > 1) sig.intervals.annot <- rbind(sig.intervals.annot, interval.annot)
    }
    sig.intervals.out <- cbind(sig.intervals.comb, sig.intervals.annot)
    sig.intervals.out$length_Mb <- (sig.intervals.out$stop.LD.pos - sig.intervals.out$start.LD.pos) / 1e6
    sig.intervals.out
    # reorder and delete some columns
    keep.cols <- c("trait", "P.type", "sig.type", "chr", "chr.num", "chr_bin", "start.LD.pos", "stop.LD.pos", "length_Mb", "sig.snps", "peak.snp", "rsId", "peak.pos", "peak.Pval", "add.Beta_peak_snp", "add.StdErr_peak_snp", "add.Zscore_peak_snp", "add.P_peak_snp", "dom.Beta_peak_snp", "dom.StdErr_peak_snp", "dom.Zscore_peak_snp", "dom.P_peak_snp", "markers", "closest_gene", "closest_gene_type", "closest_gene_description", "total_genes", "list_total_genes", "total_protein_coding_genes", "list_protein_coding_genes")
    sig.intervals.out <- sig.intervals.out[, keep.cols]
    # rename some columns
    colnames(sig.intervals.out)[which(colnames(sig.intervals.out) == "start.LD.pos")] <- "start"
    colnames(sig.intervals.out)[which(colnames(sig.intervals.out) == "stop.LD.pos")] <- "stop"
    write.csv(sig.intervals.out, file = paste0(output.name, "_gwscan_sigIntervals.annot.csv"))
  }
  rm(list = ls()[grep(output.name, ls())])
}


# STUDY-WIDE significant intervals ####

# Get significant intervals & annotate
# Plink slows down after many traits for some reason - running ~5 at a time instead of whole trait list is faster
plink <- "/Users/doms/Documents/Research/Software/plink/plink"
for (i in 1:nrow(trait_list)) {
  gws <- trait_list$Identifier[i]
  output.name <- paste0(KO, "_", gsub(" ", "", gsub("[[:punct:]]", "_", gws)))
  # output.name <- paste0(KO, "_", gsub(" ", "", gws)) # for EC level 1

  load(paste0(output.name, "_gwscan_sigSNPs.Rdata"))
  gwscan_sig <- eval(parse(text = paste0(output.name, "_gwscan_sigSNPs")))

  rm(list = paste(output.name, "_gwscan_sigSNPs", sep = ""))
  sig.cats <- "n"
  for (sig.type in c("Meff", "bon")) {
    for (P.type in c("P", "add.P", "dom.P")) {
      if (length(which(gwscan_sig[, paste0(P.type, ".sig.", sig.type, "_SW")] == "Y")) > 0) {
        sig.cats <- c(sig.cats, paste0(P.type, ".", sig.type))
        gwscan_sig_sub <- gwscan_sig[gwscan_sig[, paste0(P.type, ".sig.", sig.type, "_SW")] == "Y", ]
        sig.chrs <- unique(gwscan_sig_sub$chr.num)
        for (chr in sig.chrs) {
          gwscan_sig_chr <- gwscan_sig_sub[gwscan_sig_sub$chr.num == chr, ]
          gwscan_sig_chr$bin <- 1
          # combine sig snps <10 Mb apart into bins
          if (nrow(gwscan_sig_chr) > 1) {
            for (j in 2:nrow(gwscan_sig_chr)) {
              gwscan_sig_chr$bin[j] <- gwscan_sig_chr$bin[j - 1]
              if (gwscan_sig_chr$pos[j] - gwscan_sig_chr$pos[j - 1] > 10 * 1e6) gwscan_sig_chr$bin[j] <- gwscan_sig_chr$bin[j] + 1
            }
          }
          if (chr == sig.chrs[1]) gwscan_sigsnp_bins <- gwscan_sig_chr
          if (chr != sig.chrs[1]) gwscan_sigsnp_bins <- rbind(gwscan_sigsnp_bins, gwscan_sig_chr)
        }
        gwscan_sigsnp_bins$chr_bin <- paste(gwscan_sigsnp_bins$chr.num, gwscan_sigsnp_bins$bin, sep = "_")
        sig.bins <- data.frame(chr_bin = unique(gwscan_sigsnp_bins$chr_bin), chr.num = 0, start.snp = "a", start.pos = 0, peak.snp = "a", peak.Pval = 0, peak.pos = 0, stop.snp = "a", stop.pos = 0, sig.snps = 0)
        for (k in 1:nrow(sig.bins)) {
          snps.sub <- gwscan_sigsnp_bins[gwscan_sigsnp_bins$chr_bin == sig.bins$chr_bin[k], ]
          sig.bins$chr.num[k] <- snps.sub$chr.num[1]
          sig.bins$start.pos[k] <- min(snps.sub$pos)
          sig.bins$stop.pos[k] <- max(snps.sub$pos)
          sig.bins$start.snp[k] <- snps.sub[which.min(snps.sub$pos), ]$marker
          sig.bins$stop.snp[k] <- snps.sub[which.max(snps.sub$pos), ]$marker
          sig.bins$sig.snps[k] <- nrow(snps.sub)
          sig.bins$peak.snp[k] <- snps.sub[which.min(snps.sub[, P.type]), ]$marker
          sig.bins$peak.pos[k] <- snps.sub[which.min(snps.sub[, P.type]), ]$pos
          sig.bins$peak.Pval[k] <- snps.sub[which.min(snps.sub[, P.type]), P.type]
        }
        # go back to pre-LD-filtered snps (hybrid_f2_5 from cleaning) to extend each bin to include SNPs in LD>0.9 with start/stop of interval #####Note Shauni's script used hybrid_f2
        sig.bins$start.LD.pos <- sig.bins$start.pos
        sig.bins$stop.LD.pos <- sig.bins$stop.pos
        for (l in 1:nrow(sig.bins)) {
          snp <- sig.bins$start.snp[l]
          system(paste0(plink, " --bfile ~/Documents/Research/Experiments/Final_QTL_mapping/Scripts/Data/genotypes/hybrid_f2_5 -r2 --ld-snp ", snp, " --ld-window-kb 10000 --ld-window-r2 0.9 --ld-window 99999 --allow-extra-chr"))
          ld_set <- read.table("plink.ld", header = T)
          rownames(ld_set) <- ld_set$SNP_B
          ld_set$marker <- ld_set$SNP_B
          ld_set <- match_df(preLDpruneSNPs, ld_set, on = "marker")
          sig.bins$start.LD.pos[l] <- min(c(ld_set$pos, sig.bins$start.pos[l]))
          if (sig.bins$sig.snps[l] == 1) sig.bins$stop.LD.pos[l] <- max(c(ld_set$pos, sig.bins$stop.pos[l]))
          if (sig.bins$sig.snps[l] > 1) {
            snp <- sig.bins$stop.snp[l]
            system(paste0(plink, " --bfile ~/Documents/Research/Experiments/Final_QTL_mapping/Scripts/Data/genotypes/hybrid_f2_5 -r2 --ld-snp ", snp, " --ld-window-kb 10000 --ld-window-r2 0.9 --ld-window 99999 --allow-extra-chr"))
            ld_set <- read.table("plink.ld", header = T)
            rownames(ld_set) <- ld_set$SNP_B
            ld_set$marker <- ld_set$SNP_B
            ld_set <- match_df(preLDpruneSNPs, ld_set, on = "marker")
            sig.bins$stop.LD.pos[l] <- max(c(ld_set$pos, sig.bins$stop.pos[l]))
          }
        }
        # check to make sure all intervals are still >10Mb apart after LD expansion
        sig.bin.dist <- rep(100 * 1e6, nrow(sig.bins))
        if (nrow(sig.bins) == 1) sig.bins.all <- sig.bins
        if (nrow(sig.bins) > 1) {
          for (l in 2:nrow(sig.bins)) {
            if (sig.bins$chr.num[l] == sig.bins$chr.num[l - 1]) sig.bin.dist[l] <- sig.bins$start.LD.pos[l] - sig.bins$stop.LD.pos[l - 1]
          }
          redo.bin <- which(sig.bin.dist > 0 & sig.bin.dist < 10 * 1e6)
          if (length(redo.bin) == 0) sig.bins.all <- sig.bins
          if (length(redo.bin) > 0) {
            chr.redo <- unique(sig.bins[redo.bin, ]$chr.num)
            for (chr in sig.chrs) {
              if (chr %in% chr.redo == F) sig.bins.chr <- sig.bins[sig.bins$chr.num == chr, ]
              if (chr %in% chr.redo == T) {
                sig.bins.chr.old <- sig.bins[sig.bins$chr.num == chr, ]
                sig.bins.chr.old$bin.new <- 1
                for (m in 2:nrow(sig.bins.chr.old)) {
                  sig.bins.chr.old$bin.new[m] <- sig.bins.chr.old$bin.new[m - 1]
                  if (sig.bins.chr.old$start.LD.pos[m] - sig.bins.chr.old$stop.LD.pos[m - 1] > 10 * 1e6) sig.bins.chr.old$bin.new[m] <- sig.bins.chr.old$bin.new[m] + 1
                  sig.bins.chr.old$chr_bin.new <- paste(sig.bins.chr.old$chr.num, sig.bins.chr.old$bin.new, sep = "_")
                }
                sig.bins.chr <- data.frame(chr_bin = unique(sig.bins.chr.old$chr_bin.new), chr.num = chr, start.snp = "a", start.pos = 0, peak.snp = "a", peak.Pval = 0, peak.pos = 0, stop.snp = "a", stop.pos = 0, sig.snps = 0, start.LD.pos = 0, stop.LD.pos = 0)
                for (k in 1:nrow(sig.bins.chr)) {
                  sig.bins.chr$start.pos[k] <- min(sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ]$start.pos)
                  sig.bins.chr$stop.pos[k] <- max(sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ]$stop.pos)
                  sig.bins.chr$start.snp[k] <- sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ]$start.snp[1]
                  sig.bins.chr$stop.snp[k] <- sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ]$start.snp[nrow(sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ])]
                  sig.bins.chr$start.LD.pos[k] <- min(sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ]$start.LD.pos)
                  sig.bins.chr$stop.LD.pos[k] <- max(sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ]$stop.LD.pos)
                  sig.bins.chr$sig.snps[k] <- nrow(gwscan_sig_sub[gwscan_sig_sub$chr.num == chr & gwscan_sig_sub$pos >= sig.bins.chr$start.LD.pos[k] & gwscan_sig_sub$pos <= sig.bins.chr$stop.LD.pos[k], ])
                  peak.row <- which.min(sig.bins.chr.old[sig.bins.chr.old$chr_bin.new == sig.bins.chr$chr_bin[k], ]$peak.Pval)
                  sig.bins.chr$peak.snp[k] <- sig.bins.chr.old$peak.snp[peak.row]
                  sig.bins.chr$peak.Pval[k] <- sig.bins.chr.old$peak.Pval[peak.row]
                  sig.bins.chr$peak.pos[k] <- sig.bins.chr.old$peak.pos[peak.row]
                }
              }
              if (chr == sig.chrs[1]) sig.bins.all <- sig.bins.chr
              if (chr != sig.chrs[1]) sig.bins.all <- rbind(sig.bins.all, sig.bins.chr)
            }
          }
        }
        sig.bins.all$P.type <- P.type
        sig.bins.all$sig.type <- sig.type
        sig.bins.all$trait <- paste0(trait_list$Module[i], "_", gsub(" ", "_", gsub(pattern = "[^[:alnum:] ]", "_", trait_list$Description[i])))
        assign(paste(output.name, "_gwscan_sigIntervals.", P.type, ".", sig.type, "_SW", sep = ""), sig.bins.all)
        save(list = paste(output.name, "_gwscan_sigIntervals.", P.type, ".", sig.type, "_SW", sep = ""), file = paste(output.name, "_gwscan_sigIntervals.", P.type, ".", sig.type, "_SW", ".Rdata", sep = ""))
      }
    }
  }
  # Annotation of significant intervals
  if (length(sig.cats) > 1) {
    sig.cats <- sig.cats[-1]
    sig.intervals.comb <- get(paste0(output.name, "_gwscan_sigIntervals.", sig.cats[1], "_SW"))
    if (length(sig.cats) > 1) {
      for (k in 2:length(sig.cats)) {
        sig.intervals.comb <- rbind(sig.intervals.comb, get(paste0(output.name, "_gwscan_sigIntervals.", sig.cats[k], "_SW")))
      }
    }
    sig.intervals.comb$chr <- sig.intervals.comb$chr.num
    if (20 %in% sig.intervals.comb$chr.num == T) sig.intervals.comb$chr <- gsub("20", "X", as.character(sig.intervals.comb$chr.num))
    rownames(sig.intervals.comb) <- 1:nrow(sig.intervals.comb)
    for (j in 1:nrow(sig.intervals.comb)) {
      P.type <- sig.intervals.comb$P.type[j]
      sig.type <- sig.intervals.comb$sig.type[j]
      max_snp <- sig.intervals.comb$peak.snp[j]
      pos_max_snp <- sig.intervals.comb$peak.pos[j]
      gwscan_sig_sub <- gwscan_sig[gwscan_sig[, paste0(P.type, ".sig.", sig.type)] == "Y", ]
      # get genes in region using bioMart
      out.bm.genes.region <- getBM(
        attributes = c("start_position", "end_position", "ensembl_gene_id", "external_gene_name", "gene_biotype", "description"),
        filters = c("chromosome_name", "start", "end"),
        values = list(sig.intervals.comb$chr[j], sig.intervals.comb$start.LD.pos[j], sig.intervals.comb$stop.LD.pos[j]),
        mart = gene.ensembl
      )
      interval.annot <- data.frame(
        chr_bin = sig.intervals.comb$chr_bin[j], P.type = sig.intervals.comb$P.type[j], peak_snp = max_snp,
        rsId = preLDpruneSNPs[max_snp, "rsID"], add.Beta_peak_snp = gwscan_sig_sub[max_snp, "add.Beta"], add.StdErr_peak_snp = gwscan_sig_sub[max_snp, "add.StdErr"], add.Zscore_peak_snp = gwscan_sig_sub[max_snp, "add.T"], add.P_peak_snp = gwscan_sig_sub[max_snp, "add.P"], dom.Beta_peak_snp = gwscan_sig_sub[max_snp, "dom.Beta"], dom.StdErr_peak_snp = gwscan_sig_sub[max_snp, "dom.StdErr"], dom.Zscore_peak_snp = gwscan_sig_sub[max_snp, "dom.T"], dom.P_peak_snp = gwscan_sig_sub[max_snp, "dom.P"],
        markers = paste(gwscan_sig_sub[gwscan_sig_sub$pos >= sig.intervals.comb$start.LD.pos[j] & gwscan_sig_sub$pos <= sig.intervals.comb$stop.LD.pos[j], ]$marker, collapse = ", "),
        closest_gene = "", closest_gene_type = "",
        closest_gene_description = "",
        total_genes = nrow(out.bm.genes.region),
        list_total_genes = "",
        total_protein_coding_genes = 0,
        list_protein_coding_genes = ""
      )
      if ((dim(out.bm.genes.region)[1] > 0)) {
        absol_diff <- matrix(c(
          abs(out.bm.genes.region$start_position - pos_max_snp),
          abs(out.bm.genes.region$end_position - pos_max_snp)
        ), ncol = 2, byrow = F)
        closest_gene <- out.bm.genes.region[which(absol_diff == min(absol_diff), arr.ind = T)[1], ]
        interval.annot$closest_gene[1] <- paste(closest_gene$external_gene_name, collapse = ", ")
        interval.annot$closest_gene_type[1] <- paste(closest_gene$gene_biotype, collapse = ", ")
        interval.annot$closest_gene_description[1] <- paste(closest_gene$description, collapse = ", ")
        interval.annot$list_total_genes[1] <- paste(out.bm.genes.region$external_gene_name, collapse = ", ")
        interval.annot$total_protein_coding_genes[1] <- nrow(out.bm.genes.region[which(out.bm.genes.region$gene_biotype == "protein_coding"), ])
        interval.annot$list_protein_coding_genes[1] <- paste(out.bm.genes.region[which(out.bm.genes.region$gene_biotype == "protein_coding"), "external_gene_name"], collapse = ", ")
      }
      # For intervals not containing any genes - get closest gene within 0.5 Mb on either side of peak SNP
      if (dim(out.bm.genes.region)[1] == 0) {
        nearby.genes <- getBM(
          attributes = c("start_position", "end_position", "ensembl_gene_id", "external_gene_name", "gene_biotype", "description"),
          filters = c("chromosome_name", "start", "end"),
          values = list(sig.intervals.comb$chr[j], sig.intervals.comb$start.LD.pos[j] - 5e5, sig.intervals.comb$stop.LD.pos[j] + 5e5),
          mart = gene.ensembl
        )
        if (nrow(nearby.genes) > 0) {
          absol_diff <- matrix(c(
            abs(nearby.genes$start_position - pos_max_snp),
            abs(nearby.genes$end_position - pos_max_snp)
          ), ncol = 2, byrow = F)
          closest_gene <- nearby.genes[which(absol_diff == min(absol_diff), arr.ind = T)[1], ]
          interval.annot$closest_gene[1] <- paste(closest_gene$external_gene_name, collapse = ", ")
          interval.annot$closest_gene_type[1] <- paste(closest_gene$gene_biotype, collapse = ", ")
          interval.annot$closest_gene_description[1] <- paste(closest_gene$description, collapse = ", ")
        }
      }
      if (j == 1) sig.intervals.annot <- interval.annot
      if (j > 1) sig.intervals.annot <- rbind(sig.intervals.annot, interval.annot)
    }
    sig.intervals.out <- cbind(sig.intervals.comb, sig.intervals.annot)
    sig.intervals.out$length_Mb <- (sig.intervals.out$stop.LD.pos - sig.intervals.out$start.LD.pos) / 1e6
    sig.intervals.out
    # reorder and delete some columns
    keep.cols <- c("trait", "P.type", "sig.type", "chr", "chr.num", "chr_bin", "start.LD.pos", "stop.LD.pos", "length_Mb", "sig.snps", "peak.snp", "rsId", "peak.pos", "peak.Pval", "add.Beta_peak_snp", "add.StdErr_peak_snp", "add.Zscore_peak_snp", "add.P_peak_snp", "dom.Beta_peak_snp", "dom.StdErr_peak_snp", "dom.Zscore_peak_snp", "dom.P_peak_snp", "markers", "closest_gene", "closest_gene_type", "closest_gene_description", "total_genes", "list_total_genes", "total_protein_coding_genes", "list_protein_coding_genes")
    sig.intervals.out <- sig.intervals.out[, keep.cols]
    # rename some columns
    colnames(sig.intervals.out)[which(colnames(sig.intervals.out) == "start.LD.pos")] <- "start"
    colnames(sig.intervals.out)[which(colnames(sig.intervals.out) == "stop.LD.pos")] <- "stop"
    write.csv(sig.intervals.out, file = paste0(output.name, "_gwscan_sigIntervals_SW.annot.csv"))
  }
  rm(list = ls()[grep(output.name, ls())])
}
