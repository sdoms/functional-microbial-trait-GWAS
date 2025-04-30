setwd("~/Documents/Research/Experiments/shotgun_hybrids/mapping_results/ancestry/")
library(tidyverse)
library(patchwork)
library(ggsci)
library(paletteer)
theme_set(theme_test())

all_sig_snps <- read_csv("../all_sig_snps20250317.csv")


all_sig_loci <- read_csv("../all_sig_loci20250317.csv")


# ancestry

# load info about ancestry informative SNPs
load("snps.anc.inf.info.Rdata")
# load genotypes_ref which says if AA or BB is A1 or A2 -> AA is always minor allele
load("ref_alt_min_major.Rdata")
# add maj/min allele and ancestry informative to snp info
all_sig_snps <- merge(all_sig_snps, genotypes_ref[, 5:7], by = "marker")

all_sig_snps$anc.inf <- snps.anc.inf.info[all_sig_snps$marker, "anc.inf"]
all_sig_snps[is.na(all_sig_snps$anc.inf), ]$anc.inf <- "N"
all_sig_snps$dom.allele <- snps.anc.inf.info[all_sig_snps$marker, "dom.allele"]
all_sig_snps$mus.allele <- snps.anc.inf.info[all_sig_snps$marker, "mus.allele"]

##### add dominance ######
all_sig_snps$degree_of_dominance <- all_sig_snps$dom.T / abs(all_sig_snps$add.T)
all_sig_snps$deg_beta <- all_sig_snps$dom.Beta / abs(all_sig_snps$add.Beta)
all_sig_snps$dominance_cat <- cut(all_sig_snps$degree_of_dominance,
  breaks = c(-Inf, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, Inf),
  labels = c("underdominant", "recessive", "partially recessive", "additive", "partially dominant", "dominant", "overdominant")
)

# get DM recoded G2 genos for ancestry informative sig SNPs
sig.markers <- unique(all_sig_snps$marker)

# # read in phenotypes
phenotypes <- read.csv("../all_functional_traits.csv")

G2s <- phenotypes$X

# get sig snp genos native coded
load("clean_geno_20240327.Rdata")
colnames(clean_geno_20240327) <- gsub("/", ".", colnames(clean_geno_20240327))
sig.snps.geno.natcode <- clean_geno_20240327[sig.markers, G2s]
save(sig.snps.geno.natcode, file = "sig.snps.geno.natcode.Rdata")

# DM coded
informative_snps <- snps.anc.inf.info |>
  filter(anc.inf != "N") |>
  rownames_to_column("marker")
inform_geno <- clean_geno_20240327[rownames(clean_geno_20240327) %in% informative_snps$marker, ]

# Merge allele_table and coding_table by SNP
mapping <- right_join(genotypes_ref, informative_snps, by = "marker") |> filter(marker %in% rownames(inform_geno))

# Initialize a translated matrix
geno.anc.inf.DMcode <- inform_geno

# Iterate over SNPs
for (snp in mapping$marker) {
  # Extract mapping info
  snp_info <- mapping[mapping$marker == snp, ]
  allele_a <- snp_info$AA_min
  allele_b <- snp_info$BB_maj
  dom_allele <- snp_info$dom.allele
  mus_allele <- snp_info$mus.allele

  # Translate each entry in the genotype_matrix for the SNP
  geno.anc.inf.DMcode[snp, ] <- sapply(inform_geno[snp, ], function(geno) {
    if (geno == "H") {
      return("H") # Keep heterozygous as is
    } else if (geno == allele_a) {
      if (dom_allele == "A") {
        return("D")
      } else if (mus_allele == "A") {
        return("M")
      } else {
        return(NA)
      }
    } else if (geno == allele_b) {
      if (dom_allele == "B") {
        return("D")
      } else if (mus_allele == "B") {
        return("M")
      } else {
        return(NA)
      }
    } else {
      return(NA) # Invalid genotype
    }
  })
}

save(geno.anc.inf.DMcode, file = "geno.anc.inf.DMcode.RData")

load("geno.anc.inf.DMcode.RData")

sig.snps.geno.anc.inf.DMcode <- geno.anc.inf.DMcode[sig.markers[sig.markers %in% rownames(geno.anc.inf.DMcode)], G2s]
save(sig.snps.geno.anc.inf.DMcode, file = "sig.snps.geno.anc.inf.DMcode.Rdata")


# get pheno means/geno for each sig snp
sig.traits <- unique(all_sig_snps$trait)
phenos.map <- phenotypes


# get pheno mean by geno nativecode & D/H/M
all_sig_snps$trait.mean.AA <- NA
all_sig_snps$trait.mean.AB <- NA
all_sig_snps$trait.mean.BB <- NA
all_sig_snps$N.AA <- NA
all_sig_snps$N.AB <- NA
all_sig_snps$N.BB <- NA
all_sig_snps$trait.mean.DD <- NA
all_sig_snps$trait.mean.DM <- NA
all_sig_snps$trait.mean.MM <- NA
all_sig_snps$N.DD <- NA
all_sig_snps$N.DM <- NA
all_sig_snps$N.MM <- NA


# transpose genos for easier subsetting
sig.snps.geno.transposed <- t(sig.snps.geno.natcode)
sig.snps.anc.inf.transposed <- t(sig.snps.geno.anc.inf.DMcode)
# make rownames mouse names for easier subsetting
rownames(phenos.map) <- phenos.map$X

for (row in 1:nrow(all_sig_snps)) {
  marker <- all_sig_snps$marker[row]
  trait <- all_sig_snps$trait[row]
  pheno.mice <- rownames(phenos.map)[!is.na(phenos.map[, trait])]
  AAmin.Allele <- all_sig_snps$AA_min[row]
  BBmaj.Allele <- all_sig_snps$BB_maj[row]
  AA.mice <- rownames(sig.snps.geno.transposed)[sig.snps.geno.transposed[, marker] == AAmin.Allele]
  AA.mice <- AA.mice[AA.mice %in% pheno.mice]
  BB.mice <- rownames(sig.snps.geno.transposed)[sig.snps.geno.transposed[, marker] == BBmaj.Allele]
  BB.mice <- BB.mice[BB.mice %in% pheno.mice]
  H.mice <- rownames(sig.snps.geno.transposed)[sig.snps.geno.transposed[, marker] == "H"]
  H.mice <- H.mice[H.mice %in% pheno.mice]
  all_sig_snps$trait.mean.AA[row] <- mean(phenos.map[AA.mice, trait])
  all_sig_snps$trait.mean.AB[row] <- mean(phenos.map[H.mice, trait])
  all_sig_snps$trait.mean.BB[row] <- mean(phenos.map[BB.mice, trait])
  all_sig_snps$N.AA[row] <- length(AA.mice)
  all_sig_snps$N.AB[row] <- length(H.mice)
  all_sig_snps$N.BB[row] <- length(BB.mice)
  if (all_sig_snps$anc.inf[row] != "N") {
    D.mice <- rownames(sig.snps.anc.inf.transposed)[sig.snps.anc.inf.transposed[, marker] == "D"]
    D.mice <- D.mice[D.mice %in% pheno.mice]
    M.mice <- rownames(sig.snps.anc.inf.transposed)[sig.snps.anc.inf.transposed[, marker] == "M"]
    M.mice <- M.mice[M.mice %in% pheno.mice]
    H.mice <- rownames(sig.snps.anc.inf.transposed)[sig.snps.anc.inf.transposed[, marker] == "H"]
    H.mice <- H.mice[H.mice %in% pheno.mice]
    all_sig_snps$trait.mean.DD[row] <- mean(phenos.map[D.mice, trait])
    all_sig_snps$trait.mean.DM[row] <- mean(phenos.map[H.mice, trait])
    all_sig_snps$trait.mean.MM[row] <- mean(phenos.map[M.mice, trait])
    all_sig_snps$N.DD[row] <- length(D.mice)
    all_sig_snps$N.DM[row] <- length(H.mice)
    all_sig_snps$N.MM[row] <- length(M.mice)
  }
}




genos.minmaj <- c("A", "H", "B")
genos.anc <- c("D", "H", "M")



all_sig_snps$dominance.d <- all_sig_snps$trait.mean.AB - ((all_sig_snps$trait.mean.AA + all_sig_snps$trait.mean.BB) / 2)
all_sig_snps$additive.a.wrtAA <- (all_sig_snps$trait.mean.AA - all_sig_snps$trait.mean.BB) / 2
all_sig_snps$additive.a.wrtDD <- (all_sig_snps$trait.mean.DD - all_sig_snps$trait.mean.MM) / 2
all_sig_snps$d.a.ratio <- all_sig_snps$dominance.d / abs(all_sig_snps$additive.a.wrtAA)
all_sig_snps$low.allele <- "A"
all_sig_snps[!is.na(all_sig_snps$additive.a.wrtAA) & all_sig_snps$additive.a.wrtAA > 0, ]$low.allele <- "B"
all_sig_snps[is.na(all_sig_snps$additive.a.wrtAA), ]$low.allele <- NA
all_sig_snps$low.geno <- all_sig_snps$low.allele
all_sig_snps[!is.na(all_sig_snps$d.a.ratio) & all_sig_snps$d.a.ratio < (-1), ]$low.geno <- "H"
all_sig_snps$low.allele.anc <- "D"
all_sig_snps[!is.na(all_sig_snps$additive.a.wrtDD) & all_sig_snps$additive.a.wrtDD > 0, ]$low.allele.anc <- "M"
all_sig_snps[is.na(all_sig_snps$additive.a.wrtDD), ]$low.allele.anc <- NA
all_sig_snps$low.geno.anc <- all_sig_snps$low.allele.anc
all_sig_snps[!is.na(all_sig_snps$dom.allele) & !is.na(all_sig_snps$d.a.ratio) & all_sig_snps$d.a.ratio < (-1), ]$low.geno.anc <- "H"
all_sig_snps$high.allele <- "B" # TEMP to make sterile geno designation easier
all_sig_snps[!is.na(all_sig_snps$additive.a.wrtAA) & all_sig_snps$additive.a.wrtAA > 0, ]$high.allele <- "A"
all_sig_snps[is.na(all_sig_snps$additive.a.wrtAA), ]$high.allele <- NA
all_sig_snps$high.allele.anc <- "M" # TEMP to make sterile geno designation easier
all_sig_snps[!is.na(all_sig_snps$additive.a.wrtDD) & all_sig_snps$additive.a.wrtDD > 0, ]$high.allele.anc <- "D"
all_sig_snps[is.na(all_sig_snps$additive.a.wrtDD), ]$high.allele.anc <- NA
all_sig_snps$high.geno <- all_sig_snps$high.allele
all_sig_snps[!is.na(all_sig_snps$d.a.ratio) & all_sig_snps$d.a.ratio > 1, ]$high.geno <- "H"
all_sig_snps$high.geno.anc <- all_sig_snps$high.allele.anc
all_sig_snps[!is.na(all_sig_snps$dom.allele) & !is.na(all_sig_snps$d.a.ratio) & all_sig_snps$d.a.ratio > 1, ]$high.geno.anc <- "H"


write.csv(all_sig_snps, file = "all_sig_snps_pheno_means20250317.csv")
all_sig_snps <- read_csv("all_sig_snps_pheno_means20250317.csv")

# SW sig snps
all_sig_snps_SW_bon <- all_sig_snps |> filter(P.sig.bon_SW == "Y" | add.P.sig.bon_SW == "Y" | dom.P.sig.bon_SW == "Y")
all_sig_snps_SW_meff <- all_sig_snps |> filter(P.sig.Meff_SW == "Y" | add.P.sig.Meff_SW == "Y" | dom.P.sig.Meff_SW == "Y")

write.csv(all_sig_snps_SW_meff, file = "all_sig_snps_pheno_means20250317_SW_meff.csv")
write.csv(all_sig_snps_SW_bon, file = "all_sig_snps_pheno_means20250317_SW_bon.csv")

# To get consensus sterile allele D/M for locus make composite col listing sterile allele for each SNP

all_sig_snps <- all_sig_snps[order(all_sig_snps$trait, all_sig_snps$chr.num, all_sig_snps$pos), ]


all_sig_loci$high.allele.anc <- NA
all_sig_loci$high.geno.anc <- NA
all_sig_snps_SW <- na.omit(all_sig_snps)
for (row in 1:nrow(all_sig_loci)) {
  snps <- all_sig_snps_SW[which(all_sig_snps_SW$sig.loc.ID == all_sig_loci$sig.loc.ID[row]), ]
  all_sig_loci$high.allele.anc[row] <- paste(snps[is.na(snps$high.allele.anc) == F, ]$high.allele.anc, collapse = ".")
  all_sig_loci$high.geno.anc[row] <- paste(snps[is.na(snps$high.geno.anc) == F, ]$high.geno.anc, collapse = ".")
}

all_sig_loci <- all_sig_loci[order(all_sig_loci$trait, all_sig_loci$chr.num, all_sig_loci$start.pos), ]
write.csv(all_sig_loci, file = "all_sig_loci_SW_20250317.csv", row.names = F)
save(all_sig_loci, file = "all_sig_loci_SW_20250317.Rdata")




# Loci #####
all_peak_snps <- all_sig_snps %>% inner_join(all_sig_loci %>% select(sig.loc.ID, peak.snp), by = (c("marker" = "peak.snp", "sig.loc.ID")))

##### minor allele ##



all_peak_snps$bins <- cut(all_peak_snps$degree_of_dominance,
  breaks = c(-Inf, seq(-1.95, 2, 0.05), Inf),
  labels = c(seq(-2, 2, 0.05))
)

sig.all_cat <- all_peak_snps %>%
  group_by(dominance_cat) %>%
  tally() %>%
  mutate(percent = (n / nrow(all_peak_snps)) * 100) %>%
  drop_na() %>%
  mutate(
    start = c(-Inf, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25), stop = c(-1.25, -0.75, -0.25, 0.25, 0.75, 1.25, Inf),
    middle = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
  )

sig.all_sum <- all_peak_snps %>%
  group_by(bins) %>%
  tally() %>%
  drop_na()



ggplot() +
  geom_segment(data = sig.all_sum, aes(
    x = as.numeric(as.character(bins)),
    xend = as.numeric(as.character(bins)), y = 0, yend = n
  ), linewidth = 5, alpha = 0.9) +
  scale_x_continuous(
    breaks = c(-Inf, -2, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.99, Inf),
    labels = c(-Inf, "<-2", -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, ">2", Inf)
  ) +
  geom_rect(
    data = sig.all_cat, aes(xmin = start, xmax = stop, fill = dominance_cat),
    ymin = -Inf, ymax = Inf, alpha = 0.2
  ) +
  geom_text(data = sig.all_cat, mapping = aes(x = middle, y = 1900, label = paste(round(percent, 2), "%")), vjust = -0.5, size = 5.5) +
  geom_text(data = sig.all_cat, mapping = aes(x = middle, y = 2000, label = dominance_cat), vjust = -0.5, size = 4.5) +
  labs(x = expression(paste("Dominance", italic("(d/|a|)"))), y = "Frequency", fill = "") +
  theme(legend.position = "none") +
  scale_fill_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 7))
ggsave("dominance_classified_loci_minor_allele20250317.pdf", width = 12)

# low allele
all_peak_snps_low_allele <- all_peak_snps %>% mutate(degree_of_dominance = ifelse(low.allele == "B", degree_of_dominance * -1, degree_of_dominance))
all_peak_snps_low_allele$bins <- cut(all_peak_snps_low_allele$degree_of_dominance,
  breaks = c(-Inf, seq(-1.95, 2, 0.05), Inf),
  labels = c(seq(-2, 2, 0.05))
)
all_peak_snps_low_allele$dominance_cat_low_allele <- cut(all_peak_snps_low_allele$degree_of_dominance,
  breaks = c(-Inf, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, Inf),
  labels = c("underdominant", "recessive", "partially recessive", "additive", "partially dominant", "dominant", "overdominant")
)

sig.all_cat_low_allele <- all_peak_snps_low_allele %>%
  group_by(dominance_cat_low_allele) %>%
  tally() %>%
  mutate(percent = (n / nrow(all_peak_snps_low_allele)) * 100) %>%
  drop_na() %>%
  mutate(
    start = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), stop = c(-0.75,-0.25, 0.25, 0.75, 1.25, Inf),
    middle = c(-1,-0.5, 0, 0.5, 1, 1.5)
  )

sig.all_sum_low_allele <- all_peak_snps_low_allele %>%
  group_by(bins) %>%
  tally() %>%
  drop_na()




ggplot() +
  geom_segment(data = sig.all_sum_low_allele, aes(
    x = as.numeric(as.character(bins)),
    xend = as.numeric(as.character(bins)), y = 0, yend = n
  ), linewidth = 5, alpha = 0.9) +
  scale_x_continuous(
    breaks = c(-Inf, -2, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.99, Inf),
    labels = c(-Inf, "<-2", -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, ">2", Inf)
  ) +
  geom_rect(
    data = sig.all_cat_low_allele, aes(xmin = start, xmax = stop, fill = dominance_cat_low_allele),
    ymin = -Inf, ymax = Inf, alpha = 0.2
  ) +
  geom_text(data = sig.all_cat_low_allele, mapping = aes(x = middle, y = 2300, label = paste(round(percent, 2), "%")), vjust = -0.5, size = 5.5) +
  geom_text(data = sig.all_cat_low_allele, mapping = aes(x = middle, y = 2400, label = dominance_cat_low_allele), vjust = -0.5, size = 4.5) +
  labs(x = expression(paste("Dominance", italic("(d/|a|)"))), y = "Frequency", fill = "") +
  theme(legend.position = "none") +
  scale_fill_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 7))
ggsave("dominance_classified_loci_low_allele20250317.pdf", width = 12)
ggsave("dominance_classified_loci_low_allele20250317.png", width = 12)

# all snps #####

##### minor allele ##
all_sig_snps$bins <- cut(all_sig_snps$degree_of_dominance,
  breaks = c(-Inf, seq(-1.95, 2, 0.05), Inf),
  labels = c(seq(-2, 2, 0.05))
)

sig.all_cat.snps <- all_sig_snps %>%
  group_by(dominance_cat) %>%
  tally() %>%
  mutate(percent = (n / nrow(all_sig_snps)) * 100) %>%
  drop_na() %>%
  mutate(
    start = c(-Inf, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25), stop = c(-1.25, -0.75, -0.25, 0.25, 0.75, 1.25, Inf),
    middle = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
  )

sig.all_sum.snps <- all_sig_snps %>%
  group_by(bins) %>%
  tally() %>%
  drop_na()




ggplot() +
  geom_segment(data = sig.all_sum.snps, aes(
    x = as.numeric(as.character(bins)),
    xend = as.numeric(as.character(bins)), y = 0, yend = n
  ), linewidth = 5, alpha = 0.9) +
  scale_x_continuous(
    breaks = c(-Inf, -2, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.99, Inf),
    labels = c(-Inf, "<-2", -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, ">2", Inf)
  ) +
  geom_rect(
    data = sig.all_cat.snps, aes(xmin = start, xmax = stop, fill = dominance_cat),
    ymin = -Inf, ymax = Inf, alpha = 0.2
  ) +
  geom_text(data = sig.all_cat.snps, mapping = aes(x = middle, y = 100, label = paste(round(percent, 2), "%")), vjust = -0.5, size = 5.5) +
  geom_text(data = sig.all_cat.snps, mapping = aes(x = middle, y = 95, label = dominance_cat), vjust = -0.5, size = 4.5) +
  labs(x = expression(paste("Dominance", italic("(d/|a|)"))), y = "Frequency", fill = "") +
  theme(legend.position = "none") +
  scale_fill_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 7))
ggsave("dominance_classified_snps_minor_allele20250317.pdf", width = 12)


##### low allele ##
all_sig_snps_low_allele <- all_sig_snps %>% mutate(degree_of_dominance = ifelse(low.allele == "B", degree_of_dominance * -1, degree_of_dominance))
all_sig_snps_low_allele$bins <- cut(all_sig_snps_low_allele$degree_of_dominance,
  breaks = c(-Inf, seq(-1.95, 2, 0.05), Inf),
  labels = c(seq(-2, 2, 0.05))
)
all_sig_snps_low_allele$cat <- cut(all_sig_snps_low_allele$degree_of_dominance,
  breaks = c(-Inf, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, Inf),
  labels = c("underdominant", "recessive", "partially recessive", "additive", "partially dominant", "dominant", "overdominant")
)

sig.all_cat_low_allele.snps <- all_sig_snps_low_allele %>%
  group_by(cat) %>%
  tally() %>%
  mutate(percent = (n / nrow(all_sig_snps_low_allele)) * 100) %>%
  drop_na() %>%
  mutate(
    start = c(-Inf, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25), stop = c(-1.25, -0.75, -0.25, 0.25, 0.75, 1.25, Inf),
    middle = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
  )

sig.all_sum_low_allele.snps <- all_sig_snps_low_allele %>%
  group_by(bins) %>%
  tally() %>%
  drop_na()




ggplot() +
  geom_segment(data = sig.all_sum_low_allele.snps, aes(
    x = as.numeric(as.character(bins)),
    xend = as.numeric(as.character(bins)), y = 0, yend = n
  ), linewidth = 5, alpha = 0.9) +
  scale_x_continuous(
    breaks = c(-Inf, -2, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.99, Inf),
    labels = c(-Inf, "<-2", -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, ">2", Inf)
  ) +
  geom_rect(
    data = sig.all_cat_low_allele.snps, aes(xmin = start, xmax = stop, fill = cat),
    ymin = -Inf, ymax = Inf, alpha = 0.2
  ) +
  geom_text(data = sig.all_cat_low_allele.snps, mapping = aes(x = middle, y = 100, label = paste(round(percent, 2), "%")), vjust = -0.5, size = 5.5) +
  geom_text(data = sig.all_cat_low_allele.snps, mapping = aes(x = middle, y = 95, label = cat), vjust = -0.5, size = 4.5) +
  labs(x = expression(paste("Dominance", italic("(d/|a|)"))), y = "Frequency", fill = "") +
  theme(legend.position = "none") +
  scale_fill_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 7))
ggsave("dominance_classified_snps_low_allele20250317.pdf", width = 12)
ggsave("dominance_classified_snps_low_allele20250317.png", width = 12)



#### ancestry #####
tot_inf_trait <- all_sig_snps %>%
  group_by(trait, high.allele.anc) %>%
  # drop_na() %>%
  tally()

ggplot(tot_inf_trait, aes(x = fct_reorder(trait, n), y = n, fill = high.allele.anc)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() +
  theme_test() +
  labs(y = "number of SNPs", x = "", fill = "Ancestry informative") +
  scale_fill_manual(values = c("purple4", "dark green", "grey"))
# scale_y_continuous(limits = c(0,35), expand = c(0, 0))

ggsave("high_allele_snps_ancestry_bar_plot20250317.pdf", height = 15, width = 12)

tot_inf_trait_loci <- all_peak_snps %>%
  group_by(trait, high.allele.anc) %>%
  # drop_na() %>%
  tally()

ggplot(tot_inf_trait_loci, aes(x = fct_reorder(trait, n), y = n, fill = high.allele.anc)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() +
  theme_test() +
  labs(y = "Number of loci", x = "", fill = "Origin high allele") +
  scale_fill_manual(values = c("purple4", "dark green"))
# scale_y_continuous(limits = c(0,35), expand = c(0, 0))
ggsave("high_allele_loci_ancestry_bar_plot20250317.pdf", width = 12, height = 12)


tot_inf_cat <- all_sig_snps_low_allele %>%
  group_by(trait, low.allele.anc, cat) %>%
  drop_na() %>%
  tally()

ggplot(tot_inf_cat, aes(x = fct_reorder(trait, n), y = n, fill = cat)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() +
  theme_test() +
  facet_wrap(~low.allele.anc, scales = "free_x") +
  scale_fill_d3()
ggsave("low_allele_ancestry_snps_facet20250317.pdf", width = 12, height = 15)


tot_inf_cat_loci <- all_peak_snps_low_allele %>%
  group_by(trait, low.allele.anc, dominance_cat) %>%
  drop_na() %>%
  tally() %>%
  ungroup() %>%
  group_by(trait) %>%
  mutate(sum = sum(n)) %>%
  ungroup()

ggplot(tot_inf_cat_loci, aes(x = fct_reorder(trait, sum), y = n, fill = dominance_cat)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() +
  theme_test() +
  facet_wrap(~low.allele.anc, scales = "free_x") +
  scale_fill_d3()
ggsave("low_allele_ancestry_loci_facet20250317.pdf", width = 12, height = 15)


tot_inf_sum <- all_sig_snps_low_allele %>%
  group_by(low.allele.anc) %>%
  tally()

pie <- ggplot(tot_inf_sum, aes(x = "", y = n, fill = factor(low.allele.anc))) +
  geom_bar(stat = "identity", width = 1, color = "black", alpha = 0.8) +
  coord_polar("y", start = 0) +
  theme_void() + # remove background, grid, numeric labels
  scale_fill_manual(values = c("purple4", "dark green"))

library(scales)
pie + geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%")),
  position = position_stack(vjust = 0.5)
) + labs(fill = "Low allele")


ggsave("low.allele_mus_dom_informative_markers_snps20250317.pdf")

tot_inf_sum_loci <- all_peak_snps_low_allele %>%
  group_by(low.allele.anc) %>%
  tally()

pie <- ggplot(tot_inf_sum_loci, aes(x = "", y = n, fill = factor(low.allele.anc))) +
  geom_bar(stat = "identity", width = 1, color = "black", alpha = 0.8) +
  coord_polar("y", start = 0) +
  theme_void() + # remove background, grid, numeric labels
  scale_fill_manual(values = c("purple4", "dark green", "honeydew3"))

library(scales)
pie + geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%")),
  position = position_stack(vjust = 0.5)
) + labs(fill = "Low allele")


ggsave("low.allele_mus_dom_informative_markers_loci20250317.pdf")

tot_inf_sum_loci <- all_peak_snps_low_allele %>%
  group_by(low.geno) %>%
  tally()

pie <- ggplot(tot_inf_sum_loci, aes(x = "", y = n, fill = factor(low.geno))) +
  geom_bar(stat = "identity", width = 1, color = "white", alpha = 0.8) +
  coord_polar("y", start = 0) +
  theme_void() + # remove background, grid, numeric labels
  scale_fill_d3()

library(scales)
pie + geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%")),
  position = position_stack(vjust = 0.5)
) + labs(fill = "Low genotype")


ggsave("low.allele_genotype_min_maj_H_loci20250317.pdf")



tot_inf_sum_loci <- all_sig_snps_low_allele %>%
  group_by(low.geno) %>%
  tally()

pie <- ggplot(tot_inf_sum_loci, aes(x = "", y = n, fill = factor(low.geno))) +
  geom_bar(stat = "identity", width = 1, color = "white", alpha = 0.8) +
  coord_polar("y", start = 0) +
  theme_void() + # remove background, grid, numeric labels
  scale_fill_d3()

library(scales)
pie + geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%")),
  position = position_stack(vjust = 0.5)
) + labs(fill = "Low genotype")


ggsave("low.allele_genotype_min_maj_H_snps20250317.pdf")


# all snps SW #####

##### minor allele ##
all_sig_snps_SW_meff$bins <- cut(all_sig_snps_SW_meff$degree_of_dominance,
  breaks = c(-Inf, seq(-1.95, 2, 0.05), Inf),
  labels = c(seq(-2, 2, 0.05))
)

sig.all_cat.snps_SW <- all_sig_snps_SW_meff %>%
  group_by(dominance_cat) %>%
  tally() %>%
  mutate(percent = (n / nrow(all_sig_snps_SW_meff)) * 100) %>%
  drop_na() %>%
  mutate(
    start = c(-Inf, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25), stop = c(-1.25, -0.75, -0.25, 0.25, 0.75, 1.25, Inf),
    middle = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
  )

sig.all_sum.snps_SW <- all_sig_snps_SW_meff %>%
  group_by(bins) %>%
  tally() %>%
  drop_na()




ggplot() +
  geom_segment(data = sig.all_sum.snps_SW, aes(
    x = as.numeric(as.character(bins)),
    xend = as.numeric(as.character(bins)), y = 0, yend = n
  ), linewidth = 5, alpha = 0.9) +
  scale_x_continuous(
    breaks = c(-Inf, -2, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.99, Inf),
    labels = c(-Inf, "<-2", -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, ">2", Inf)
  ) +
  geom_rect(
    data = sig.all_cat.snps_SW, aes(xmin = start, xmax = stop, fill = dominance_cat),
    ymin = -Inf, ymax = Inf, alpha = 0.2
  ) +
  geom_text(data = sig.all_cat.snps_SW, mapping = aes(
    x = middle, y = 700,
    label = paste(round(percent, 2), "%")
  ), vjust = -0.5, size = 5.5) +
  geom_text(data = sig.all_cat.snps_SW, mapping = aes(
    x = middle, y = 680,
    label = dominance_cat
  ), vjust = -0.5, size = 4.5) +
  labs(x = expression(paste("Dominance", italic("(d/|a|)"))), y = "Frequency", fill = "") +
  theme(legend.position = "none") +
  scale_fill_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 7))
ggsave("dominance_classified_snps_minor_allele_SW_20250317.pdf", width = 12)


##### low allele ##
all_sig_snps_low_allele_SW <- all_sig_snps_SW_meff %>% mutate(degree_of_dominance = ifelse(low.allele == "B", degree_of_dominance * -1, degree_of_dominance))
all_sig_snps_low_allele_SW$bins <- cut(all_sig_snps_low_allele_SW$degree_of_dominance,
  breaks = c(-Inf, seq(-1.95, 2, 0.05), Inf),
  labels = c(seq(-2, 2, 0.05))
)
all_sig_snps_low_allele_SW$cat <- cut(all_sig_snps_low_allele_SW$degree_of_dominance,
  breaks = c(-Inf, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, Inf),
  labels = c("underdominant", "recessive", "partially recessive", "additive", "partially dominant", "dominant", "overdominant")
)

sig.all_cat_low_allele.snps_SW <- all_sig_snps_low_allele_SW %>%
  group_by(cat) %>%
  tally() %>%
  mutate(percent = (n / nrow(all_sig_snps_low_allele_SW)) * 100) %>%
  drop_na() %>%
  mutate(
    start = c(-Inf, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25), stop = c(-1.25, -0.75, -0.25, 0.25, 0.75, 1.25, Inf),
    middle = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
  )

sig.all_sum_low_allele.snps_SW <- all_sig_snps_low_allele_SW %>%
  group_by(bins) %>%
  tally() %>%
  drop_na()




ggplot() +
  geom_segment(data = sig.all_sum_low_allele.snps_SW, aes(
    x = as.numeric(as.character(bins)),
    xend = as.numeric(as.character(bins)), y = 0, yend = n
  ), linewidth = 5, alpha = 0.9) +
  scale_x_continuous(
    breaks = c(-Inf, -2, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.99, Inf),
    labels = c(-Inf, "<-2", -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, ">2", Inf)
  ) +
  geom_rect(
    data = sig.all_cat_low_allele.snps_SW, aes(xmin = start, xmax = stop, fill = cat),
    ymin = -Inf, ymax = Inf, alpha = 0.2
  ) +
  geom_text(data = sig.all_cat_low_allele.snps_SW, mapping = aes(x = middle, y = 900, label = paste(round(percent, 2), "%")), vjust = -0.5, size = 5.5) +
  geom_text(data = sig.all_cat_low_allele.snps_SW, mapping = aes(x = middle, y = 850, label = cat), vjust = -0.5, size = 4.5) +
  labs(x = expression(paste("Dominance", italic("(d/|a|)"))), y = "Frequency", fill = "") +
  theme(legend.position = "none") +
  scale_fill_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 7))
ggsave("dominance_classified_snps_low_allele_SW_20250317.pdf", width = 12)
ggsave("dominance_classified_snps_low_allele_SW_20250317.png", width = 12)




##### high allele ##
all_sig_snps_high_allele_SW <- all_sig_snps_SW_meff %>%
  mutate(degree_of_dominance = ifelse(high.allele == "B", degree_of_dominance * -1, degree_of_dominance))
all_sig_snps_high_allele_SW$bins <- cut(all_sig_snps_high_allele_SW$degree_of_dominance,
  breaks = c(-Inf, seq(-1.95, 2, 0.05), Inf),
  labels = c(seq(-2, 2, 0.05))
)
all_sig_snps_high_allele_SW$cat <- cut(all_sig_snps_high_allele_SW$degree_of_dominance,
  breaks = c(-Inf, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, Inf),
  labels = c("underdominant", "recessive", "partially recessive", "additive", "partially dominant", "dominant", "overdominant")
)

sig.all_cat_high_allele.snps_SW <- all_sig_snps_high_allele_SW %>%
  group_by(cat) %>%
  tally() %>%
  mutate(percent = (n / nrow(all_sig_snps_high_allele_SW)) * 100) %>%
  drop_na() %>%
  mutate(
    start = c(-Inf, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25), stop = c(-1.25, -0.75, -0.25, 0.25, 0.75, 1.25, Inf),
    middle = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
  )

sig.all_sum_high_allele.snps_SW <- all_sig_snps_high_allele_SW %>%
  group_by(bins) %>%
  tally() %>%
  drop_na()




ggplot() +
  geom_segment(data = sig.all_sum_high_allele.snps_SW, aes(
    x = as.numeric(as.character(bins)),
    xend = as.numeric(as.character(bins)), y = 0, yend = n
  ), linewidth = 5, alpha = 0.9) +
  scale_x_continuous(
    breaks = c(-Inf, -2, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.99, Inf),
    labels = c(-Inf, "<-2", -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, ">2", Inf)
  ) +
  geom_rect(
    data = sig.all_cat_high_allele.snps_SW, aes(xmin = start, xmax = stop, fill = cat),
    ymin = -Inf, ymax = Inf, alpha = 0.2
  ) +
  geom_text(data = sig.all_cat_high_allele.snps_SW, mapping = aes(x = middle, y = 900, label = paste(round(percent, 2), "%")), vjust = -0.5, size = 5.5) +
  geom_text(data = sig.all_cat_high_allele.snps_SW, mapping = aes(x = middle, y = 850, label = cat), vjust = -0.5, size = 4.5) +
  labs(x = expression(paste("Dominance", italic("(d/|a|)"))), y = "Frequency", fill = "") +
  theme(legend.position = "none") +
  scale_fill_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 7))
ggsave("dominance_classified_snps_high_allele_SW_20250317.pdf", width = 12)
ggsave("dominance_classified_snps_high_allele_SW_20250317.png", width = 12)



#### ancestry #####
tot_inf_trait_SW <- all_sig_snps_SW_meff %>%
  group_by(trait, high.allele.anc) %>%
  # drop_na() %>%
  tally()

ggplot(tot_inf_trait_SW, aes(x = fct_reorder(trait, n), y = n, fill = high.allele.anc)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() +
  theme_test() +
  labs(y = "number of SNPs", x = "", fill = "Ancestry informative") +
  scale_fill_manual(values = c("purple4", "dark green", "grey"))
# scale_y_continuous(limits = c(0,35), expand = c(0, 0))

ggsave("high_allele_snps_ancestry_bar_plot20250317_SW.pdf", height = 15, width = 12)

tot_inf_trait_loci_SW <- all_peak_snps %>%
  group_by(trait, high.allele.anc) %>%
  # drop_na() %>%
  tally()

ggplot(tot_inf_trait_loci_SW, aes(x = fct_reorder(trait, n), y = n, fill = high.allele.anc)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() +
  theme_test() +
  labs(y = "Number of loci", x = "", fill = "Origin high allele") +
  scale_fill_manual(values = c("purple4", "dark green"))
# scale_y_continuous(limits = c(0,35), expand = c(0, 0))
ggsave("high_allele_loci_ancestry_bar_plot_SW20250317.pdf", width = 12, height = 12)

tot_inf_trait_loci_SW_geno <- all_peak_snps %>%
  group_by(trait, high.geno.anc) %>%
  # drop_na() %>%
  tally()
ggplot(tot_inf_trait_loci_SW_geno, aes(x = fct_reorder(trait, n), y = n, fill = high.geno.anc)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() +
  theme_test() +
  labs(y = "Number of loci", x = "", fill = "Origin high allele") +
  scale_fill_manual(values = c("purple4",  "red","dark green"))
# scale_y_continuous(limits = c(0,35), expand = c(0, 0))
ggsave("high_allele_geno_loci_trait_ancestry_bar_plot_SW20250317.pdf", width = 12, height = 12)

library(RColorBrewer)
n <- 61
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

plot_high_geno <- ggplot(na.omit(tot_inf_trait_loci_SW_geno), aes(x =high.geno.anc, y=n, fill = fct_reorder(trait, n))) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
scale_fill_manual(values=col_vector)

# Using the cowplot package
legend <- cowplot::get_legend(plot_high_geno)



plot_high_geno + theme(legend.position="none")

ggsave("high_allele_geno_loci_ancestry_bar_plot_SW20250317.pdf", width = 12, height = 12)
grid::grid.draw(legend)


pdf("high_allele_geno_loci_ancestry_bar_plot_SW20250317_legend.pdf", width=20, height=12)
grid::grid.newpage()
grid::grid.draw(legend)
dev.off()


tot_inf_cat_SW <- all_sig_snps_low_allele_SW %>%
  group_by(trait, low.allele.anc, cat) %>%
  drop_na() %>%
  tally()

ggplot(tot_inf_cat_SW, aes(x = fct_reorder(trait, n), y = n, fill = cat)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() +
  theme_test() +
  facet_wrap(~low.allele.anc, scales = "free_x") +
  scale_fill_d3()
ggsave("low_allele_ancestry_snps_SW_facet2041115.pdf", width = 12, height = 15)

all_peak_snps_low_allele_SW <- all_peak_snps_low_allele|> filter(P.sig.Meff_SW == "Y" | add.P.sig.Meff_SW == "Y" | dom.P.sig.Meff_SW == "Y")

tot_inf_cat_loci_SW <- all_peak_snps_low_allele_SW %>%
  group_by(trait, low.allele.anc, dominance_cat) %>%
  drop_na() %>%
  tally() %>%
  ungroup() %>%
  group_by(trait) %>%
  mutate(sum = sum(n)) %>%
  ungroup()

ggplot(tot_inf_cat_loci_SW, aes(x = fct_reorder(trait, sum), y = n, fill = dominance_cat)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() +
  theme_test() +
  facet_wrap(~low.allele.anc, scales = "free_x") +
  scale_fill_d3()
ggsave("low_allele_ancestry_loci_facet20250317_SW.pdf", width = 12, height = 15)


tot_inf_sum_SW <- all_sig_snps_low_allele_SW %>%
  group_by(low.allele.anc) %>%
  tally()

pie <- ggplot(tot_inf_sum_SW, aes(x = "", y = n, fill = factor(low.allele.anc))) +
  geom_bar(stat = "identity", width = 1, color = "black", alpha = 0.8) +
  coord_polar("y", start = 0) +
  theme_void() + # remove background, grid, numeric labels
  scale_fill_manual(values = c("purple4", "dark green"))

library(scales)
pie + geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%")),
  position = position_stack(vjust = 0.5)
) + labs(fill = "Low allele")


ggsave("low.allele_mus_dom_informative_markers_snps_SW_20250317.pdf")

tot_inf_sum_loci <- all_peak_snps_low_allele_SW %>%
  group_by(low.allele.anc) %>%
  tally()

pie <- ggplot(tot_inf_sum_loci, aes(x = "", y = n, fill = factor(low.allele.anc))) +
  geom_bar(stat = "identity", width = 1, color = "black", alpha = 0.8) +
  coord_polar("y", start = 0) +
  theme_void() + # remove background, grid, numeric labels
  scale_fill_manual(values = c("purple4", "dark green", "honeydew3"))

library(scales)
pie + geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%")),
  position = position_stack(vjust = 0.5)
) + labs(fill = "Low allele")


ggsave("low.allele_mus_dom_informative_markers_loci20250317.pdf")

tot_inf_sum_loci_SW <- all_peak_snps_low_allele_SW %>%
  group_by(low.geno) %>%
  tally()

pie <- ggplot(tot_inf_sum_loci_SW, aes(x = "", y = n, fill = factor(low.geno))) +
  geom_bar(stat = "identity", width = 1, color = "white", alpha = 0.8) +
  coord_polar("y", start = 0) +
  theme_void() + # remove background, grid, numeric labels
  scale_fill_d3()

library(scales)
pie + geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%")),
  position = position_stack(vjust = 0.5)
) + labs(fill = "Low genotype")


ggsave("low.allele_genotype_min_maj_H_loci20250317_SW.pdf")



tot_inf_sum_loci_SW <- all_sig_snps_low_allele_SW %>%
  group_by(low.geno) %>%
  tally()

pie <- ggplot(tot_inf_sum_loci_SW, aes(x = "", y = n, fill = factor(low.geno))) +
  geom_bar(stat = "identity", width = 1, color = "white", alpha = 0.8) +
  coord_polar("y", start = 0) +
  theme_void() + # remove background, grid, numeric labels
  scale_fill_d3()

library(scales)
pie + geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%")),
  position = position_stack(vjust = 0.5)
) + labs(fill = "Low genotype")


ggsave("low.allele_genotype_min_maj_H_snps_SW20250317.pdf")

# high allele
tot_inf_sum_loci_SW <- all_sig_snps_high_allele_SW %>%
  group_by(high.geno) %>%
  tally()

pie <- ggplot(tot_inf_sum_loci_SW, aes(x = "", y = n, fill = factor(high.geno))) +
  geom_bar(stat = "identity", width = 1, color = "white", alpha = 0.8) +
  coord_polar("y", start = 0) +
  theme_void() + # remove background, grid, numeric labels
  scale_fill_d3()

library(scales)
pie + geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%")),
  position = position_stack(vjust = 0.5)
) + labs(fill = "High genotype")


ggsave("high.allele_genotype_min_maj_H_snps_SW20250317.pdf")


tot_inf_sum_SW <- all_sig_snps_high_allele_SW %>%
  group_by(high.allele.anc) %>%
  tally()

pie <- ggplot(tot_inf_sum_SW, aes(x = "", y = n, fill = factor(high.allele.anc))) +
  geom_bar(stat = "identity", width = 1, color = "black", alpha = 0.8) +
  coord_polar("y", start = 0) +
  theme_void() + # remove background, grid, numeric labels
  scale_fill_manual(values = c("purple4", "dark green"))

library(scales)
pie + geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%")),
  position = position_stack(vjust = 0.5)
) + labs(fill = "High allele")


ggsave("high.allele_mus_dom_informative_markers_snps_SW_20250317.pdf")
