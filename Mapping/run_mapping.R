#########################################################################################################################
##                                Mapping Script for KEGG Identifiers                                                  ##
## This script performs association mapping using a linear mixed model with fixed additive and dominance SNP effects. ##
## Kinship is included as a random effect to account for relatedness.                                                  ##
## For each trait (e.g., KEGG ortholog abundances), it outputs GWAS statistics per chromosome.                        ##
#########################################################################################################################

# Set working directory and library paths
setwd("/home/doms/glm")
.libPaths("/home/doms/R/x86_64-pc-linux-gnu-library/3.5")

# Load required libraries
library(lme4qtl) # Linear mixed model with custom relatedness matrices
library(BEDMatrix) # Efficient genotype data handling
library(MASS)
library(parallel) # Parallel processing
library(gtools)
library(lme4)
library(plyr)
library(lmerTest) # P-values in mixed models
library(car) # ANOVA functions

# Parse input arguments
args <- commandArgs(TRUE)
tx <- as.integer(args[1]) # Index of the trait to be analyzed
input_file <- args[2] # Path to input file containing traits

# Create output folder based on input file name
output_folder <- gsub(input_file, pattern = ".csv", replacement = "_out")
output_folder <- gsub(output_folder, pattern = "input/", replacement = "")
dir.create(paste0("./out/", output_folder), showWarnings = FALSE)

## ---------------------------------------------------------------
## 1. Load Phenotype and Trait Data                             --
## ---------------------------------------------------------------
cat("Reading in phenotypes and covariates. \n")
pheno <- read.csv("./input/Phenotypes_new.csv", sep = ";", header = TRUE)
rownames(pheno) <- pheno$Mouse_Name
pheno$id <- pheno$Mouse_Name

taxa <- read.csv(input_file)
rownames(taxa) <- taxa$Identifier
taxa <- taxa[, -c(1)] # Remove identifier column
taxa <- t(taxa) # Transpose so individuals are rows

# Merge phenotypes with traits
da <- merge(taxa, pheno, by = "row.names")
rownames(da) <- da$Row.names
da$Row.names <- NULL

## ---------------------------------------------------------------
## 2. Load and Process Genotypes                                --
## ---------------------------------------------------------------
cat("Reading in genotypes. \n")
geno <- BEDMatrix("input/clean_f2")
rownames(geno) <- substr(rownames(geno), 4, 18)
rownames(geno) <- gsub("/", ".", rownames(geno)) # Sanitize row names
colnames(geno) <- substr(colnames(geno), 1, nchar(colnames(geno)) - 2)

# Load SNP metadata
load("input/clean_snps.Rdata")

# Remove mismatched individuals to preserve genotype order
geno <- geno[-c(73, 231), ]
indi_all <- rownames(geno)
data <- da[indi_all, ] # Reorder phenotype data to match genotype

stopifnot(rownames(data) == rownames(geno))

## ---------------------------------------------------------------
## 3. Create Additive and Dominance Genotype Matrices          --
## ---------------------------------------------------------------

# Additive: homozygous minor = -1, heterozygous = 0, homozygous major = 1
add.mat <- ifelse(geno[, ] == 0, 1, ifelse(geno[, ] == 1, 0, ifelse(geno[, ] == 2, -1, NA)))

# Dominance: heterozygous = 1, others = 0
dom.mat <- ifelse(abs(geno[, ]) == 2, 0, ifelse(geno[, ] == 1, 1, ifelse(geno[, ] == 0, 0, NA)))
rownames(dom.mat) <- rownames(add.mat)
colnames(dom.mat) <- colnames(add.mat)

## ---------------------------------------------------------------
## 4. Select Trait to Map                                      --
## ---------------------------------------------------------------
tax <- colnames(taxa)[tx]
taxa2 <- data[tax]
colnames(data)[which(colnames(data) == tax)] <- "tax"
cat("Running model for", tax, ".\n")

## ---------------------------------------------------------------
## 5. Run Association Model per Chromosome                     --
## ---------------------------------------------------------------

gwasResults <- data.frame()

for (chr in 1:19) {
  cat("Running model on chromosome", chr, ".\n")

  # Load chromosome-specific kinship matrix
  kinship <- as.matrix(read.table(paste0("./kinship/kinship_chr", chr, ".cXX.txt")))
  rownames(kinship) <- indi_all
  colnames(kinship) <- indi_all

  # Subset SNPs on this chromosome
  marker_chr <- snps[which(snps$chr == chr), 1]
  gts <- add.mat[indi_all, marker_chr]

  # Filter low MAF and non-missing
  sub <- colnames(gts)[
    colMeans(gts, na.rm = TRUE) / 2 > 0.025 &
      colMeans(gts, na.rm = TRUE) / 2 < 0.975 &
      !is.na(colMeans(gts, na.rm = TRUE))
  ]
  gts <- gts[, sub]

  # Prepare output data frame for this chromosome
  out <- data.frame(snps[sub, 1:6],
    tax = NA, n = NA, AA = NA, AB = NA, BB = NA,
    add.Beta = NA, add.StdErr = NA, add.T = NA,
    dom.Beta = NA, dom.StdErr = NA, dom.T = NA,
    P = NA, add.P = NA, dom.P = NA
  )
  out$tax <- tax

  # Fit null model once for full dataset
  df <- data.frame(data)
  null_model <- relmatLmer(tax ~ (1 | mating.pair) + (1 | id), df, relmat = list(id = kinship))
  size <- nrow(df)

  # Fit full model for each SNP
  f <- mclapply(as.list(sub), function(snp) {
    df <- data.frame(data, ad = add.mat[, snp], dom = dom.mat[, snp])
    df <- df[complete.cases(df[, c("tax", "ad", "dom")]), ]
    if (nrow(df) != size) {
      null_model <- relmatLmer(tax ~ (1 | mating.pair) + (1 | id), df, relmat = list(id = kinship))
    }
    model <- relmatLmer(tax ~ ad + dom + (1 | mating.pair) + (1 | id), df, relmat = list(id = kinship))

    res <- c(
      nrow(df),
      table(factor(df$ad, levels = c(-1, 0, 1))),
      tryCatch(summary(model)$coefficients[2, ], error = function(x) rep(NA, 3)),
      tryCatch(summary(model)$coefficients[3, ], error = function(x) rep(NA, 3)),
      tryCatch(anova(null_model, model)[2, 8], error = function(x) NA),
      tryCatch(Anova(model)[, 3], error = function(x) rep(NA, 2))
    )
    names(res) <- c(
      "n", "AA", "AB", "BB",
      "add.Beta", "add.StdErr", "add.T",
      "dom.Beta", "dom.StdErr", "dom.T",
      "P", "add.P", "dom.P"
    )
    return(res)
  }, mc.cores = getOption("mc.cores", 20))

  # Combine SNP-wise results
  out[, c(
    "n", "AA", "AB", "BB", "add.Beta", "add.StdErr", "add.T",
    "dom.Beta", "dom.StdErr", "dom.T", "P", "add.P", "dom.P"
  )] <- data.frame(do.call(rbind, f))

  out$index <- 1:nrow(out)
  gwasResults <- rbind(gwasResults, out)

  # Save results per chromosome
  saveRDS(out, paste0("./out/", output_folder, "/", tax, "_chr_", chr, "with_add_dom.rds"))
  head(out)
}

# Save combined results across chromosomes
saveRDS(gwasResults, paste0("./out/", tax, "_with_add_dom.rds"))
save(gwasResults, file = paste0("./out/", tax, "_with_add_dom.RData"))
