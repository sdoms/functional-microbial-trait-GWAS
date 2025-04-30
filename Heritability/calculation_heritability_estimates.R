# Set working directory to location of heritability mapping results
setwd("~/Documents/Research/Experiments/shotgun_hybrids/mapping_results/heritability/")

# Load required packages
library(lme4qtl) # For linear mixed models with kinship matrices
library(BEDMatrix) # For efficient access to PLINK binary genotype data
library(MASS) # General statistical methods
library(parallel) # Parallel computing tools
library(gtools) # Utility functions for data manipulation
library(lme4) # Fitting linear and generalized linear mixed-effects models
library(plyr) # Tools for splitting, applying, and combining data
library(lmerTest) # Add p-values to lme4 mixed models
library(car) # Companion to Applied Regression (for modeling)
library(EnvStats) # Environmental statistics tools
library(RLRsim) # Restricted likelihood ratio tests
library(tidyverse) # Collection of data science packages
library(argyle) # Tools for SNP array data (used here for kinship calculation)

# Path to GEMMA executable (used to calculate the kinship matrix)
gemma <- "/opt/homebrew/bin/gemma"

##### Load phenotype data #####
pheno.in <- read.csv("../../../Final_QTL_mapping/SterilityPhenosG2sMap20210218.csv", header = TRUE)
rownames(pheno.in) <- pheno.in$mouse_name
pheno.in$id <- pheno.in$mouse_name

# Load functional trait data (e.g. metagenomic modules, EC numbers, etc.)
modules <- read.csv("../all_functional_traits.csv")

# Merge phenotype and functional module data by mouse name
pheno <- modules %>%
  inner_join(pheno.in, by = c(X = "mouse_name")) |>
  as.data.frame()
rownames(pheno) <- pheno$X

##### Load genotype data #####
geno <- BEDMatrix("../../../Final_QTL_mapping/Cleaning_snps/clean_f2")
rownames(geno) <- substr(rownames(geno), 5, 19) # Clean mouse names
rownames(geno) <- gsub(x = rownames(geno), pattern = "\\/", replacement = ".")
colnames(geno) <- substr(colnames(geno), 1, nchar(colnames(geno)) - 2)

# Keep only individuals with phenotype data
indi_all <- rownames(geno)
load("../../../Final_QTL_mapping/Cleaning_snps/clean_snps.Rdata") # Loads `snps` object
individuals <- pheno$id[pheno$id %in% indi_all]
pheno <- pheno[individuals, ]
geno <- geno[individuals, ]
indi_all <- rownames(geno)

# Check matching individuals between genotype and phenotype
stopifnot(all(rownames(geno) == rownames(pheno)))

##### Generate kinship matrix (centered) using GEMMA #####
# Read PLINK genotype data and convert to numeric genotype matrix
F2 <- read.plink("../../../Final_QTL_mapping/Cleaning_snps/clean_f2")
F2 <- argyle::recode(F2, "relative")
colnames(F2) <- gsub(x = colnames(F2), pattern = "\\/", replacement = ".")

F2 <- as.matrix(as.data.frame(F2))
F2_part <- apply(F2[, 7:ncol(F2)], 2, as.numeric)
F2_part <- F2_part[, individuals] # Subset to matched individuals
F2 <- cbind(snps[, c(1, 5, 6)], F2_part) # Recombine SNP info and genotype matrix

# Write genotype and dummy phenotype files for GEMMA
write.table(F2, "geno.txt", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(data.frame(rep(1, nrow(pheno))), "pheno.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Run GEMMA to calculate centered kinship matrix
system(paste0(gemma, " -g geno.txt -gk 1 -p pheno.txt -o kinship_GMM"))

##### Estimate heritability using centered kinship matrix #####
results <- data.frame()
for (i in 2:(ncol(pheno) - 1)) { # Skip first column (likely mouse ID)
  trait <- colnames(pheno)[i]
  covariate <- "" # Modify if covariates are available

  cat("Running model for ", trait, ".\n")

  # Load kinship matrix from GEMMA output
  kinship <- read.table("./output/kinship_GMM.cXX.txt") |> as.matrix()
  rownames(kinship) <- colnames(F2)[4:ncol(F2)]
  colnames(kinship) <- colnames(F2)[4:ncol(F2)]

  df <- data.frame(pheno)
  df$trait <- df[, trait]
  if (covariate == "") {
    null_model <- relmatLmer(trait ~ (1 | id), df, relmat = list(id = kinship))
  } else {
    df$covar <- df[, covariate]
    null_model <- relmatLmer(trait ~ covar + (1 | id), df, relmat = list(id = kinship))
  }

  # Extract variance components (heritability estimate)
  h2 <- VarProp(null_model)
  pval <- lmerTest::ranova(null_model)

  results[trait, "h2_id"] <- h2[1, 6]
  results[trait, "h2_residual"] <- h2[2, 6]
  results[trait, "pval"] <- pval$`Pr(>Chisq)`[2]
}

# Adjust for multiple testing (Benjamini-Hochberg)
results$p.adj <- p.adjust(results$pval, method = "BH")
write.csv(results, "heritability_estimates_all_functional_traits_centered_matrix.csv")

##### Repeat with standardized kinship matrix #####
# Run GEMMA with standardized genetic relationship matrix
system(paste0(gemma, " -g geno.txt -gk 2 -p pheno.txt -o kinship_GMM"))

results <- data.frame()
for (i in 2:(ncol(pheno) - 1)) {
  trait <- colnames(pheno)[i]
  covariate <- ""

  cat("Running model for ", trait, ".\n")

  # Load standardized kinship matrix
  kinship <- read.table("./output/kinship_GMM.sXX.txt") |> as.matrix()
  rownames(kinship) <- colnames(F2)[4:ncol(F2)]
  colnames(kinship) <- colnames(F2)[4:ncol(F2)]

  df <- data.frame(pheno)
  df$trait <- df[, trait]
  if (covariate == "") {
    null_model <- relmatLmer(trait ~ (1 | id), df, relmat = list(id = kinship))
  } else {
    df$covar <- df[, covariate]
    null_model <- relmatLmer(trait ~ covar + (1 | id), df, relmat = list(id = kinship))
  }

  h2 <- VarProp(null_model)
  pval <- lmerTest::ranova(null_model)

  results[trait, "h2_id"] <- h2[1, 6]
  results[trait, "h2_residual"] <- h2[2, 6]
  results[trait, "pval"] <- pval$`Pr(>Chisq)`[2]
}

results$p.adj <- p.adjust(results$pval, method = "BH")
write.csv(results, "heritability_estimates_all_functional_traits_standardised_matrix.csv")
