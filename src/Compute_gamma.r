#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript Compute_gamma.r <alignment_file> <output_csv>")
}

alignment_file <- args[1]
output_file <- args[2]
alpha <- 8

library(ape)
library(phangorn)
library(tools)

alignment <- read.FASTA(alignment_file, type = "AA")
alignment_matrix <- as.character(as.matrix(alignment))

var_sites <- apply(alignment_matrix, 2, function(column) {
  unique_chars <- unique(column[column != "X"])
  length(unique_chars) > 1
})

alignment_matrix_var <- alignment_matrix[, var_sites, drop = FALSE]
phyDat_alignment_var <- phyDat(alignment_matrix_var, type = "AA")

dist_matrix_var <- dist.ml(phyDat_alignment_var, model = "JTT")
tree_nj_var <- nj(dist_matrix_var)

fit_var <- pml(tree_nj_var, data = phyDat_alignment_var)
fit_gamma <- update(fit_var, k = alpha)
fit_gamma <- optim.pml(fit_gamma, model = "JTT", optGamma = TRUE)

site_indices <- attr(phyDat_alignment_var, "index")
expanded_likelihoods <- fit_gamma$siteLik[site_indices]

normalized_likelihoods <- (expanded_likelihoods - min(expanded_likelihoods)) /
  (max(expanded_likelihoods) - min(expanded_likelihoods))

site_rates <- 1 - normalized_likelihoods
site_rates <- site_rates * max(fit_gamma$g)

breaks <- c(0, sort(fit_gamma$g))
gamma_categories <- cut(
  site_rates,
  breaks = breaks,
  labels = 1:(length(breaks) - 1),
  include.lowest = TRUE
)

site_table <- data.frame(
  Position = which(var_sites),
  SiteLik = expanded_likelihoods,
  SubstitutionRate = site_rates,
  GammaCategory = as.numeric(gamma_categories)
)

write.csv(site_table, output_file, row.names = FALSE)
