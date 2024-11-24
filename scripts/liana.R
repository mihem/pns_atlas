# ===============================================================================
# Cell-Cell Communication Analysis using LIANA
# ===============================================================================
# Purpose: Analyze cell-cell communication patterns focusing on CXCL14 interactions
# using multiple methods (natmi, connectome, logfc, sca) via the LIANA framework.
#
# Input:
#   - Preprocessed single-cell data (sc_merge.qs)
# Output:
#   - LIANA results object (liana_results.qs)
#   - CXCL14 interaction dotplot (cxcl14_interactions_dotplot.pdf)
# ===============================================================================

# Load required libraries ----
library(liana)
library(tidyverse)
library(Seurat)
library(BPCells)
library(qs)
library(liana)

# Load and preprocess data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

# Convert counts and data (necessary for LIANA) to dgCMatrix
sc_merge$RNA$counts <- as(sc_merge$RNA$counts, "dgCMatrix")
sc_merge$RNA$data <- as(sc_merge$RNA$data, "dgCMatrix")

# Interaction analysis setup ----
# Get consensus interactions and filter for CXCL-related pairs
consensus_interactions <- select_resource("Consensus")$Consensus

cxcl_pairs <- consensus_interactions |>
  filter(grepl("CXCL", source_genesymbol) | grepl("CXCL", target_genesymbol))

# LIANA analysis ----
# Perform cell-cell communication analysis using multiple methods
# using only CXCL interactions
liana_results <- liana_wrap(
  sc_merge,
  method = c("natmi", "connectome", "logfc", "sca"),
  resource = "custom",
  external_resource = cxcl_pairs,
  expression_pct = 0.1,
  parallel = TRUE,
  interactions = cxcl14_pairs
)

liana_prep

qs::qsave(liana_results, file.path("objects", "liana_results.qs"))

# Results processing and visualization ----
liana_results$natmi |>
  dplyr::filter(ligand == "CXCL14") |>
  print(width = Inf)

liana_results_aggregate <-
  liana_results |>
  liana_aggregate()

liana_results_aggregate |>
  dplyr::filter(ligand.complex == "CXCL14") |>
  print(width = Inf)

liana_plot_cxcl14 <-
  liana_results_aggregate |>
  dplyr::filter(ligand.complex == "CXCL14") |>
  liana_dotplot()

# Save results ----
# Generate and save dotplot visualization
ggsave(
  file.path("results", "liana", "cxcl14_interactions_dotplot.pdf"),
  liana_plot_cxcl14,
  width = 15,
  height = 5
)

# Visualize CXCR4 expression ----
dplot_cxcr4 <-
  DotPlot(sc_merge, features = c("CXCR4"), dot.min = 0.01) +
  viridis::scale_color_viridis(option = "viridis") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic")) +
  xlab(NULL) +
  ylab(NULL) +
  scale_size(range = c(0, 5))

ggsave(
  file.path("results", "liana", "cxcr4_dotplot.pdf"),
  dplot_cxcr4,
  width = 3.5,
  height = 5
)
