# enrichment analysis using enrichR

#load libraries
library(enrichR)
library(qs)
library(Seurat)
library(BPCells)
library(SeuratObject)
library(tidyverse)
library(patchwork)
library(conflicted)
library(scMisc)
library(readxl)
library(writexl)

# general settings  ----
options(warn = 0)
conflicts_prefer(base::setdiff)

# load preprocessed data ----
ic <- qread(file.path("objects", "ic.qs"), nthread = 4)

# get enrichment dbs ----
dbs <- enrichR::listEnrichrDbs()
dbs <- c("TF_Perturbations_Followed_by_Expression", "Transcription_Factor_PPIs", "WikiPathways_2023_Human", "KEGG_2021_Human", "Reactome_2022", "Panther_2016", "NCI-Nature_2016", "GO_Biological_Process_2023", "GO_Molecular_Function_2023")

# get de of macro ----
macro_de_top <- 
    lapply(c("Macro2", "Macro17", "Macro18"), 
        FUN = function(x) {
            read_excel(file.path("results", "de", "topmarkers_ic.xlsx"), sheet = x) |>
                dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.001)
                }) |>
        setNames(c("Macro2", "Macro17", "Macro18")) 


# get de of main cluster ---
cluster_de <- c("repairSC", "mySC", "nmSC", "PC2", "ven_capEC1", "artEC")

cluster_de_top_pos <-
    lapply(cluster_de,
        FUN = function(x) {
            read_excel(file.path("results", "de", "pnp_ctrl_pseudobulk_sig.xlsx"), sheet = x) |>
            dplyr::filter(avg_logFC > 2)
        } 
    ) |>
    setNames(cluster_de)

cluster_de_top_neg <-
    lapply(cluster_de,
        FUN = function(x) {
            read_excel(file.path("results", "de", "pnp_ctrl_pseudobulk_sig.xlsx"), sheet = x) |>
            dplyr::filter(avg_logFC < 2)
        } 
    ) |>
    setNames(cluster_de)

# perform enrichment with macro ----
macro_enrichr <-
    lapply(c("Macro2", "Macro17", "Macro18"),
        FUN = function(x) {
            enrichr(macro_de_top[[x]]$gene, dbs)
        }
    ) |>
    setNames(c("Macro2", "Macro17", "Macro18"))


lapply(
    c("Macro2", "Macro17", "Macro18"),
    FUN = function(x) {
        write_xlsx(macro_enrichr[[x]], file.path("results", "enrichr", paste0("enrichr_", x, ".xlsx")))
    })

# perform enrichment with main clusters ----
clusters_enrichr_pos <-
    lapply(cluster_de,
        FUN = function(x) {
            enrichr(cluster_de_top_pos[[x]]$gene, dbs)
        }
    ) |>
    setNames(cluster_de)

clusters_enrichr_neg <-
    lapply(cluster_de,
        FUN = function(x) {
            enrichr(cluster_de_top_neg[[x]]$gene, dbs)
        }
    ) |>
    setNames(cluster_de)

lapply(
    cluster_de,
    FUN = function(x) {
        write_xlsx(clusters_enrichr_pos[[x]], file.path("results", "enrichr", paste0("enrichr_pos_", x, ".xlsx")))
    }
)

lapply(
    cluster_de,
    FUN = function(x) {
        write_xlsx(clusters_enrichr_neg[[x]], file.path("results", "enrichr", paste0("enrichr_neg_", x, ".xlsx")))
    }
)

# function to plot enrichment results ----
plotEnrichrFun <- function(filename, sheet, width, height) {
    colors <- RColorBrewer::brewer.pal(5, "Set2")
    enrichr <- readxl::read_excel(file.path("results", "enrichr", glue::glue("enrichr_{filename}.xlsx")), sheet = sheet) |>
        dplyr::slice_min(order_by = Adjusted.P.value, n = 10, with_ties = FALSE) |>
        tidyr::separate(Overlap, into = c("overlap1", "overlap2")) |> # separate overlap in two columns
        dplyr::mutate(overlap = as.numeric(overlap1) / as.numeric(overlap2)) |> # calculcate overlap
        ggplot(aes(y = reorder(Term, -log10(Adjusted.P.value)), x = -log10(Adjusted.P.value))) +
        geom_col(fill = scales::hue_pal()(5)[1]) +
        labs(
            x = "-Log10 Adjusted P value",
            y = ""
        ) +
        theme_classic() +
        theme(legend.position = "none")
    ggsave(file.path("results", "enrichr", glue::glue("barplot_enrichr_{filename}_{sheet}.pdf")), width = width, height = height)
}

# plot enrichment of macro ----
lapply(
    c("Macro2", "Macro17", "Macro18"),
    FUN = function(x) {
        plotEnrichrFun(x, sheet = "GO_Biological_Process_2023", width = 7, height = 3)
    })

# plot enrichment of main clusters ----
lapply(
    c(
        paste0("pos_", cluster_de),
        paste0("neg_", cluster_de)
    ),
    FUN = function(x) {
        plotEnrichrFun(x, sheet = "GO_Biological_Process_2023", width = 7, height = 3)
    }
)
