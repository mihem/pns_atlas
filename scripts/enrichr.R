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

#macro 18
macro18_de_top <- 
    read_excel(file.path("results", "de", "topmarkers_ic.xlsx"), sheet = "Macro18") |>
    dplyr::filter(avg_log2FC > 2 & p_val_adj < 0.001)

macro2_de_top <- 
    read_excel(file.path("results", "de", "topmarkers_ic.xlsx"), sheet = "Macro2") |>
    dplyr::filter(avg_log2FC > 2 & p_val_adj < 0.001)

macro17_de_top <- 
    read_excel(file.path("results", "de", "topmarkers_ic.xlsx"), sheet = "Macro17") |>
    dplyr::filter(avg_log2FC > 2 & p_val_adj < 0.001)

macro_de_top <- 
    lapply(c("Macro2", "Macro17", "Macro18"), 
        FUN = function(x) {
            read_excel(file.path("results", "de", "topmarkers_ic.xlsx"), sheet = x) |>
                dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.001)
                }) |>
        setNames(c("Macro2", "Macro17", "Macro18")) 

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

lapply(
    c("Macro2", "Macro17", "Macro18"),
    FUN = function(x) {
        plotEnrichrFun(x, sheet = "GO_Biological_Process_2023", width = 7, height = 3)
    })
