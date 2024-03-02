# comparison with milbrandt paper https://www.nature.com/articles/s41593-021-01005-1

# libraries  ----
library(Seurat)
library(BPCells)
library(SeuratObject)
library(tidyverse)
library(conflicted)
library(qs)
library(scMisc)
library(writexl)
library(readxl)
library(RColorBrewer)
library(homologene)
library(ggvenn)

# general settings  ----
options(warn = 0)
future::plan("multicore", workers = 6)
conflicts_prefer(base::setdiff)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

# seurat map to reference datasets ----
pns_sn_sciatic_milbrandt <- qs::qread("/home/mischko/Documents/beruf/forschung/scRNA_reference/pns_atlas_milbrandt/pns_sn_sciatic_GSE182098.qs", nthreads = 4)

# markers ----
# find markers function
findMarkers <- function(ident1, ident2 = NULL, object, only_pos, min_pct, logfc_threshold, assay = assay) {
  result <- Seurat::FindMarkers(object, ident.1 = ident1, ident.2 = ident2, min.pct = min_pct, logfc.threshold = logfc_threshold, only.pos = only_pos, assay = assay) |>
    tibble::rownames_to_column("gene") |>
    dplyr::filter(p_val_adj < 0.05) |>
    dplyr::relocate(gene, avg_log2FC, p_val, p_val_adj) |>
    dplyr::arrange(desc(avg_log2FC))
  return(result)
}

topmarkers <-
  lapply(
    unique(sc_merge@misc$cluster_order),
    function(x) {
      message("Processing cluster ", x)
      try(findMarkers(ident1 = x, object = sc_merge, only_pos = TRUE, min_pct = 0.1, logfc_threshold = 0.25, assay = "RNA"))
    }
  )

# or read in previously calculated markers for our dataset
topmarkers <-
    lapply(
        sc_merge@misc$cluster_order,
        function(x) {
            read_excel(file.path("results", "de", "topmarkers_final.xlsx"), sheet = x)
        }
    ) |>
    setNames(sc_merge@misc$cluster_order)

sc_merge$cluster <- as.character(sc_merge$cluster)
sc_merge$cluster[sc_merge$cluster == "periC1"] <- "periC"
sc_merge$cluster[sc_merge$cluster == "periC2"] <- "periC"
sc_merge$cluster[sc_merge$cluster == "periC3"] <- "periC"

topmarkers_milbrandt <-
    lapply(
        unique(pns_sn_sciatic_milbrandt$cluster),
        function(x) {
            message("Processing cluster ", x)
            try(findMarkers(ident1 = x, object = pns_sn_sciatic_milbrandt, only_pos = TRUE, min_pct = 0.1, logfc_threshold = 0.25, assay = "RNA"))
        }
    ) |>
    setNames(unique(pns_sn_sciatic_milbrandt$cluster))

write_xlsx(topmarkers_milbrandt, file.path("results", "de", "topmarkers_milbrandt.xlsx"))

# make names consistent
names(topmarkers_milbrandt)[names(topmarkers_milbrandt) == "EnFibro"] <- "endoC"
names(topmarkers_milbrandt)[names(topmarkers_milbrandt) == "EpC"] <- "epiC"
names(topmarkers_milbrandt)[names(topmarkers_milbrandt) == "PnC"] <- "periC"

Idents(sc_merge) <- sc_merge$cluster
topmarkers[["periC"]] <- findMarkers(
    ident1 = "periC",
    object = sc_merge,
    only_pos = TRUE,
    min_pct = 0.1,
    logfc_threshold = 0.25,
    assay = "RNA"
)


vennPlot <- function(cluster) {
    top_rodent <- topmarkers_milbrandt[[cluster]]$gene |>
        homologene::mouse2human(db = homologene::homologeneData2) |>
        # dplyr::slice(1:25) |>
        pull(humanGene)
    top_human <- topmarkers[[cluster]]$gene[1:25]
    # top_human <- topmarkers[[cluster]]$gene
    # top_human <- markers_xenium$transcript[markers_xenium$cluster == cluster]
    plot <- ggvenn(
        list(
            rodent = top_rodent,
            human = top_human
        ),
        fill_color = brewer.pal(name = "Set2", n = 3),
        show_elements = TRUE,
        # show_elements = FALSE,
        label_sep = "\n",
        text_size = 1,
    )
    ggsave(file.path("results", "venn", paste0(cluster, "_milbrandt_top25.pdf")),
        width = 5, height = 5,
        plot = plot
    )
}

lapply(
    c("mySC", "nmSC", "endoC", "epiC", "periC"),
    FUN = vennPlot
)

vennPlot("mySC")
vennPlot("nmSC")


markers_xenium <- read_xlsx(file.path("lookup", "xenium_list_final.xlsx"))

test1 <- 
markers_xenium |>
    dplyr::filter(cluster == "mySC") |>
    pull(transcript)

test2 <- topmarkers$mySC$gene[1:25]

union(test1, test2)
