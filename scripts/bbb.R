# compare with BBB reference dataset https://www.nature.com/articles/s41586-022-04521-7

# libraries
library(Seurat)
library(qs)
library(tidyverse)
library(BiocParallel)
library(BPCells)
library(ProjecTILs)
library(scMisc)

# load data
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"))
pals::cols25()

artEC <- subset(sc_merge, subset  = cluster %in% c("artEC"))

DimPlot(artEC, label = TRUE, reduction = "umap.scvi.full")

# load reference
bbb <- qs::qread("/home/mischko/Documents/beruf/forschung/scRNA_reference/bbb/bbb_vascular_projectil.qs")
rosmap <- qs::qread("/home/mischko/Documents/beruf/forschung/scRNA_reference/bbb/rosmap_projectil.qs")


rosmap_plot <-
 DimPlot(rosmap, group.by = "cluster", raster = FALSE, pt.size = .5, alpha = .3, cols = pals::cols25(), label = TRUE) +
  theme_rect() + 
  NoLegend() +
  ggtitle("")
ggsave(plot = rosmap_plot, file.path("results", "projectil", "rosmap_reference.png"), width = 8, height = 8)



bbb_plot <-
 DimPlot(bbb, group.by = "cluster", raster = FALSE, pt.size = .5, alpha = .3, cols = pals::cols25(), label = TRUE) +
  theme_rect() + 
  NoLegend() +
  ggtitle("")
ggsave(plot = bbb_plot, file.path("results", "projectil", "bbb_reference.png"), width = 8, height = 8)

# function to plot projection from reference using projectil
PlotProjectil <- function(query_cluster, ref) {
    ref_name <- deparse(substitute(ref))
    query <- subset(sc_merge, subset = cluster %in% query_cluster)
    projectil <- Run.ProjecTILs(query = query, ref = ref)
    p1 <- plot.projection(
        ref = ref,
        query = projectil,
        linesize = 0.2,
        pointsize = 0.2,
        cols = pals::cols25()
    ) +
        scMisc::theme_rect() +
        ggtitle(query_cluster)
    ggsave(
        filename = file.path("results", "projectil", paste0(ref_name, "_", query_cluster, ".pdf")),
        plot = p1,
        width = 8,
        height = 8
    )
}

# select clusters to project separately
sel_clusters <- c("artEC", "ven_capEC1", "ven_capEC2", "venEC", "PC1", "PC2", "VSMC")

# project selected clusters on rosmap
lapply(
    sel_clusters,
    FUN = purrr::possibly(
        function(x) {
            PlotProjectil(
                query_cluster = x,
                ref = rosmap
            )
        }
    )
)

# project selected clusters on bbb
lapply(
    sel_clusters,
    FUN = purrr::possibly(
        function(x) {
            PlotProjectil(
                query_cluster = x,
                ref = bbb
            )
        }
    )
)
