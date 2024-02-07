#libraries ---
library(Seurat)
library(tidyverse)
library(factoextra)
library(FactoMineR)
library(viridis)

# general settings ----
conflicts_prefer(base::setdiff)
conflicts_prefer(base::as.data.frame)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"))
ic <- qs::qread(file.path("objects", "ic.qs"))
ic_cidp_vn_ctrl_ciap <- subset(ic, level2 %in% c("CIDP", "VN", "CTRL", "CIAP"))

table(ic$level2)
table(ic_cidp_vn_ctrl_ciap$level2)

# function to plot pca of Seurat abundancies
pcaSeurat <- function(object, label1, label2, label3) {
  object_parse <- deparse(substitute(object))
  cl_size <-
    base::as.data.frame.matrix(table(object@meta.data[[label1]], object@meta.data[[label2]])) |>
    t()

  colnames(cl_size) <- levels(object@meta.data[[label1]])

  pca_result <-FactoMineR::PCA(cl_size, scale.unit = TRUE, ncp = 30, graph = FALSE)
  factoextra::fviz_eig(pca_result, addlabels = TRUE, ylim = c(0,50), ncp = 7)
  ggsave(paste0("./results/pca/pca_", object_parse, "_", label3, "_eigen.pdf"))

pca_var_plot <-
  factoextra::fviz_pca_var(
    pca_result,
    col.var = "contrib",
    gradient.cols = viridis::viridis(100),
    repel = TRUE,
   select.var = list(contrib = 10)) + 
    labs(title = "")+
    theme_classic()

  lookup_pre <-
    data.frame(
      label1 = object@meta.data[[label2]],
      label3 = object@meta.data[[label3]]
    ) |>
    distinct()

  lookup <-
    data.frame(label1 = rownames(cl_size)) |>
    left_join(lookup_pre)


  pca_plot_ind <-
    factoextra::fviz_pca_ind(
      pca_result,
      pointsize = 3,
      pointshape = 21,
      fill = "#E7B800",
      col.ind = "black",
      palette = "Set2",
      axes.linetype = "solid"
    )

  pca_ggplot_ind <-
    ggpubr::ggpar(
    pca_plot_ind,
    title = "",
    xlab = "PC1",
    ylab = "PC2",
    ggtheme = theme_bw() +
      theme(axis.title.x = element_text(size=15),
            axis.title.y = element_text(size=15),
            plot.title = element_text(size=25))) +
    theme_classic() +
    coord_fixed()

   pca_plot_group <-
    factoextra::fviz_pca_ind(
      pca_result,
      pointsize = 3,
      pointshape = 21,
      geom.ind = "point",
      fill.ind = lookup$label3,
      col.ind = "black",
      palette = "Set2",
      addEllipses = FALSE,
      ellipse.type = "confidence",
      legend.title = "group",
      axes.linetype = "solid",
      alpha.ind = .5
    )

  pca_ggplot_group <-
    ggpubr::ggpar(
      pca_plot_group,
      title = "",
      xlab = "PC1",
      ylab = "PC2",
      ggtheme = theme_classic() +
        theme(axis.title.x = element_text(size=15),
              axis.title.y = element_text(size=15),
              plot.title = element_text(size=25))) +
    # geom_point(size = 2, alpha = .5, shape = 16) + 
    # scale_fill_manual(values = object@misc$level2_cols) + 
    theme(
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        aspect.ratio = 1
    ) 


if(is.numeric(lookup$label3)) {
  pca_ggplot_group <- pca_ggplot_group + viridis::scale_fill_viridis(option = "magma")
}

  pca_plots <- patchwork::wrap_plots(pca_var_plot, pca_ggplot_ind, pca_ggplot_group, ncol = 3)
  ggsave(paste0("./results/pca/pca_", object_parse, "_", label3,  ".pdf"), width = 18, height = 6,
         plot = pca_plots)

}

debugonce(pcaSeurat)

# run pca abundance function
pcaSeurat(
  object = sc_merge,
  label1 = "cluster",
  label2 = "sample",
  label3 = "level1"
)

pcaSeurat(
  object = ic_cidp_vn_ctrl_ciap,
  label1 = "ic_cluster",
  label2 = "sample",
  label3 = "level2"
)

pcaSeurat(
  object = sc_merge,
  label1 = "cluster",
  label2 = "sample",
  label3 = "level2"
)

pcaSeurat(
  object = sc_merge_cidp_vn_ctrl_ciap,
  label1 = "cluster",
  label2 = "sample",
  label3 = "level2"
)

pcaSeurat(
  object = sc_merge,
  label1 = "cluster",
  label2 = "sample",
  label3 = "center"
)

pcaSeurat(
  object = sc_merge,
  label1 = "cluster",
  label2 = "sample",
  label3 = "log_axon_normal"
)

ic_cidp_vn_ctrl_ciap$incat <- as.numeric(ic_cidp_vn_ctrl_ciap$incat)

pcaSeurat(
  object = ic_cidp_vn_ctrl_ciap,
  label1 = "ic_cluster",
  label2 = "sample",
  label3 = "incat"
)

# without scDamage
sc_merge_no_damageSC <- subset(sc_merge, subset = cluster %in% c("damageSC"), invert = TRUE)
sc_merge_no_damageSC$cluster <- droplevels(sc_merge_no_damageSC$cluster)

pcaSeurat(
  object = sc_merge_no_damageSC,
  label1 = "cluster",
  label2 = "sample",
  label3 = "level2"
)

# only intraneuronal cells
sc_merge_neuronal <- subset(sc_merge, subset = cluster %in% c("LEC", "artEC", "venEC", "epiC", "periC3", "periC2", "periC1", "Adipo", "damageSC"), invert = TRUE)

sc_merge_neuronal$cluster <- droplevels(sc_merge_neuronal$cluster)

dplyr::count(sc_merge_neuronal@meta.data, cluster)

pcaSeurat(
  object = sc_merge_neuronal,
  label1 = "cluster",
  label2 = "sample",
  label3 = "level2"
)

# only sc cluster
sc_merge_sc <- subset(sc_merge, subset = cluster %in% c("nmSC", "mySC", "repairSC", "damageSC"))

sc_merge_sc$cluster <- droplevels(sc_merge_sc$cluster)

dplyr::count(sc_merge_sc@meta.data, cluster)

pcaSeurat(
  object = sc_merge_sc,
  label1 = "cluster",
  label2 = "sample",
  label3 = "level2"
)

# NMF not superior to PCA
# library(NMF)

# cl_size_sc <-
#   as.data.frame.matrix(table(sc_merge@meta.data[["cluster"]], sc_merge@meta.data[["sample"]])) |>
#   t()
# colnames(cl_size_sc) <- levels(sc_merge@meta.data[["cluster"]])

# nmf_sc <- NMF::nmf(cl_size_sc, 2, method = "Frobenius")

# coldata <- readr::read_csv(file.path("lookup", "sample_lookup.csv")) 

# nmf_df <-
#  as.data.frame(NMF::basis(nmf_sc)) |>
#  mutate(level2 = coldata$level2) |>
#  tibble::rownames_to_column("sample")


# ggplot(nmf_df, aes(x = V1, V2, color = sample)) +
#   geom_point()
