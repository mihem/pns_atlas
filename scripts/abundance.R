#libraries ---
library(Seurat)
library(tidyverse)
library(scMisc)
library(qs)
library(pheatmap)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)
ic <- qs::qread(file.path("objects", "ic.qs"), nthread = 4)

sample_lookup <- 
  readr::read_csv(file.path("lookup", "sample_lookup.csv")) |>
  janitor::clean_names() |>
  mutate(age_calc = lubridate::time_length(difftime(nerve_date, birth_date), "years")) |>
  mutate(age_calc = floor(age_calc)) |>
  mutate(age = coalesce(age_calc, age)) |>
  dplyr::select(-age_calc) |>
  mutate(level0 = if_else(level1 == "CTRL", "CTRL", "PNP"))

  # abundance table
scMisc::abundanceTbl(sc_merge, "cluster", "sample")
scMisc::abundanceTbl(sc_merge, "cluster", "level2")
scMisc::abundanceTbl(ic, "ic_cluster", "sample")
scMisc::abundanceTbl(ic, "ic_cluster", "level2")

# stacked plots --
scMisc::stackedPlot(
  object = sc_merge,
  x_axis = "sample",
  y_axis = "cluster",
  x_order = unique(sc_merge$sample),
  y_order = sc_merge@misc$cluster_order,
  color = sc_merge@misc$cluster_col,
  width = 10
)

scMisc::stackedPlot(
  object = sc_merge,
  x_axis = "level1",
  y_axis = "cluster",
  x_order = unique(sc_merge$level1),
  y_order = sc_merge@misc$cluster_order,
  color = sc_merge@misc$cluster_col,
  width = 4
)

scMisc::stackedPlot(
  object = sc_merge,
  x_axis = "level2",
  y_axis = "cluster",
  x_order = unique(sc_merge$level2),
  y_order = sc_merge@misc$cluster_order,
  color = sc_merge@misc$cluster_col,
  width = 5
)

scMisc::stackedPlot(
  object = sc_merge,
  x_axis = "incat",
  y_axis = "cluster",
  x_order = as.character(1:6),
  y_order = sc_merge@misc$cluster_order,
  color = sc_merge@misc$cluster_col,
  width = 5
)

# immune cells ---
scMisc::stackedPlot(
  object = ic,
  x_axis = "level2",
  y_axis = "ic_cluster",
  x_order = unique(ic$level2),
  y_order = ic@misc$ic_cluster_order,
  color = ic@misc$ic_cluster_col,
  width = 5
)


str(sc_merge@meta.data)

# propeller abundance analysis ---

# PNP vs CTRL
propeller_PNP_CTRL <-
  scMisc::propellerCalc(
    seu_obj1 = sc_merge,
    condition1 = "PNP",
    condition2 = "CTRL",
    cluster_col = "cluster",
    meta_col = "level0",
    lookup = sample_lookup,
    sample_col = "sample",
    formula = "~0 + level0",
    min_cells = 30
  )

scMisc::plotPropeller(
  data = propeller_PNP_CTRL,
  color = sc_merge@misc$cluster_col,
  filename = "PNP_CTRL",
  FDR = 0.1
)

scMisc::dotplotPropeller(
    data = propeller_PNP_CTRL,
    color = sc_merge@misc$cluster_col,
    filename = "PNP_CTRL",
)

# immune cells
propeller_PNP_CTRL_ic <-
  scMisc::propellerCalc(
    seu_obj1 = ic,
    condition1 = "PNP",
    condition2 = "CTRL",
    cluster_col = "ic_cluster",
    meta_col = "level0",
    lookup = sample_lookup,
    sample_col = "sample",
    formula = "~0 + level0",
    min_cells = 30
  )

scMisc::plotPropeller(
  data = propeller_PNP_CTRL_ic,
  color = ic@misc$ic_cluster_col,
  filename = "PNP_CTRL_ic",
  FDR = 0.1
)

scMisc::dotplotPropeller(
    data = propeller_PNP_CTRL_ic,
    color = ic@misc$ic_cluster_col,
    filename = "PNP_CTRL_ic",
    height = 6
)

# CIDP vs CTRL
propeller_CIDP_CTRL <-
  scMisc::propellerCalc(
    seu_obj1 = sc_merge,
    condition1 = "CIDP",
    condition2 = "CTRL",
    cluster_col = "cluster",
    meta_col = "level2",
    lookup = sample_lookup,
    sample_col = "sample",
    formula = "~0 + level2",
    # formula = "~0 + level2 + sex + age + center",
    min_cells = 30
  )

scMisc::plotPropeller(
  data = propeller_CIDP_CTRL,
  color = sc_merge@misc$cluster_col,
#   filename = "CIDP_CTRL_sex_age_center",
  filename = "CIDP_CTRL",
  FDR = 0.1
)

scMisc::dotplotPropeller(
    data = propeller_CIDP_CTRL,
    color = sc_merge@misc$cluster_col,
    filename = "CIDP_CTRL",
)

# immune cells ----
propeller_CIDP_CTRL_ic <-
  scMisc::propellerCalc(
    seu_obj1 = ic,
    condition1 = "CIDP",
    condition2 = "CTRL",
    cluster_col = "ic_cluster",
    meta_col = "level2",
    lookup = sample_lookup,
    sample_col = "sample",
    formula = "~0 + level2",
    # formula = "~0 + level2 + sex + age + center",
    min_cells = 30
  )

scMisc::plotPropeller(
  data = propeller_CIDP_CTRL_ic,
  color = ic@misc$ic_cluster_col,
  filename = "CIDP_CTRL_ic",
  FDR = 0.1
)

scMisc::dotplotPropeller(
    data = propeller_CIDP_CTRL_ic,
    color = ic@misc$ic_cluster_col,
    filename = "CIDP_CTRL_ic",
    height = 6
)

#  VN vs CTRL
propeller_VN_CTRL <-
    scMisc::propellerCalc(
        seu_obj1 = sc_merge,
        condition1 = "VN",
        condition2 = "CTRL",
        cluster_col = "cluster",
        meta_col = "level2",
        lookup = sample_lookup,
        sample_col = "sample",
        formula = "~0 + level2",
        min_cells = 30
    )

scMisc::plotPropeller(
  data = propeller_VN_CTRL,
  color = sc_merge@misc$cluster_col,
  filename = "VN_CTRL",
  FDR = 0.1
)

scMisc::dotplotPropeller(
  data = propeller_VN_CTRL,
  color = sc_merge@misc$cluster_col,
  filename = "VN_CTRL"
)

# immune cells
propeller_VN_CTRL_ic <-
    scMisc::propellerCalc(
        seu_obj1 = ic,
        condition1 = "VN",
        condition2 = "CTRL",
        cluster_col = "ic_cluster",
        meta_col = "level2",
        lookup = sample_lookup,
        sample_col = "sample",
        formula = "~0 + level2",
        min_cells = 30
    )

scMisc::plotPropeller(
  data = propeller_VN_CTRL_ic,
  color = ic@misc$ic_cluster_col,
  filename = "VN_CTRL_ic",
  FDR = 0.1
)

scMisc::dotplotPropeller(
  data = propeller_VN_CTRL_ic,
  color = ic@misc$ic_cluster_col,
  filename = "VN_CTRL_ic",
  height = 6
)

#  CIAP vs CTRL
propeller_CIAP_CTRL <-
    scMisc::propellerCalc(
        seu_obj1 = sc_merge,
        condition1 = "CIAP",
        condition2 = "CTRL",
        cluster_col = "cluster",
        meta_col = "level2",
        lookup = sample_lookup,
        sample_col = "sample",
        formula = "~0 + level2",
        min_cells = 30
    )

scMisc::plotPropeller(
  data = propeller_CIAP_CTRL,
  color = sc_merge@misc$cluster_col,
  filename = "CIAP_CTRL",
  FDR = 0.1
)

scMisc::dotplotPropeller(
  data = propeller_CIAP_CTRL,
  color = sc_merge@misc$cluster_col,
  filename = "CIAP_CTRL"
)

# immune cells
propeller_CIAP_CTRL_ic <-
    scMisc::propellerCalc(
        seu_obj1 = ic,
        condition1 = "CIAP",
        condition2 = "CTRL",
        cluster_col = "ic_cluster",
        meta_col = "level2",
        lookup = sample_lookup,
        sample_col = "sample",
        formula = "~0 + level2",
        min_cells = 30
    )

scMisc::plotPropeller(
  data = propeller_CIAP_CTRL_ic,
  color = ic@misc$ic_cluster_col,
  filename = "CIAP_CTRL_ic",
  FDR = 0.1
)

scMisc::dotplotPropeller(
  data = propeller_CIAP_CTRL_ic,
  color = ic@misc$ic_cluster_col,
  filename = "CIAP_CTRL_ic",
  height = 6
)

# abundance of mrVI groups ----
mrvi_lookup <- read_csv(file.path("lookup", "mrvi_lookup.csv"))

sc_merge@meta.data <-
    sc_merge@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(mrvi_lookup, by = "sample") |>
    tibble::column_to_rownames(var = "barcode")

vn_cidp_ciap_ctrl <- subset(sc_merge, level2 %in% c("VN", "CIDP", "CIAP", "CTRL"))

scMisc::abBoxPlot(
  object = vn_cidp_ciap_ctrl,
  cluster_idents = "cluster",
  sample = "sample",
  cluster_order = vn_cidp_ciap_ctrl@misc$cluster_order,
  group_by =  "mrvi_cluster",
  group_order = paste0("cl", 1:5),
  color = pals::cols25()
)

phmap_mrvi_cluster <-
  table(vn_cidp_ciap_ctrl$cluster, vn_cidp_ciap_ctrl$mrvi_cluster) |>
  pheatmap(
    scale = "column",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = viridis::magma(100),
    cellwidth = 10,
    cellheight = 10,
    treeheight_row = 15,
    treeheight_col = 15,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    border_color = NA,
    cutree_rows = 5,
    main = "mrVI cluster"
  )

pdf(file.path("results", "abundance", "mrvi_heatmap_abundance_cluster.pdf"), width = 5, height = 7)
print(phmap_mrvi_cluster)
dev.off()


# g ratio and axon diameter boxplots ----
g_ratio_axon_diameter_mrvi <-
  g_ratio |>
  group_by(sample) |>
  mutate(
    g_ratio = mean(g_ratio),
    axon_diameter = mean(axon_diameter),
  ) |>
  ungroup() |>
  distinct() |>
  left_join(mrvi_lookup, join_by(sample)) |>
  dplyr::filter(!is.na(mrvi_cluster))  

# axon diameter
axon_diameter_mrvi_stats <- scMisc:::compStat(x_var = "axon_diameter", group = "mrvi_cluster", data = g_ratio_axon_diameter_mrvi, paired = FALSE)

axon_diameter_mrvi_plot <-
  g_ratio_axon_diameter_mrvi |>
  ggplot(aes(x = mrvi_cluster, y = axon_diameter)) +
  ggsignif::geom_signif(comparisons = axon_diameter_mrvi_stats$comparisons, annotation = axon_diameter_mrvi_stats$annotation, textsize = 5, step_increase = 0.05, vjust = 0.7)  +
  geom_boxplot(aes(fill = mrvi_cluster)) +
  geom_point() +
  theme_bw() +
  xlab("") +
  ylab("") +
  ggtitle("axon diameter") +
  scale_fill_manual(values = pals::cols25()) +
  theme(legend.position = "none")
  
ggsave(file.path("results", "abundance", "boxplot_axon_diameter_mrvi.pdf"), plot = axon_diameter_mrvi_plot, width = 5, height = 5)

# axon counts
axon_count_mrvi <-
  axon_count_mean |>
  left_join(mrvi_lookup, join_by(sample)) |>
  dplyr::filter(!is.na(mrvi_cluster))

axon_count_mrvi_stats <- scMisc:::compStat(x_var = "log_axon_normal", group = "mrvi_cluster", data = axon_count_mrvi, paired = FALSE)

axon_count_mrvi_plot <-
  axon_count_mrvi |>
  ggplot(aes(x = mrvi_cluster, y = axon_normal)) +
  ggsignif::geom_signif(comparisons = axon_count_mrvi_stats$comparisons, annotation = axon_count_mrvi_stats$annotation, textsize = 5, step_increase = 0.05, vjust = 0.7)  +
  geom_boxplot(aes(fill = mrvi_cluster)) +
  geom_point() +
  theme_bw() +
  xlab("") +
  ylab("") +
  ggtitle("axon count") +
  scale_fill_manual(values = pals::cols25()) +
  theme(legend.position = "none")
  
ggsave(file.path("results", "abundance", "boxplot_axon_count_mrvi.pdf"), plot = axon_count_mrvi_plot, width = 5, height = 5)

# incat -----
incat_mrvi <-
  sample_lookup |>
  left_join(mrvi_lookup, join_by(sample)) |>
  dplyr::filter(!is.na(mrvi_cluster)) |>
  mutate(incat = as.numeric(incat))

incat_mrvi_stats <- scMisc:::compStat(x_var = "incat", group = "mrvi_cluster", data = incat_mrvi, paired = FALSE)

incat_mrvi_plot <-
  incat_mrvi |>
  ggplot(aes(x = mrvi_cluster, y = incat)) +
  ggsignif::geom_signif(comparisons = incat_mrvi_stats$comparisons, annotation = incat_mrvi_stats$annotation, textsize = 5, step_increase = 0.05, vjust = 0.7)  +
  geom_boxplot(aes(fill = mrvi_cluster)) +
  geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  xlab("") +
  ylab("") +
  ggtitle("INCAT score") +
  scale_fill_manual(values = pals::cols25()) +
  theme(legend.position = "none")
  
ggsave(file.path("results", "abundance", "boxplot_incat_mrvi.pdf"), plot = incat_mrvi_plot, width = 5, height = 5)
