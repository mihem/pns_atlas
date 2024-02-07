# libraries  ----
library(Seurat)
library(BPCells)
library(SeuratObject)
library(tidyverse)
library(writexl)
library(patchwork)
library(conflicted)
library(qs)
library(Polychrome)
library(pals)
library(scMisc)
library(limma)
library(DESeq2)
library(viridis)
library(Libra)
library(readxl)

# general settings  ----
options(warn = 0)
options(Seurat.object.assay.version = "v5")
future::plan("multicore", workers = 6)
my_cols_25 <- pals::cols25()
my_cols_50 <- unname(Polychrome::createPalette(50, pals::cols25()))
conflicts_prefer(base::setdiff)
conflicts_prefer(base::as.data.frame)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)


#subset seurat object
sc_merge_cidp_vn_ctrl_ciap <- subset(sc_merge, subset = level2 %in% c("CIDP", "VN", "CTRL", "CIAP"))


# pseudobulk expression 
bulk <- AverageExpression(
    sc_merge,
    method = "aggregate",
    return.seurat = FALSE,
    slot = "counts",
    assays = "RNA",
    group.by = "sample"
)

bulk_cidp_vn_ctrl_ciap <- AverageExpression(
    sc_merge_cidp_vn_ctrl_ciap,
    method = "aggregate",
    return.seurat = FALSE,
    slot = "counts",
    assays = "RNA",
    group.by = "sample"
)


# pseudobulk MDS with DESeq ---
sc_merge$log_axon_normal <- log(sc_merge$axon_normal)

# get sample metadata
axon_count_lookup <-
    tibble(sample = sc_merge$sample, axon_normal = sc_merge$axon_normal, log_axon_normal = sc_merge$log_axon_normal) |>
    distinct() |>
    mutate()

coldata <- readr::read_csv(file.path("lookup", "sample_lookup.csv")) |>
    left_join(axon_count_lookup, by = "sample")

# create DESeq2 object
dds <- DESeq2::DESeqDataSetFromMatrix(bulk$RNA,
    colData = coldata,
    design = ~sample
)

coldata_cidp_vn_ctrl_ciap <-
    coldata |>
    dplyr::filter(level2 %in% c("CIDP", "VN", "CTRL", "CIAP"))

dds_cidp_vn_ctrl_ciap <- DESeq2::DESeqDataSetFromMatrix(bulk_cidp_vn_ctrl_ciap$RNA,
    colData = coldata_cidp_vn_ctrl_ciap,
    design = ~sample
)

# normalize
rld <- DESeq2::vst(dds, blind = TRUE)
rld_cidp_vn_ctrl_ciap <- DESeq2::vst(dds_cidp_vn_ctrl_ciap, blind = TRUE)

# save object
qs::qsave(rld, file.path("objects", "rld.qs"))

rld <- qs::qread(file.path("objects", "rld.qs"))

rld$level2 <- factor(rld$level2, sc_merge@misc$level2_order)


# plot PCA using DESeq2
DESeq2::plotPCA(rld, intgroup = "level2") +
    scale_color_manual(values = sc_merge@misc$level2_cols) + 
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        aspect.ratio = 1
    ) +
    geom_text(aes(label = coldata$sample), vjust = 2)

ggsave(file = file.path("results", "pca", "deseq2_pca_level2.pdf"), width = 10, height = 8)


pca_cidp_vn_ctrl_ciap <- DESeq2::plotPCA(rld_cidp_vn_ctrl_ciap, intgroup = "level2", ntop = 50, returnData = TRUE)

pca_cidp_vn_ctrl_ciap |>
    ggplot(aes(x = PC1, y = PC2, color = level2)) + 
    geom_point(size = 2, alpha = .5, shape = 16) + 
    scale_color_manual(values = sc_merge@misc$level2_cols) + 
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        aspect.ratio = 1
    ) +
    geom_text(aes(label = coldata_cidp_vn_ctrl_ciap$sample), vjust = 2)

ggsave(file = file.path("results", "pca", "deseq2_pca_level2.pdf"), width = 3, height = 4)
ggsave(file = file.path("results", "pca", "deseq2_pca_level2.pdf"), width = 3, height = 4)


library(factoextra)
data(decathlon2)
decathlon2_active <- decathlon2[1:12, 1:10]
res_pca <- prcomp(decathlon2_active, scale = TRUE)

# deseq2 pca manually ----
rv <- MatrixGenerics::rowVars(assay(rld_cidp_vn_ctrl_ciap))
select <- order(rv, decreasing = TRUE)[seq_len(50)]
pca_cidp_vn_ctrl_ciap <- prcomp(t(assay(rld_cidp_vn_ctrl_ciap)[select,]))

fviz_pca_var(pca_cidp_vn_ctrl_ciap, col.var = "contrib", select.var = list(contrib = 25), repel = TRUE)
ggsave(file = file.path("results", "pca", "deseq2_pca_level2_contrib.pdf"), width = 5, height = 5)

pca_cidp_vn_ctrl_ciap$x |>
    as_tibble() |>
    mutate(level2 = rld_cidp_vn_ctrl_ciap$level2) |>
    ggplot(aes(x = PC1, y = PC2, color = level2)) +
    geom_point(size = 2, alpha = .5, shape = 16) +
    scale_color_manual(values = sc_merge@misc$level2_cols) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        aspect.ratio = 1
    ) +
    geom_text(aes(label = coldata_cidp_vn_ctrl_ciap$sample), vjust = 2)

ggsave(file = file.path("results", "pca", "deseq2_pca_level2.pdf"), width = 3, height = 4)
ggsave(file = file.path("results", "pca", "deseq2_pca_level2_label.pdf"), width = 3, height = 4)

# plot PCA using DESeq2 color by normal_axon
DESeq2::plotPCA(rld, intgroup = "log_axon_normal") +
    viridis::scale_color_viridis(option = "magma") +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        aspect.ratio = 1
    ) +
    geom_text(aes(label = coldata$sample), vjust = 2) + 
    labs(color = "log normal axon")

ggsave(file = file.path("results", "pca", "deseq2_pca_axon_normal.pdf"), width = 10, height = 8)


# pseudobulk DGE but clusterwise approach ---
bulk_sample <- AverageExpression(
    sc_merge,
    method = "aggregate",
    return.seurat = FALSE,
    slot = "counts",
    assays = "RNA",
    group.by = "sample"
)

# pseudobulk MDS with limma ---
dge <- edgeR::DGEList(counts = bulk_sample$RNA, group = colnames(bulk_sample$RNA))
count_check <- edgeR::cpm(dge) > 1
keep <- which(rowSums(count_check) > 2)
dge <- dge[keep,]

str(bulk_sample)

str(sc_merge@meta.data)

# pseudobulk cells by stimulation condition AND cell type

# function to calculate pseuodbulk using Libra
dePseudo <- function(seu_obj, cell_type_col, label_col) {
    seu_obj_parse <- deparse(substitute(seu_obj))
    res <- Libra::run_de(seu_obj,
        replicate_col = "sample",
        cell_type_col = cell_type_col,
        label_col = label_col,
        min_cells = 3,
        min_reps = 2,
        min_feature = 0,
        de_family = "pseudobulk",
        de_method = "edgeR",
        de_type = "LRT",
        n_threads = 6
    )
    # res_count <-
    #     res |>
    #     dplyr::count(gene) |>
    #     arrange(desc(n)) |>
    #     dplyr::filter(n < 10)

    # res_sig <-
    #     res |>
    #     dplyr::filter(p_val < 0.05) |>
    #     # dplyr::filter(gene %in% res_count$gene) |>
    #     arrange(desc(avg_logFC))


    res <- arrange(res, cell_type, desc(avg_logFC))
    res_split <- split(res, res$cell_type)
    # res_split_sig <- split(res_sig, res_sig$cell_type)
    write_xlsx(res_split, path = file.path("results", "de", paste0(seu_obj_parse, ".xlsx")))
    # write_xlsx(res_split_sig, path = file.path("results", "de", paste0(seu_obj_parse, "_sig.xlsx")))
}

#needs to be dgCMatrix for Libra not BPCellsMatrix
# sc_merge$RNA$counts <- as(object = sc_merge[["RNA"]]$counts, Class = "dgCMatrix")

# PNP vs CTRL pseudobulk
dePseudo(sc_merge, cell_type_col = "cluster", label_col = "level0")


# function to calculate DEG for each cluster with Seurat wilcox
findMarkersWilcox <- function(cluster, condition1, condition2) {
    ident1 <- paste0(condition1, "_", cluster)
    ident2 <- paste0(condition2, "_", cluster)
    deg <-
        FindMarkers(sc_merge, ident.1 = ident1, ident.2 = ident2, only.pos = TRUE, logfc.threshold = 0.5) |>
        tibble::rownames_to_column("gene") |>
        dplyr::filter(p_val_adj < 0.05) |>
        dplyr::relocate(gene, avg_log2FC, p_val, p_val_adj) |>
        dplyr::arrange(desc(avg_log2FC))
    message("DEG found for ", condition1, " vs ", condition2, " in ", cluster)
    return(deg)
}

dplyr::count(sc_merge@meta.data, level0_cluster, sort = TRUE) |>
    write_xlsx(path = file.path("results", "de", "level0_cluster_count.xlsx"))

# PNP vs CTRL wilcox ----
level0_cluster_order <-
    expand_grid(
        level0 = c("CTRL", "PNP"),
        cluster = sc_merge@misc$cluster_order
    ) |>
    mutate(
        level0_cluster = paste0(level0, "_", cluster)
    )

sc_merge$level0_cluster <- 
    factor(paste0(sc_merge$level0, "_", sc_merge$cluster),
        levels = level0_cluster_order$level0_cluster
    )

Idents(sc_merge) <- sc_merge$level0_cluster

pnp_ctrl_deg_wilcox <-
    lapply(
        sc_merge@misc$cluster_order,
        function(cluster) {
            findMarkersWilcox(cluster, "PNP", "CTRL")
        }
    ) |>
    setNames(sc_merge@misc$cluster_order)

write_xlsx(pnp_ctrl_deg_wilcox, path = file.path("results", "de", "pnp_ctrl_wilcox.xlsx"))

sheets <- readxl::excel_sheets(path = file.path("results", "de", "pnp_ctrl_pseudobulk.xlsx"))

    pseudobulk_de <-
        lapply(
            sheets,
            function(sheet) {
                read_xlsx(path = file.path("results", "de", "pnp_ctrl_pseudobulk.xlsx"), sheet = sheet) |>
                    dplyr::filter(p_val_adj < 0.05)
            }
        ) |>
        setNames(sheets)

write_xlsx(pseudobulk_de, path = file.path("results", "de", "pnp_ctrl_pseudobulk_sig.xlsx"))

# plot number of DEG per cluster ---
plotDE <- function(name, title) {
    sheets <- readxl::excel_sheets(path = file.path("results", "de", paste0(name, ".xlsx")))
    cl_sig <-
        lapply(
            sheets,
            function(sheet) {
                read_xlsx(path = file.path("results", "de", paste0(name, ".xlsx")), sheet = sheet) |>
                    dplyr::filter(p_val_adj < 0.05) |>
                    nrow()
            }
        )
    result <- tibble(
        cluster = sheets,
        n = unlist(cl_sig)
    )
    plot <-
        result |>
        mutate(cluster = fct_reorder(cluster, n)) |>
        ggplot(aes(x = cluster, y = n, fill = cluster)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = sc_merge@misc$cluster_col) +
        theme_bw() +
        theme(legend.position = "none") +
        labs(
            x = "",
            y = "",
            title = title
        )
    ggsave(
        plot = plot,
        filename = file.path("results", "de", paste0(name, ".pdf")),
        width = 3,
        height = 3
    )
}

plotDE("pnp_ctrl_wilcox", title = "PNP vs CTRL")
plotDE("pnp_ctrl_pseudobulk", title = "PNP vs CTRL")


# VN vs CTRL
vn_ctrl <- subset(sc_merge, level2 %in% c("VN", "CTRL"))
vn_ctrl$level2 <- factor(vn_ctrl$level2, levels = c("VN", "CTRL"))

#sanity check
table(sc_merge$level2)
table(vn_ctrl$level2)

dePseudo(vn_ctrl, cell_type_col = "cluster", label_col = "level2")

# CIDP vs CTRL
cidp_ctrl <- subset(sc_merge, level2 %in% c("CIDP", "CTRL"))
cidp_ctrl$level2 <- factor(cidp_ctrl$level2, levels = c("CIDP", "CTRL"))

#sanity check
table(sc_merge$level2)
table(cidp_ctrl$level2)

dePseudo(cidp_ctrl, cell_type_col = "cluster", label_col = "level2")


# # using a linear model to adjust for sex and age is probably too conservative
# sample_lookup <- 
#   readr::read_csv(file.path("lookup", "sample_lookup.csv")) |>
#   janitor::clean_names() |>
#   mutate(age_calc = lubridate::time_length(difftime(nerve_date, birth_date), "years")) |>
#   mutate(age_calc = floor(age_calc)) |>
#   mutate(age = coalesce(age_calc, age)) |>
#   dplyr::select(-age_calc)

# sc_merge@misc$cluster_order
# runLimma(seurat_object = vn_ctrl, cluster = "mySC", condition1 = "VN", condition2 = "CTRL")

# runLimma <- function(seurat_object, cluster, condition1, condition2) {
#     pseudobulk_data <- Libra::to_pseudobulk(seurat_object,
#         cell_type_col = "cluster",
#         label_col = "level2",
#         min_cells = 10,
#         min_features = 3,
#         replicate_col = "sample"
#     )

#     dge <- edgeR::DGEList(counts = pseudobulk_data[[cluster]], group = colnames(pseudobulk_data[[cluster]]))
#     count_check <- edgeR::cpm(dge) > 1
#     keep <- which(rowSums(count_check) > 2)
#     dge <- dge[keep, ]
#     dge <- edgeR::calcNormFactors(dge, method = "TMM")
#     meta_limma <-
#         tibble(sample = str_extract(colnames(dge), pattern = "[^:]+")) |>
#         dplyr::left_join(sample_lookup, by = "sample") |>
#         dplyr::distinct(sample, .keep_all = TRUE)
#     designMat <- model.matrix(~ 0 + level2 + sex + age, data = meta_limma)
#     # designMat <- model.matrix(~ 0 + cohort + age + sex,  data = meta_limma)
#     my_contrasts <- glue::glue("level2{condition1}-level2{condition2}")
#     my_args <- list(my_contrasts, levels = designMat)
#     my_contrasts <- do.call(limma::makeContrasts, my_args)
#     dge_voom <- limma::voomWithQualityWeights(dge, designMat, plot = FALSE)

#     dge_voom <- dge_voom |>
#         limma::lmFit(design = designMat, block = NULL) |>
#         limma::contrasts.fit(my_contrasts) |>
#         limma::eBayes(robust = TRUE)

#     topgenes_sig <- limma::topTable(dge_voom, n = Inf, adjust.method = "BH") |>
#         dplyr::filter(adj.P.Val < 0.05) |>
#         tibble::rownames_to_column("gene") |>
#         tibble::tibble() |>
#         dplyr::arrange(desc(logFC)) |>
#         dplyr::rename(
#             avg_log2FC = logFC,
#             p_val_adj = adj.P.Val
#         )

#     readr::write_csv(topgenes_sig, file.path("results", "de", glue::glue("de_{condition1}_{condition2}_{cluster}_sig.csv")))

#     topgenes_all <- limma::topTable(dge_voom, n = Inf, adjust.method = "BH") |>
#         tibble::rownames_to_column("gene") |>
#         tibble::tibble() |>
#         dplyr::arrange(desc(logFC)) |>
#         dplyr::rename(
#             avg_log2FC = logFC,
#             p_val_adj = adj.P.Val
#         )

#     readr::write_csv(topgenes_all, file.path("results", "de", glue::glue("de_{condition1}_{condition2}_{cluster}_all.csv")))
# }
