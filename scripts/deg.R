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
library(EnhancedVolcano)


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

#volcano plot
volcanoPlot <- function(filename, sheet, FCcutoff = 0.5, selectLab = NULL, drawConnectors = TRUE) {
    input <- readxl::read_excel(file.path("results", "de", paste0(filename, ".xlsx")), sheet = sheet)
    if(nrow(input) != 0) { 
      volcano <- EnhancedVolcano::EnhancedVolcano(
        data.frame(input),
        lab = input[["gene"]],
        x = 'avg_logFC',
        y = 'p_val_adj',
        pCutoff = 0.05,
        FCcutoff = FCcutoff,
        axisLabSize = 25,
        pointSize = 5,
        labSize = 5,
        subtitle = NULL,
        caption = NULL,
        drawConnectors = drawConnectors,
        lengthConnectors = unit(0.0001, "npc"),
        title = NULL,
        #    title = paste(x, input),
        boxedLabels = FALSE,
        selectLab = selectLab[[x]],
        #xlim=c(0, 2),
        #  ylim =c(0,50),
        xlab = "Log2 fold change",
        ylab = "-Log10 adjusted pvalue",
        legendLabels = c('NS', "avg logFC",
                         'p-value', "p-value and avg logFC"))
      #    legendPosition = "bottom")
      pdf(file.path("results", "de", paste0(filename, "_", x, ".pdf")), width = 8, height = 12)
      print(volcano)
      dev.off()
    }
}

lab_blood <- list("actCD4" = c("NKG7", "GNLY", "GZMB", "KLDR1", "CCL5", "HLA-DRB1", "HLA-DRB5", "HLA-DRA", "LTB", "CCR7"), "naiveBc" = c("FOS", "JUNB", "FCER1G", "IFITM3"))

debugonce(volcanoPlot)

volcanoPlot(
    filename = "pnp_ctrl_pseudobulk",
    sheet = "repairSC",
    # selectLab = lab_blood
)


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
