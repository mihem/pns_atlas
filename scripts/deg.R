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
    write_xlsx(res_split, path = file.path("results", "de", paste0(seu_obj_parse, ".xlsx")))
}

#needs to be dgCMatrix for Libra not BPCellsMatrix
# sc_merge$RNA$counts <- as(object = sc_merge[["RNA"]]$counts, Class = "dgCMatrix")

# PNP vs CTRL pseudobulk ----
dePseudo(sc_merge, cell_type_col = "cluster", label_col = "level0")

# VN vs CTRL pseudobulk ---
vn_ctrl <- subset(sc_merge, level2 %in% c("VN", "CTRL"))
vn_ctrl$level2 <- factor(vn_ctrl$level2, levels = c("VN", "CTRL"))


##sanity check
table(sc_merge$level2)
table(vn_ctrl$level2)
AggregateExpression(vn_ctrl, assay = "RNA", group.by = "level2", features = c("F2RL3", "CXCL14", "XIST", "TSIX"))

## perform DE
dePseudo(vn_ctrl, cell_type_col = "cluster", label_col = "level2")

# CIDP vs CTRL pseudoublk ----
cidp_ctrl <- subset(sc_merge, level2 %in% c("CIDP", "CTRL"))
cidp_ctrl$level2 <- factor(cidp_ctrl$level2, levels = c("CIDP", "CTRL"))

# #sanity check
table(sc_merge$level2)
table(cidp_ctrl$level2)

# perform DE
dePseudo(cidp_ctrl, cell_type_col = "cluster", label_col = "level2")


# CIAP vs CTRL pseudoublk ----
ciap_ctrl <- subset(sc_merge, level2 %in% c("CIAP", "CTRL"))
ciap_ctrl$level2 <- factor(ciap_ctrl$level2, levels = c("CIAP", "CTRL"))

# #sanity check
table(sc_merge$level2)
table(ciap_ctrl$level2)

# perform DE
dePseudo(ciap_ctrl, cell_type_col = "cluster", label_col = "level2")

deSig <- function(name) {
    sheets <- readxl::excel_sheets(path = file.path("results", "de", paste0(name, ".xlsx")))
    de <-
        lapply(
            sheets,
            function(sheet) {
                read_xlsx(path = file.path("results", "de", paste0(name, ".xlsx")), sheet = sheet) |>
                    dplyr::filter(p_val_adj < 0.1) |>
                    dplyr::filter(abs(avg_logFC) > 2)
            }
        )
    de <- setNames(de, sheets)
    write_xlsx(de, path = file.path("results", "de", paste0(name, "_sig.xlsx")))
}

deSig("pnp_ctrl_pseudobulk")
deSig("vn_ctrl")
deSig("cidp_ctrl")
deSig("ciap_ctrl")

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

# plotDE("pnp_ctrl_wilcox", title = "PNP vs CTRL")
plotDE("pnp_ctrl_pseudobulk", title = "PNP vs CTRL")

# volcano plot
volcanoPlot <- function(filename, sheet, FCcutoff = 2, selectLab = NULL, drawConnectors = TRUE, condition1, condition2) {
    input <- readxl::read_excel(file.path("results", "de", paste0(filename, ".xlsx")), sheet = sheet)
    if (nrow(input) != 0) {
        volcano <- EnhancedVolcano::EnhancedVolcano(
            data.frame(input),
            lab = input[["gene"]],
            x = "avg_logFC",
            y = "p_val_adj",
            pCutoff = 0.1,
            FCcutoff = FCcutoff,
            axisLabSize = 15,
            pointSize = 1,
            labSize = 2,
            subtitle = NULL,
            caption = NULL,
            border = "full",
            gridlines.major = FALSE,
            gridlines.minor = FALSE,
            drawConnectors = drawConnectors,
            lengthConnectors = unit(0.0001, "npc"),
               title = paste(condition1, "vs", condition2, "in ", sheet),
            boxedLabels = TRUE,
            selectLab = selectLab,
            xlab = "Log2 fold change",
            ylab = "-Log10 adjusted pvalue",
            legendLabels = c(
                "NS", "logFC",
                "p-val", "p-val + logFC"
            ),
            legendPosition = "right"
        )
        pdf(file.path("results", "de", paste0(filename, "_", sheet, ".pdf")), width = 8, height = 6)
        print(volcano)
        dev.off()
    }
}

# lab_blood <- list("actCD4" = c("NKG7", "GNLY", "GZMB", "KLDR1", "CCL5", "HLA-DRB1", "HLA-DRB5", "HLA-DRA", "LTB", "CCR7"), "naiveBc" = c("FOS", "JUNB", "FCER1G", "IFITM3"))

cluster_de <- c("repairSC", "mySC", "nmSC", "PC2", "ven_capEC1", "artEC")

# PNP vs CTRL
lapply(
    cluster_de,
    function(cluster) {
        volcanoPlot(
            filename = "pnp_ctrl_pseudobulk",
            sheet = cluster,
            FCcutoff = 2,
            condition1 = "PNP",
            condition2 = "CTRL"
        )
    }
)

# VN vs CTRL
lapply(
    cluster_de,
    function(cluster) {
        volcanoPlot(
            filename = "vn_ctrl",
            sheet = cluster,
            FCcutoff = 2,
            condition1 = "VN",
            condition2 = "CTRL"
        )
    }
)

# ciap vs CTRL
lapply(
    cluster_de,
    function(cluster) {
        volcanoPlot(
            filename = "ciap_ctrl",
            sheet = cluster,
            FCcutoff = 2,
            condition1 = "ciap",
            condition2 = "CTRL"
        )
    }
)
# CIAP vs CTRL
lapply(
    cluster_de,
    function(cluster) {
        volcanoPlot(
            filename = "ciap_ctrl",
            sheet = cluster,
            FCcutoff = 2,
            condition1 = "CIAP",
            condition2 = "CTRL"
        )
    }
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
