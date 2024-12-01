# ===============================================================================
# Differential Expression Analysis (DEG) Script
# ===============================================================================
# Purpose: Analyze differential gene expression between disease conditions:
# - PNP vs CTRL
# - VN vs CTRL
# - CIDP vs CTRL
# - CIAP vs CTRL
#
# Methods: Uses pseudobulk differential expression via Libra package with edgeR
# ===============================================================================

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

# function to calculate pseuodbulk differential expression analysis using Libra
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
    res <- arrange(res, cell_type, desc(avg_logFC))
    res_split <- split(res, res$cell_type)
    write_xlsx(res_split, path = file.path("results", "de", paste0(seu_obj_parse, ".xlsx")))
}

# PNP vs CTRL pseudobulk ----
sc_merge$level0 <- factor(sc_merge$level0, levels = c("PNP", "CTRL"))
dePseudo(sc_merge, cell_type_col = "cluster", label_col = "level0")

nmSC <- subset(sc_merge, cluster %in% c("nmSC"))
AggregateExpression(nmSC, assay = "RNA", group.by = "level0", features = c("IFI44L"))

# VN vs CTRL pseudobulk ---
vn_ctrl <- subset(sc_merge, level2 %in% c("VN", "CTRL"))
vn_ctrl$level2 <- factor(vn_ctrl$level2, levels = c("VN", "CTRL"))

## sanity check
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

# deSig() filters significant DEGs based on thresholds:
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
        theme_classic() +
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
            lab = paste0("italic('", input[["gene"]], "')"),
            x = "avg_logFC",
            y = "p_val_adj",
            xlim = c(min(input[["avg_logFC"]], max(input[["avg_logFC"]]))),
            ylim = c(0, max(-log10(input[["p_val_adj"]]))),
            pCutoff = 0.1,
            FCcutoff = FCcutoff,
            axisLabSize = 15,
            pointSize = 2,
            labSize = 5,
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
            xlab = bquote(~ Log[2] ~ "fold change"),
            ylab = bquote(~ -Log[10] ~ "adjusted p-value"),
            parseLabels = TRUE,
            legendLabels = c(
                "NS", "logFC",
                "p-val", "p-val + logFC"
            ),
            legendPosition = "right",
        )
        pdf(file.path("results", "de", paste0(filename, "_", sheet, ".pdf")), width = 8, height = 6)
        print(volcano)
        dev.off()
    }
}

cluster_de <- c("repairSC", "mySC", "nmSC", "PC2")

lab_pnp_ctrl <- list(
    "mySC" = paste0("italic('", c("DCN", "TNXB", "COL1A1", "COL15A1", "CD53", "IL4R", "CD74"), "')"),
    "nmSC" = paste0("italic('", c("IL10RA", "IL13RA1", "CSF2RA", "TGFBI"), "')"),
    "repairSC" = paste0("italic('", c("GALR1", "TMEM47"), "')"),
    "PC2" = paste0("italic('", c("MFAP5", "NLGN4Y", "PCDH11Y", "IFIT3", "OASL", "MX1"), "')")
)

paste0("italic('", c("DCN", "TNXB"), "')")

# PNP vs CTRL
lapply(
    cluster_de,
    function(cluster) {
        volcanoPlot(
            filename = "pnp_ctrl_pseudobulk",
            sheet = cluster,
            FCcutoff = 2,
            condition1 = "PNP",
            condition2 = "CTRL",
            selectLab = lab_pnp_ctrl[[cluster]]
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

# function to compare gene expression between conditions using VlnPlot -----
compareGeneExpression <- function(seu_obj, gene, seu_obj_name) {
    plot <- VlnPlot(
        seu_obj,
        features = gene,
        group.by = "level2",
        cols = seu_obj@misc$level2_cols,
        pt.size = 0
    ) +
        NoLegend() + 
        xlab("") + 
        ylab("") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(
        plot = plot,
        filename = file.path("results", "de", paste0(seu_obj_name, "_", gene, ".pdf")),
        width = 3,
        height = 3
    )
}

# compare CXCL14 expression in pericytes between conditions
pericytes <- subset(
    sc_merge,
    cluster %in% c("periC1", "periC2", "periC3") &
        level2 %in% c("CTRL", "CIDP", "VN", "CIAP")
)
pericytes$level2 <- factor(pericytes$level2, levels = c("CTRL", "CIDP", "VN", "CIAP"))
compareGeneExpression(pericytes, "CXCL14", "pericytes")

# compare GRIK3 and PRIMA1 expression in nmSC between conditions
nmSC <- subset(
    sc_merge,
    cluster %in% c("nmSC") &
        level2 %in% c("CTRL", "CIDP", "VN", "CIAP")
)
nmSC$level2 <- factor(nmSC$level2, levels = c("CTRL", "CIDP", "VN", "CIAP"))

lapply(
    c("GRIK3", "PRIMA1"),
    function(gene) compareGeneExpression(nmSC, gene, "nmSC")
)

# compare MLIP expression in mySC between conditions
mySC <- subset(
    sc_merge,
    cluster %in% c("mySC") &
        level2 %in% c("CTRL", "CIDP", "VN", "CIAP")
)
mySC$level2 <- factor(mySC$level2, levels = c("CTRL", "CIDP", "VN", "CIAP"))
compareGeneExpression(mySC, "MLIP", "mySC")

