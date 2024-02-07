# description: calculcate module/gene scores

# libraries ---
library(Seurat)
library(tidyverse)
library(qs)
library(BPCells)

my_cols_25 <- pals::cols25()

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"))
ic <- qs::qread(file.path("objects", "ic.qs"))

ic_cidp_vn_ctrl_ciap <- subset(ic, level2 %in% c("CIDP", "VN", "CTRL", "CIAP"))

markers <- read_csv(file.path("lookup", "markers.csv"))

str(ic_cidp_vn_ctrl_ciap@meta.data)

modules <- c(
    "szabo_proliferation",
    "szabo_cd8_cytotoxic",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_FATTY_ACID_METABOLISM",
    "HALLMARK_CHOLESTEROL_HOMEOSTASIS"
)

# calculate modules
for (module in modules) {
    ic_cidp_vn_ctrl_ciap <-
        AddModuleScore(
            object = ic_cidp_vn_ctrl_ciap,
            features = list(markers[[module]]),
            name = module,
            assay = "RNA",
        )
}

# better names for modules
for (module in modules) {
    module1 <- paste0(module, "1")
    ic_cidp_vn_ctrl_ciap@meta.data[[module]] <- ic_cidp_vn_ctrl_ciap@meta.data[[module1]]
    ic_cidp_vn_ctrl_ciap@meta.data[[module1]] <- NULL
}

# plot modules
module_plots <- lapply(modules,
    function(x) {
        scMisc::ModulePlot(
            object = ic_cidp_vn_ctrl_ciap,
            x_var = "level2",
            module = x,
            color = my_cols_25
        )
    })


module_plots_patch <- patchwork::wrap_plots(module_plots, ncol = 1)

ggsave(
    plot = module_plots_patch,
    file.path("results", "module", "ic_cidp_vn_ctrl_ciap_modules.pdf"), width = 5, height = 20
)

#UCell
renv::hydrate("UCell")
library(UCell)

ic_cidp_vn_ctrl_ciap

markers <- read_csv(file.path("lookup", "markers.csv"))

str(ic_cidp_vn_ctrl_ciap@meta.data)

modules <- c(
    "szabo_proliferation",
    "szabo_cd8_cytotoxic",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_FATTY_ACID_METABOLISM",
    "HALLMARK_CHOLESTEROL_HOMEOSTASIS"
)

markers_ucell <- list()

for (module in modules) {
    markers_ucell[[module]] <- markers[[module]]
    markers_ucell[[module]] <- markers_ucell[[module]][!is.na(markers_ucell[[module]])]
}

ic_cidp_vn_ctrl_ciap <- AddModuleScore_UCell(ic_cidp_vn_ctrl_ciap, features = markers_ucell)

names_ucell <- paste0(modules, "_UCell")

# VlnPlot(ic_cidp_vn_ctrl_ciap, features = names_ucell, group.by = "level2", ncol = 1)
# ggsave(file.path("results", "module", "ic_cidp_vn_ctrl_ciap_modules_UCell.pdf"), width = 5, height = 20)

# VlnPlot(ic_cidp_vn_ctrl_ciap, features = modules, group.by = "level2", ncol = 1)
# ggsave(file.path("results", "module", "ic_cidp_vn_ctrl_ciap_modules.pdf"), width = 7, height = 20)


# plot modules u cell
module_plots_ucell <- lapply(names_ucell,
    function(x) {
        scMisc::ModulePlot(
            object = ic_cidp_vn_ctrl_ciap,
            x_var = "level2",
            module = x,
            color = my_cols_25
        )
    })


module_plots_ucell_patch <- patchwork::wrap_plots(module_plots_ucell, ncol = 1)

ggsave(
    plot = module_plots_ucell_patch,
    file.path("results", "module", "ic_cidp_vn_ctrl_ciap_ucell_modules.pdf"), width = 5, height = 20
)

qs::qsave("ic_cidp_vn_ctrl_ciap", file.path("objects", "ic_cidp_vn_ctrl_ciap.qs"))
