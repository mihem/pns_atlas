# ===============================================================================
# Analysis of perineurium
# ===============================================================================
# Purpose: Analyze manual measurements of perineurium thickness
# ===============================================================================

# load libraries ----
library(tidyverse)
library(readxl)

# define metadata ----
# use the same order and color palette as in the rest of the analysis
diagnosis_order <- c(
    "CTRL",
    "VN",
    "CIDP",
    "CIAP",
    "PPN",
    "DPN",
    "OIN",
    "ONIN",
    "MNC",
    "IN"
)
diagnosis_col <- setNames(
    pals::cols25(length(diagnosis_order)),
    diagnosis_order
)
# we only compare CTRL vs CIDP here
diagnosis_order <- c("CTRL", "CIDP")

# read perineurium lookup and measurement data ----
perineurium_lookup <- read_csv(file.path("lookup", "perineurium_lookup.csv"))

perineurium <-
    read_csv(file.path("lookup", "perineurium_measurement.csv")) |>
    dplyr::filter(is.na(remove) | !remove == "yes") |>
    left_join(perineurium_lookup) |>
    mutate(diagnosis = factor(diagnosis, levels = diagnosis_order)) |>
    mutate(sample_ordered = reorder(sample, as.numeric(diagnosis))) |>
    dplyr::mutate(
        outer_inner_diameter_ratio = outer_diameter / inner_diameter
    ) |>
    dplyr::mutate(outer_inner_area_ratio = outer_area / inner_area)

# sanity check
all(perineurium$inner_area < perineurium$outer_area)
all(perineurium$inner_diameter < perineurium$outer_diameter)
all(!duplicated(count(perineurium, sample, diagnosis)$sample))
perineurium[perineurium$inner_diameter > perineurium$outer_diameter, ]
print(dplyr::count(perineurium, sample, diagnosis), n = Inf)
select(
    arrange(perineurium, desc(outer_inner_area_ratio)),
    sample,
    outer_area,
    inner_area,
    outer_inner_area_ratio
)

# plot parameters per sample ----
plot_peri_sample <- function(param) {
    plot <- ggplot(
        perineurium,
        aes(x = sample_ordered, y = .data[[param]], fill = diagnosis)
    ) +
        geom_boxplot() +
        geom_point() +
        theme_classic() +
        scale_fill_manual(values = diagnosis_col)

    ggsave(
        file.path(
            "results",
            "perineurium",
            paste0("perineurium_", param, "_sample.pdf")
        ),
        plot = plot,
        width = 13,
        height = 5
    )
}

lapply(
    c(
        "outer_inner_diameter_ratio",
        "outer_inner_area_ratio"
    ),
    plot_peri_sample
)

# plot perineurial outer-inner area ratio per diagnosis ---
outer_inner_area_ratio_resid <- resid(lm(
    outer_inner_area_ratio ~ sex + age,
    data = perineurium
))

outer_inner_area_ratio_stats <- broom::tidy(wilcox.test(
    outer_inner_area_ratio_resid ~ perineurium$diagnosis
))


outer_inner_diameter_ratio_resid <- resid(lm(
    outer_inner_diameter_ratio ~ sex + age,
    data = perineurium
))

outer_inner_diameter_ratio_stats <- broom::tidy(wilcox.test(
    outer_inner_diameter_ratio_resid ~ perineurium$diagnosis
))


plot_peri_diagnosis <- function(param, y_lab) {
    # calculcate residuals in a first step to adjust for age and sex
    formula <- as.formula(paste0(param, " ~ sex + age"))
    resid_val <- resid(lm(
        formula,
        data = perineurium
    ))
    # perform Wilcoxon test on residuals
    stats <- broom::tidy(wilcox.test(
        resid_val ~ perineurium$diagnosis
    ))
    # plot boxplot
    plot <- ggplot(
        perineurium,
        aes(x = diagnosis, y = .data[[param]], fill = diagnosis)
    ) +
        geom_boxplot() +
        geom_jitter(size = 0.3, width = 0.3) +
        theme_classic() +
        ggsignif::geom_signif(
            comparisons = list(c("CTRL", "CIDP")),
            annotation = signif(stats$p.value, 3)
        ) +
        scale_fill_manual(values = diagnosis_col) +
        theme(legend.position = "none") +
        xlab("") +
        ylab(paste0("Perineurium outer-inner ratio ", y_lab))

    # save plot
    ggsave(
        file.path(
            "results",
            "perineurium",
            paste0("perineurium_", param, "_diagnosis.pdf")
        ),
        plot = plot,
        width = 2,
        height = 3.5
    )
}


# Define the parameters and labels
params <- c("outer_inner_diameter_ratio", "outer_inner_area_ratio")
labels <- c("diameter", "area")
map2(params, labels, plot_peri_diagnosis)
