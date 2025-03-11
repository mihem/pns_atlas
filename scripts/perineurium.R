# ===============================================================================
# Analysis of perineurium
# ===============================================================================
# Purpose: Analyze manual measurements of perineurium thickness
# ===============================================================================

# load libraries ----
library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)

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


# Count the number of samples for each diagnosis
sample_counts <- perineurium |>
    group_by(diagnosis) |>
    summarise(n = n_distinct(sample))

print(sample_counts)

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
perineurium_median <- perineurium |>
    group_by(sample) |>
    mutate(
        outer_inner_area_ratio = median(outer_inner_area_ratio),
        outer_inner_diameter_ratio = median(outer_inner_diameter_ratio)
    ) |>
    ungroup() |>
    distinct(sample, .keep_all = TRUE)

plot_peri_diagnosis_mixed <- function(param, y_lab) {
    # Fit a linear mixed model
    formula <- as.formula(paste0(
        param,
        " ~ diagnosis + sex + age + (1|sample) + (1|center_sample) + (1|center_stain)"
    ))

    mixed_model <- lmer(formula, data = perineurium)

    # Get model summary and p-value for diagnosis effect
    model_summary <- summary(mixed_model)
    diagnosis_p_value <- coef(model_summary)[2, "Pr(>|t|)"]

    # Get fixed effects to calculate adjusted values
    fixed_effects <- fixef(mixed_model)

    # First adjust all individual measurements
    adjusted_perineurium <- perineurium %>%
        mutate(
            # Remove effects of all covariates except diagnosis
            covariate_effect = (sex == "male") *
                fixed_effects["sexmale"] +
                age * fixed_effects["age"],

            # Calculate adjusted value (raw value minus covariate effects)
            adjusted_value = .data[[param]] - covariate_effect
        )

    # Then calculate median of adjusted values per sample
    sample_data <- adjusted_perineurium |>
        group_by(sample, diagnosis) |>
        summarise(
            adjusted_value = median(adjusted_value),
            .groups = "drop"
        )

    # Plot
    plot <- ggplot(
        sample_data,
        aes(x = diagnosis, y = adjusted_value, fill = diagnosis)
    ) +
        geom_boxplot() +
        geom_jitter(width = 0.2, size = 1) +
        theme_classic() +
        ggsignif::geom_signif(
            comparisons = list(c("CTRL", "CIDP")),
            annotation = paste("p =", signif(diagnosis_p_value, 3))
        ) +
        scale_fill_manual(values = diagnosis_col) +
        theme(legend.position = "none") +
        xlab("") +
        ylab(paste0("Adjusted perineurium outer-inner ratio (", y_lab, ")"))

    # Save plot
    ggsave(
        file.path(
            "results",
            "perineurium",
            paste0("perineurium_mixed_model_", param, "_diagnosis.pdf")
        ),
        plot = plot,
        width = 2,
        height = 3.7
    )
}

params <- c("outer_inner_diameter_ratio", "outer_inner_area_ratio")
labels <- c("diameter", "area")
map2(params, labels, plot_peri_diagnosis_mixed)
