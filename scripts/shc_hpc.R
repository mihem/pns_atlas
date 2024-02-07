
message("Running scSHC")

library(scSHC)
message("Libraries loaded")

counts_small <- qs::qread("counts_small.qs")
message("Counts loaded")

meta_data <- qs::qread("sc_small_meta_data.qs")
message("Meta data loaded")

message("Running scSHC")
new_clusters <-
    scSHC::testClusters(
        counts_small,
        cluster_ids = as.character(meta_data$RNA_snn_res.0.4),
        batch = meta_data$sample,
        parallel = TRUE,
        cores = 6
    )
message("ScSHC done")


message("Save object")
qs::qsave(new_clusters, "sc_shc_small.qs")

message("Done")
