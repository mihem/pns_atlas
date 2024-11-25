# ===============================================================================
# Export Xenium Data to .h5 for Tissuumap
# ===============================================================================
# Purpose: Convert Xenium data to AnnData format and save it for visualization
# in Tissuumap.
# ===============================================================================

# Import necessary libraries
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
import anndata as ad
import pyarrow.parquet as pq

# Function to read in Xenium data and convert to AnnData


def create_adata(i):
    directory = xenium_names["File"][i]
    name = xenium_names["name"][i]
    # Read the 10x HDF5 file
    adata = sc.read_10x_h5(
        filename=f"xenium/raw/{directory}/cell_feature_matrix.h5")
    # Read the cell metadata
    df = pd.read_csv(f"xenium/raw/{directory}/cells.csv.gz")
    # Read the predictions
    predictions = pq.read_table(
        f"results/xenium/xenium_predictions_{name}.parquet")
    df.set_index(adata.obs_names, inplace=True)
    adata.obs = df.copy()
    # Filter cells with minimum counts
    sc.pp.filter_cells(adata, min_counts=10)
    # Add sample and condition metadata
    adata.obs["sample"] = xenium_names["sample"][i]
    adata.obs["condition"] = xenium_names["level2"][i]
    # Add spatial coordinates
    adata.obsm["spatial"] = adata.obs[[
        "x_centroid", "y_centroid"]].copy().to_numpy()
    # Add predictions to the AnnData object
    adata.obs["sn_predictions"] = predictions["sn_predictions"]
    adata.obs["sn_predictions_group"] = predictions["sn_predictions_group"]

    return adata


# List all directories in the /xenium/raw path
directories = sorted(os.listdir("xenium/raw"))
directories_df = pd.DataFrame({"File": directories})

# Read Xenium metadata
xenium_meta = pd.read_excel("lookup/xenium_meta.xlsx")
xenium_meta["File"] = xenium_meta["File"].str.replace(".tar", "")

# Read sample lookup table
sample_lookup = pd.read_csv("lookup/sample_lookup.csv")

# Merge directories with metadata
xenium_names = pd.merge(directories_df, xenium_meta, on="File", how='inner')
xenium_names["sample"] = xenium_names["Name"].str.extract(r"(S\d+)")
xenium_names = pd.merge(xenium_names, sample_lookup, on="sample", how='inner')
xenium_names["name"] = xenium_names["sample"] + "_" + xenium_names["level2"]

# Process data for each directory using the function
adata_list = [create_adata(i) for i in range(len(xenium_names))]

# Function to drop unwanted columns and save for Tissuumap
def adata_write(i):
    adata_name = xenium_names["name"][i]
    adata = adata_list[i]
    # Drop unnecessary columns
    adata.obs.drop(["sample", "condition", "cell_id", "x_centroid", "y_centroid",
                    "cell_area", "control_codeword_counts", "control_probe_counts",
                    "deprecated_codeword_counts", "nucleus_area", "total_counts",
                    "transcript_counts", "unassigned_codeword_counts", "n_counts"], axis=1, inplace=True)
    # Save the AnnData object
    adata.write("objects/tissuumap_" + adata_name + ".h5ad")


# Save each processed AnnData object
for i in range(len(adata_list)):
    adata_write(i)
