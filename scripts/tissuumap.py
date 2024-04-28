# export Xenium to .h5 for tissuumap

# import libraries
import numpy as np
import pandas as pd

import os
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq

import anndata as ad
import pyarrow.parquet as pq

# adata list

# function to read in Xenium and convert to anndata
def create_adata(i):
    directory = xenium_names["File"][i]
    name = xenium_names["name"][i]
    adata = sc.read_10x_h5(filename=f"xenium/raw/{directory}/cell_feature_matrix.h5")
    df = pd.read_csv(f"xenium/raw/{directory}/cells.csv.gz")
    predictions = pq.read_table(f"results/xenium/xenium_predictions_{name}.parquet")
    df.set_index(adata.obs_names, inplace=True)
    adata.obs = df.copy()
    adata.obs["sample"] = xenium_names["sample"][i]
    adata.obs["condition"] = xenium_names["level2"][i]
    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()
    adata.obs["sn_predictions"] = predictions["sn_predictions"]
    adata.obs["sn_predictions_group"] = predictions["sn_predictions_group"]
    sc.pp.filter_cells(adata, min_counts = 10)
    
    return adata

# List all directories in the /xenium/raw path
directories = sorted(os.listdir("xenium/raw"))
directories_df = pd.DataFrame({"File": directories})

# read xenium meta data
xenium_meta = pd.read_excel("lookup/xenium_meta.xlsx")
xenium_meta["File"] = xenium_meta["File"].str.replace(".tar", "")

sample_lookup = pd.read_csv("lookup/sample_lookup.csv")

xenium_names = pd.merge(directories_df, xenium_meta, on = "File", how='inner')
xenium_names["sample"] = xenium_names["Name"].str.extract(r"(S\d+)")
xenium_names = pd.merge(xenium_names, sample_lookup, on = "sample", how='inner')
xenium_names["name"] = xenium_names["sample"] + "_" + xenium_names["level2"]

# Process data for each directory using the function
adata_list = [create_adata(i) for i in range(len(xenium_names))]

#  drop unwanted columns and save for tissuumap
def adata_write(i):
    adata_name = xenium_names["name"][i]
    adata = adata_list[i]
    adata.obs.drop("sample", axis=1, inplace=True)
    adata.obs.drop("condition", axis=1, inplace=True)
    adata.obs.drop("cell_id", axis=1, inplace=True)
    adata.obs.drop("x_centroid", axis=1, inplace=True)
    adata.obs.drop("y_centroid", axis=1, inplace=True)
    adata.obs.drop("cell_area", axis=1, inplace=True)
    adata.obs.drop("control_codeword_counts", axis=1, inplace=True)
    adata.obs.drop("control_probe_counts", axis=1, inplace=True)
    adata.obs.drop("deprecated_codeword_counts", axis=1, inplace=True)
    adata.obs.drop("nucleus_area", axis=1, inplace=True)
    adata.obs.drop("total_counts", axis=1, inplace=True)
    adata.obs.drop("transcript_counts", axis=1, inplace=True)
    adata.obs.drop("unassigned_codeword_counts", axis=1, inplace=True)
    adata.obs.drop("n_counts", axis=1, inplace=True)

    adata.write("objects/tissuumap_" + adata_name + ".h5ad")

for i in range(len(adata_list)):
    adata_write(i)