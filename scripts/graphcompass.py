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

import graphcompass as gc

# adata list

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
    
    return adata

# List all directories in the /xenium/raw path
directories = sorted(os.listdir("xenium/raw"))
directories_df = pd.DataFrame({"File": directories})

# read xenium meta data
xenium_meta = pd.read_excel("lookup/xenium_meta.xlsx")
xenium_meta["File"] = xenium_meta["File"].str.replace(".tar", "")


sample_lookup = pd.read_csv("lookup/sample_lookup.csv")

xenium_names = pd.merge(directories_df, xenium_meta, on = "File", how='inner')
xenium_names["sample"] =  xenium_names["Name"].str.extract(r"(S\d+)")
xenium_names = pd.merge(xenium_names, sample_lookup, on = "sample", how='inner')
xenium_names["name"] = xenium_names["sample"] + "_" + xenium_names["level2"]

# Process data for each directory using the function
adata_list = [create_adata(i) for i in range(len(xenium_names))]

# plot predictions
sq.pl.spatial_scatter(
    adata,
    library_key = "sample",
    library_id="S01",
    shape=None,
    color=[
        "sn_predictions_group",
    ],
    wspace=0.4,
)

plt.show()

adata = ad.concat(adata_list)
adata.obs["sample"] = adata.obs["sample"].astype("category")
adata.obs["condition"] = adata.obs["condition"].astype("category")
adata.obs["sn_predictions_group"] = adata.obs["sn_predictions_group"].astype("category")
adata.obs["sn_predictions"] = adata.obs["sn_predictions"].astype("category")

adata.obsm["spatial"]

adata.obs

adata.obs["sample"].value_counts()

# sq.gr.spatial_autocorr(
#     adata,
#     mode="moran",
#     n_perms=100,
#     n_jobs=1,
# )

# adata.uns["moranI"].head(10)
# adata.uns["moranI"].tail(10)

# sq.pl.spatial_scatter(
#     adata,
#     library_id="spatial",
#     color=[
#         "ACTA2",
#         "RGS5",
#     ],
#     shape=None,
#     size=2,
#     img=False,
# )

# plt.show()

# graphcompass
library_key="sample"
cluster_key="sn_predictions_group"
condition_key = "condition"

gc.tl.wlkernel.compare_conditions(
   adata=adata,
   library_key=library_key,
   cluster_key=cluster_key,
   compute_spatial_graphs=True,
   kwargs_spatial_neighbors={
        'coord_type': 'generic',
        'delaunay': True,  
  }  
)

adata.write("objects/graphcompass_adata.h5ad")
adata = ad.read_h5ad("objects/graphcompass_adata.h5ad")

adata1 = ad.read_h5ad("objects/graphcompass_adata.h5ad")

adata1.write("objects/graphcompass_adata1.h5ad")
adata.uns["filtration_curves"]
adata1.obs[["cell_id"]]

adata.uns["wl_kernel"]
adata.uns

help(gc.pl.wlkernel.compare_conditions)

# define necessary params
control_group="CTRL" # reference group
metric_key="wasserstein_distance" 
method="wl_kernel"

gc.pl.wlkernel.compare_conditions(
    adata=adata,
    library_key=library_key,
    condition_key=condition_key,
    control_group=control_group,
    metric_key=metric_key,
    method=method,
    figsize=(3,5),
    dpi=300,
    save="results/xenium/graphcompass/wwlkerenl.pdf"
)

# filtration curves
gc.tl.filtration_curves.compare_conditions(
    adata=adata,
    library_key=library_key,
    cluster_key=cluster_key,
    condition_key=condition_key,
    compute_spatial_graphs=False
    )
# Error in igraph._igraph.InternalError: Error at src/graph/type_indexededgelist.c:1414: Cannot get edge ID, no such edge. -- Invalid value
    
# define necessary params
node_labels=["EC", "IC", "PC", "SC", "VSMC", "endoC", "periC", "epiC"] # node labels (e.g. cell types) we are intrested in visualising
metric_key="filtration_curves"

gc.pl.filtration_curves.compare_conditions(
    adata=adata,
    node_labels=node_labels,
    metric_key=metric_key,
    palette="Set2",
    dpi=100,
    figsize=(30,5),
    save="results/xenium/graphcompass/filtration_curves.pdf"
)

# specifc cell type subgraph comparisons
gc.tl.distance.compare_conditions(
    adata=adata,
    library_key=library_key,
    cluster_key=cluster_key,
    method="portrait",
    compute_spatial_graphs=False,
)


import inspect
source_code = inspect.getsource(gc.tl.distance.compare_conditions)
print(source_code)

source_code = inspect.getsource(gc.tl_calculate_graph_distances)


adata.uns["pairwise_similarities"].to_csv("similarity_scores.csv", index=False)
adata.uns["wl_kernel"]["wasserstein_distance"].to_csv("wasserstein_distance.csv", index=True)

data_dict = adata.uns["wl_kernel"]
df = pd.DataFrame(list(data_dict.items()), columns=['Key', 'Value'])
df.to_csv("wl_kernel.csv", index=False)

gc.pl.distance.compare_conditions(
    adata=adata,
    library_key=library_key,
    condition_key=condition_key,
    control_group=control_group,
    add_ncells_and_density_plots=True,
    palette="Reds",
    dpi=300,
    figsize=(8,3),
    save="results/xenium/graphcompass/portrait.pdf"
)

import inspect
source_code = inspect.getsource(gc.pl.distance.compare_conditions)
print(source_code)

adata.obs["condition"].values.unique()


# # remove S22 because technical artifacts
# adata_subset = adata[~adata.obs["sample"].isin(["S22"])]

# gc.tl.distance.compare_conditions(
#     adata=adata_subset,
#     library_key=library_key,
#     cluster_key=cluster_key,
#     method="portrait",
#     compute_spatial_graphs=False,
# )

# gc.pl.distance.compare_conditions(
#     adata=adata_subset,
#     library_key=library_key,
#     condition_key=condition_key,
#     control_group=control_group,
#     add_ncells_and_density_plots=True,
#     palette="Reds",
#     dpi=300,
#     figsize=(8,3),
#     save="results/xenium/graphcompass/portrait_without_s22.pdf"
# )

adata = ad.read_h5ad("objects/graphcompass_adata_s01_s04_s24_s30.h5ad")


from graphcompass.tl._filtration_curves import _compute_edge_weights

# assumes you have computed the spatial graph before
adj_matrix = _compute_edge_weights(gene_expression_matrix=adata.X, adjacency_matrix=adata.obsp['spatial_connectivities'])
weights = adj_matrix.tocoo()

# check if there are any nans
has_nans = np.isnan(weights.data).any()
print(has_nans)

# check if there are any 0
has_zeros = np.isclose(weights.data, 0).any()
print(has_zeros)

## testing
breast = ad.read_h5ad("objects/mibitof_breast_cancer.h5ad")

sq.pl.spatial_scatter(breast, color = "phenotype", library_key = "Point_Num", library_id = "2203", shape = None)
sq.pl.spatial_scatter(breast, color = "phenotype", library_id = "Point_Num", shape = None)
sq.pl.spatial_scatter(breast, color = "phenotype", library_id = "sample", shape = None)

plt.show()

breast.obs
# sq.pl.spatial_scatter(breast, color = "cell_type_original", library_id = "sample", shape = None)
