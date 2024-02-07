import mrvi
import anndata

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pickle
import pandas as pd

print("Libraries loaded")

adata = anndata.read_h5ad("objects/sc_diet.h5ad")
print("Data loaded")

type(adata.obs[["sex"]])

mrvi.MrVI.setup_anndata(adata, sample_key="sample", categorical_nuisance_keys=["center"])
mrvi_model = mrvi.MrVI(adata)
print("Preprocessing done")

mrvi_model.train()
print("Model trained")

adata.obsm["X_mrvi_z"] = mrvi_model.get_latent_representation(give_z = True)
adata.obsm["X_mrvi_u"] = mrvi_model.get_latent_representation(give_z = False)

adata.write("sc_diet_mrvi.h5ad")
print("adata saved")

cell_sample_representations = mrvi_model.get_local_sample_representation()
cell_sample_sample_distances = mrvi_model.get_local_sample_representation(return_distances=True)

mrvi_model.save("mrvi_model")
print("Model saved")

pickle.dump(cell_sample_representations, file = open("cell_sample_representations.pickle", "wb"))
pickle.dump(cell_sample_sample_distances, file = open("cell_sample_sample_distances.pickle", "wb"))
print("sample representation and distances saved")

# anaylsis after hpc run
adata = anndata.read_h5ad("objects/sc_diet_mrvi.h5ad")

adata.obsm["X_mrvi_u"]

import scanpy as sc

sc.pp.neighbors(adata, use_rep = "X_mrvi_u")

adata.obsp["distances"]

sc.tl.umap(adata, min_dist = 0.1)

sc.pl.umap(
    adata,
    color=["level2"],
    frameon=False,
)

sc.pl.heatmap(
    adata,
    markers,
    groupby="cell_type",
    layer="scvi_normalized",
    standard_scale="var",
    dendrogram=True,
    figsize=(8, 12),
)

adata = anndata.read_h5ad("objects/sc_diet_mrvi.h5ad")
adata

# Load the pickle object
with open("objects/cell_sample_sample_distances.pkl", "rb") as file:
    cell_sample_sample_distances = pickle.load(file)

cell_sample_sample_distances.shape
type(cell_sample_sample_distances)

sample_order = adata.obs.loc[
        lambda x: ~x["sample"].duplicated(keep="first")
    ].sort_values("_scvi_sample")["sample"].values

# Read in the CSV file from the lookup folder
sample_lookup = pd.read_csv("lookup/sample_lookup.csv")

# Compute the average over the first dimension
averaged_array = np.mean(cell_sample_sample_distances, axis=0)

# Create a DataFrame from the array
averaged_df = pd.DataFrame(averaged_array, columns = sample_order, index = sample_order)

#write csv of averaged_df
averaged_df.to_csv("results/mrvi/mrvi_average_all.csv")

# average for each cluster separately ---
# Get the unique clusters from adata.obs["cluster"]
unique_clusters = np.unique(adata.obs["cluster"])

# Create a dictionary to hold the separate DataFrames for each cluster
cluster_dataframes = {}

# Iterate over the clusters
for cluster in unique_clusters:
    # Filter the distances for the current cluster
    cluster_distances = cell_sample_sample_distances[adata.obs["cluster"] == cluster]
    
    # Compute the mean for the current cluster
    mean_distances = np.mean(cluster_distances, axis=0)
    
    # Create a DataFrame for the mean distances
    df_mean_distances = pd.DataFrame(mean_distances, columns = sample_order, index = sample_order)
    
    # Add the DataFrame to the dictionary with the cluster name as the key
    cluster_dataframes[cluster] = df_mean_distances

# save the cluster_dataframes as csv
for cluster, df in cluster_dataframes.items():
    df.to_csv(f"results/mrvi/mrvi_cluster_{cluster}_average.csv")





# # toy_df = pd.DataFrame(reshaped_array, index=row_labels, columns=col_labels)

# lut = dict(zip(species.unique(), "rbg"))
# row_colors = species.map(lut)

# # map the colors
# row_colors = sample_lookup.level1
# sample_lookup["sample"].map(sample_lookup.level1)

# sample_lookup.level2.unique()

# lut = dict(zip(sample_lookup.level2.unique(), "rgbrgbrgbrgb"))
# row_color = sample_lookup.level2.map(lut)

# import itertools

# # Create a dictionary to map each unique item in sample_lookup.level2 to a color
# color_mapping = {level: color for level, color in zip(sample_lookup.center.unique(), ['r', 'g', 'b'])}
# color_mapping = dict(zip(sample_lookup.center.unique(), "rbg"))
# color_mapping = {level: color for level, color in zip(sample_lookup.level2.unique(), itertools.cycle(['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w']))}

# # Map the colors to the rows in the clustermap based on sample_lookup.level2
# row_colors = sample_lookup.level2.map(color_mapping)
# row_colors = sample_lookup.center.map(color_mapping)

# #change columnames of row_colors
# row_colors.index = sample_order

# # Plot the clustermap with colored labels
# sns.clustermap(toy_df, cmap = "viridis", row_colors=row_colors)
# sns.clustermap(toy_df, row_colors=row_colors)
# plt.show()
# plt.savefig("output1.png")


# iris = sns.load_dataset("iris")
# species = iris.pop("species")
# sns.clustermap(iris)
# plt.savefig("output.png")

# lut = dict(zip(species.unique(), "rbg"))
# row_colors = species.map(lut)
# sns.clustermap(iris, row_colors=row_colors, cmap = "viridis")
# plt.savefig("output.png")

# iris.shape
# row_colors.shape

# # Create a toy array
# toy_array = np.array([[[1, 2], [3, 4]],
#                       [[5, 6], [7, 8]],
#                       [[9, 10], [11, 12]],
#                       [[13, 14], [15, 16]]])

# toy_array.shape
# toy_array

# # Reshape the array to 2-D
# reshaped_array = np.reshape(toy_array, (4, -1))
# reshaped_array.shape

# # Create row and column labels
# row_labels = [f"Cell {i+1}" for i in range(4)]
# col_labels = ["Sample A", "Sample B", "Sample A", "Sample B"]

# toy_df = pd.DataFrame(reshaped_array, index=row_labels, columns=col_labels)
# sns.clustermap(toy_df)
# plt.savefig("toy1.png")