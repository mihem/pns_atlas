#===============================================================================
# MrVI Model Analysis Script
#===============================================================================
# Purpose: Run the MrVI model on preprocessed single-cell data to obtain latent 
# representations and analyze sample distances.
#===============================================================================

# Load Libraries ----
import mrvi
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pickle
import pandas as pd

print("Libraries loaded")

# Load Preprocessed Data ----
adata = anndata.read_h5ad("objects/sc_diet.h5ad")
print("Data loaded")

# Setup the MrVI model with the appropriate keys ---
mrvi.MrVI.setup_anndata(adata, sample_key="sample", categorical_nuisance_keys=["center"])
mrvi_model = mrvi.MrVI(adata)
print("Preprocessing done")

# Train the MrVI model
mrvi_model.train()
print("Model trained")

# Extract Latent Representations ----
# Get latent representations from the trained MrVI model
adata.obsm["X_mrvi_z"] = mrvi_model.get_latent_representation(give_z=True)
adata.obsm["X_mrvi_u"] = mrvi_model.get_latent_representation(give_z=False)

# Save the updated AnnData object
adata.write("sc_diet_mrvi.h5ad")
print("adata saved")

# Extract Sample Representations and Distances ----
# Get local sample representations and distances
cell_sample_representations = mrvi_model.get_local_sample_representation()
cell_sample_sample_distances = mrvi_model.get_local_sample_representation(return_distances=True)

# Save the MrVI model and extracted data
mrvi_model.save("mrvi_model")
print("Model saved")

pickle.dump(cell_sample_representations, file=open("cell_sample_representations.pickle", "wb"))
pickle.dump(cell_sample_sample_distances, file=open("cell_sample_sample_distances.pickle", "wb"))
print("Sample representation and distances saved")

# Reload Data for Further Analysis ----
adata = anndata.read_h5ad("objects/sc_diet_mrvi.h5ad")

# Perform UMAP Visualization ----
import scanpy as sc

# Compute neighbors using the latent representation
sc.pp.neighbors(adata, use_rep="X_mrvi_u")

# Compute UMAP embedding
sc.tl.umap(adata, min_dist=0.1)

# Plot UMAP colored by level2
sc.pl.umap(
    adata,
    color=["level2"],
    frameon=False,
)

# Plot heatmap of marker genes
sc.pl.heatmap(
    adata,
    markers,
    groupby="cell_type",
    layer="scvi_normalized",
    standard_scale="var",
    dendrogram=True,
    figsize=(8, 12),
)

# Compute and Save Average Sample Distances ----
# Load the pickle object
with open("objects/cell_sample_sample_distances.pkl", "rb") as file:
    cell_sample_sample_distances = pickle.load(file)

# Compute the average over the first dimension
averaged_array = np.mean(cell_sample_sample_distances, axis=0)

# Create a DataFrame from the array
sample_order = adata.obs.loc[
    lambda x: ~x["sample"].duplicated(keep="first")
].sort_values("_scvi_sample")["sample"].values
averaged_df = pd.DataFrame(averaged_array, columns=sample_order, index=sample_order)

# Save the averaged distances as a CSV file
averaged_df.to_csv("results/mrvi/mrvi_average_all.csv")

# Compute and Save Cluster-Specific Average Distances ----
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
    df_mean_distances = pd.DataFrame(mean_distances, columns=sample_order, index=sample_order)
    
    # Add the DataFrame to the dictionary with the cluster name as the key
    cluster_dataframes[cluster] = df_mean_distances

# Save the cluster-specific average distances as CSV files
for cluster, df in cluster_dataframes.items():
    df.to_csv(f"results/mrvi/mrvi_cluster_{cluster}_average.csv")
