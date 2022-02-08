# %%
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
# %%

sample_one = scv.read("../test/G328E2L2_scRNAseq_G328E2L3_CITEseq.loom", cache=False)
# sample_one = sample_one.var_names_make_unique
# .... 
# sample_n = anndata.read_loom("sample_n.loom")
# %%
sample_obs = pd.read_csv("../test/cellID_obs.csv")
cell_clusters = pd.read_csv("../test/cell_clusters.csv")

#%%
sample_one.obs = sample_one.obs.rename(index= lambda x:x.split(":")[-1].replace("x","-1"))
sample_one.obs.head()
# %%
sample_one = sample_one[np.isin(sample_one.obs.index, sample_obs["x"])]
sample_one.obs.head()
# %%
# Now that we have our Velocity file filtered based upon our Seurat object, we can go ahead and add UMAP coordinates. We'll first upload them:
umap = pd.read_csv("../test/cell_embeddings.csv")
#%%
# With the coordinates, we will need to make sure we add them so they match the order of the Cell IDs in our anndata object. Our Cell IDs are rownames in the observation layer of our object, so we can view them by using the following:

sample_one.obs.index
# Let's cast our index as a data frame and change the column name
#%%
sample_one_index = pd.DataFrame(sample_one.obs.index)
sample_one_index = sample_one_index.rename(columns = {0:'CellID'})
# Let's also change the first column of our UMAP data frame to the same name:
#%%
umap = umap.rename(columns = {'Unnamed: 0':'CellID'})
cell_clusters = cell_clusters.rename(columns = {'Unnamed: 0':'CellID'})
# Now if we merge our index dataframe with our UMAP, the order will match our anndata object.

umap_ordered = sample_one_index.merge(umap, on = "CellID")
clusters_ordered = sample_one_index.merge(cell_clusters, on = "CellID")
# Since we're certain the orders are the same, we can remove the first column of the data frame and add the UMAP coordinates to our anndata object.
#%%
umap_ordered = umap_ordered.iloc[:,1:]
sample_one.obsm['X_umap'] = umap_ordered.values
# Clusters and their cluster colors can be added in the same fashion (and again, they must match the order of the Cell IDs.) Instead of adding them as an multidimensional observation ('obsm'), we'd add them under the unstructured annotation 'uns.'
clusters_ordered = clusters_ordered.iloc[:,1:]
# celltype_ordered = celltype_ordered.iloc[:,1:]

sample_one.uns['clusters'] = clusters_ordered.values
# sample_one.obs['celltype'] = celltype_ordered.values
# Running RNA Velocity
# At this point, we can now run the scVelo commands and generate our RNA Velocity plot based upon our Seurat UMAP coordinates.
#%%
scv.pp.filter_and_normalize(sample_one, min_shared_counts=30, n_top_genes=2000)
#%%
scv.pp.moments(sample_one, n_pcs=30, n_neighbors=30)
scv.tl.velocity(sample_one, mode = "stochastic")
#%%
scv.tl.velocity_graph(sample_one)
#%%
scv.pl.velocity_embedding(sample_one, basis = 'X_umap', arrow_size=5)
scv.pl.velocity_embedding_stream(sample_one,basis="umap",color="Clusters")
# If you want to incorporate your clusters and cluster colors in the embedding plot, the parameters you would add would be:
#%%
scv.pl.scatter(sample_one, color = "seurat_clusters")

# color = sample_one.uns['Cluster_colors']

#%%
sample_one
#%%