import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import anndata
import scipy.sparse as sp
import numpy as np

sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
parser.add_argument('annotations')
args = parser.parse_args()

myObject =  args.myObject
annot_file = args.annotations

newObject = "annotated.h5ad"

sample = "WT"
combined_adata = sc.read_h5ad(myObject)

cluster_to_celltype_dict = {}
with open(annot_file, "r") as f:
    for line in f:
        cluster, celltype = line.strip().split(',')
        cluster_to_celltype_dict[cluster] = celltype

cluster_to_celltype_dict = {str(key): value for key, value in cluster_to_celltype_dict.items()}

# Remove unwanted clusters
clusters_to_remove = ["47","33","50","55","22","34","20","10","56"]
combined_adata = combined_adata[~combined_adata.obs["leiden"].isin(clusters_to_remove)].copy()

combined_adata.obs["celltype"] = combined_adata.obs["leiden"].map(cluster_to_celltype_dict)

figure_name = "figures/"+sample +"_annotationsON.png"
combined_adata.obs_names_make_unique()
fig = sc.pl.umap(combined_adata, color='celltype', legend_loc="on data", show=False, return_fig=True)
fig.savefig(figure_name, dpi=600, bbox_inches='tight')
plt.close(fig)

figure_name = "figures/"+sample +"_annotations.png"
fig = sc.pl.umap(combined_adata, color='celltype', show=False, return_fig=True)
fig.savefig(figure_name, dpi=600, bbox_inches='tight')
plt.close(fig)

combined_adata.write(newObject, compression="gzip")

