#!/usr/bin/env python3
import scanpy as sc
import sys
import importlib_metadata
import argparse
import numpy as np

# Fix importlib issue
sys.modules['importlib.metadata'] = importlib_metadata

# Parse input arguments
parser = argparse.ArgumentParser()
parser.add_argument('myObject')

args = parser.parse_args()

myObject = args.myObject
newObject = "analysed.h5ad" 
basename = "analysed" 

# Load data
combined_adata = sc.read(myObject)

# Normalize and log-transform
sc.pp.normalize_total(combined_adata, target_sum=1e4)
sc.pp.log1p(combined_adata)

# Highly variable genes
sc.pp.highly_variable_genes(combined_adata, flavor='seurat', n_top_genes=2000)

# Scale
sc.pp.scale(combined_adata, max_value=10)

# PCA and neighbors
sc.tl.pca(combined_adata, svd_solver='arpack')
sc.pp.neighbors(combined_adata)

# UMAP
sc.tl.umap(combined_adata)
sc.pl.umap(combined_adata, color='sample', size=2, save=f'_{basename}.png')
# Save processed object
combined_adata.write(newObject)

