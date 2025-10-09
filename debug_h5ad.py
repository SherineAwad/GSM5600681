#!/usr/bin/env python3
import scanpy as sc
import argparse
import pandas as pd
import numpy as np
import os

# --- Parse input argument ---
parser = argparse.ArgumentParser(description="Inspect contents of an h5ad file")
parser.add_argument("h5ad_file", help="Path to the .h5ad file")
args = parser.parse_args()

file_path = args.h5ad_file

# --- Load AnnData ---
print(f"Reading h5ad file: {file_path}")
adata = sc.read(file_path)

# --- Basic info ---
print("\n=== Basic Info ===")
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")
print(f"Shape of X: {adata.X.shape}")
print(f"Memory usage (approx.): {adata.n_obs * adata.n_vars * 8 / 1e6:.2f} MB")

# --- .obs information ---
print("\n=== obs (cell metadata) ===")
print(f"Columns: {adata.obs.columns.tolist()}")
print("First 10 rows:")
print(adata.obs.head(10))

# --- .var information ---
print("\n=== var (gene metadata) ===")
print(f"Columns: {adata.var.columns.tolist()}")
print("First 10 rows:")
print(adata.var.head(10))

# --- .obsm / .varm keys ---
print("\n=== obsm keys (embeddings) ===")
print(list(adata.obsm.keys()))

print("\n=== varm keys (gene embeddings) ===")
print(list(adata.varm.keys()))

# --- .uns keys (misc info) ===
print("\n=== uns keys ===")
print(list(adata.uns.keys()))

# --- Gene list ---
print("\n=== Gene Names ===")
print(f"First 50 genes: {adata.var_names[:50].tolist()}")

# --- Check for duplicated gene names ---
dupes = adata.var_names[adata.var_names.duplicated()]
print(f"Duplicated genes: {dupes.tolist() if len(dupes) > 0 else 'None'}")

# --- Check for NaN or Inf in X ---
if hasattr(adata, "X"):
    X = adata.X
    has_nan = np.any(np.isnan(X))
    has_inf = np.any(np.isinf(X))
    print(f"\nContains NaNs in X? {has_nan}")
    print(f"Contains Infs in X? {has_inf}")

# --- Summary of expression stats ---
print("\n=== Expression summary (first 10 genes) ===")
if hasattr(adata, "X"):
    if hasattr(X, "toarray"):  # sparse matrix
        X_dense = X.toarray()
    else:
        X_dense = X
    df_summary = pd.DataFrame({
        "mean": np.mean(X_dense[:, :10], axis=0),
        "std": np.std(X_dense[:, :10], axis=0),
        "min": np.min(X_dense[:, :10], axis=0),
        "max": np.max(X_dense[:, :10], axis=0),
        "nonzero_counts": np.sum(X_dense[:, :10] != 0, axis=0)
    }, index=adata.var_names[:10])
    print(df_summary)

# --- Optional: Save lists for inspection ---
gene_list_file = os.path.splitext(file_path)[0] + "_genes.txt"
cell_list_file = os.path.splitext(file_path)[0] + "_cells.txt"
adata.var_names.to_series().to_csv(gene_list_file, index=False)
adata.obs_names.to_series().to_csv(cell_list_file, index=False)
print(f"\nSaved gene list to {gene_list_file}")
print(f"Saved cell list to {cell_list_file}")

print("\n=== Debugging complete ===")

