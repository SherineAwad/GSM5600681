#!/usr/bin/env python3
import scanpy as sc
import pandas as pd
import argparse

# -------------------------------
# PARSE ARGUMENTS
# -------------------------------
parser = argparse.ArgumentParser(
    description="Export mean expression per celltype from h5ad object"
)
parser.add_argument("h5ad", help="Path to h5ad object")
parser.add_argument("output", help="Output CSV file")
parser.add_argument("--samples", nargs='+', required=True, help="Sample names to include")
args = parser.parse_args()

# -------------------------------
# LOAD DATA
# -------------------------------
print("🔹 Loading data...")
adata = sc.read_h5ad(args.h5ad)

# -------------------------------
# SUBSET SAMPLES
# -------------------------------
print("🔹 Subsetting samples...")
adata_subset = adata[adata.obs['sample'].isin(args.samples)].copy()
print(f"✅ Retained {adata_subset.n_obs} cells from samples: {args.samples}")

# -------------------------------
# AVERAGE EXPRESSION PER CELLTYPE
# -------------------------------
print("🔹 Calculating mean expression per cell type...")
expr_ct = adata_subset.to_df().groupby(adata_subset.obs['celltype'], observed=False).mean()

# -------------------------------
# EXPORT - TRANSPOSE to get genes as rows
# -------------------------------
print(f"🔹 Exporting to {args.output}...")
expr_ct.T.to_csv(args.output)
print(f"✅ Exported {expr_ct.T.shape[0]} genes × {expr_ct.T.shape[1]} celltypes")
print("✅ Done!")
