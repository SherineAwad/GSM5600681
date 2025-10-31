#!/usr/bin/env python3
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

# -------------------------------
# PARSE ARGUMENTS
# -------------------------------
parser = argparse.ArgumentParser(
    description="Compare celltypes from two h5ad objects using Pearson correlation "
                "with gene filtering by average expression and presence."
)
parser.add_argument("h5ad1", help="Path to first h5ad object")
parser.add_argument("h5ad2", help="Path to second h5ad object")
parser.add_argument("--samples1", nargs='+', required=True, help="Sample names from first object")
parser.add_argument("--sample2", required=True, help="Sample name from second object")
parser.add_argument("--gene_cutoff", type=float, default=0.1,
                    help="Minimum average expression threshold for genes in any celltype (default: 0.1)")
parser.add_argument("--output", default="correlation_celltypes.png", help="Output heatmap file")
args = parser.parse_args()

# -------------------------------
# LOAD DATA
# -------------------------------
print("🔹 Loading data...")
adata1 = sc.read_h5ad(args.h5ad1)
adata2 = sc.read_h5ad(args.h5ad2)

# -------------------------------
# USE RAW DATA THEN NORMALIZE
# -------------------------------
print("🔹 Using raw data and normalizing...")
if adata1.raw is not None:
    adata1 = adata1.raw.to_adata()
if adata2.raw is not None:
    adata2 = adata2.raw.to_adata()

sc.pp.normalize_total(adata1, target_sum=1e4)
sc.pp.log1p(adata1)
sc.pp.normalize_total(adata2, target_sum=1e4)
sc.pp.log1p(adata2)

# -------------------------------
# SUBSET SAMPLES
# -------------------------------
print("🔹 Subsetting samples...")
adata1_subset = adata1[adata1.obs['sample'].isin(args.samples1)].copy()
adata2_subset = adata2[adata2.obs['sample'] == args.sample2].copy()

# -------------------------------
# ALIGN GENES
# -------------------------------
print("🔹 Aligning genes between datasets...")
common_genes = adata1_subset.var_names.intersection(adata2_subset.var_names)
adata1_subset = adata1_subset[:, common_genes]
adata2_subset = adata2_subset[:, common_genes]

# -------------------------------
# AVERAGE EXPRESSION PER CELLTYPE
# -------------------------------
print("🔹 Averaging expression per cell type...")
expr1_ct = adata1_subset.to_df().groupby(adata1_subset.obs['celltype'], observed=False).mean()
expr2_ct = adata2_subset.to_df().groupby(adata2_subset.obs['celltype'], observed=False).mean()

# -------------------------------
# GENE FILTERING - FIXED
# -------------------------------
print("🔹 Filtering genes based on expression cutoff and presence...")

expressed_in_any_ct = (
    (expr1_ct > args.gene_cutoff).any(axis=0) |
    (expr2_ct > args.gene_cutoff).any(axis=0)
)
present_in_any_ct = (
    (expr1_ct > 0).any(axis=0) |
    (expr2_ct > 0).any(axis=0)
)
genes_to_keep = expressed_in_any_ct & present_in_any_ct
genes_to_keep = genes_to_keep[genes_to_keep].index

print(f"✅ Retained {len(genes_to_keep)} genes out of {len(common_genes)} total.")

expr1_ct = expr1_ct.loc[:, genes_to_keep]
expr2_ct = expr2_ct.loc[:, genes_to_keep]

# -------------------------------
# PEARSON CORRELATION (without z-scoring)
# -------------------------------
print("🔹 Computing Pearson correlation between celltypes...")
corr_matrix = np.corrcoef(expr1_ct.values, expr2_ct.values)[:expr1_ct.shape[0], expr1_ct.shape[0]:]
corr_matrix = pd.DataFrame(corr_matrix, index=expr1_ct.index, columns=expr2_ct.index)

# -------------------------------
# PLOT HEATMAP (blue to red, no white)
# -------------------------------
print(f"🔹 Plotting heatmap and saving to {args.output}...")
plt.figure(figsize=(8, 6))
vmin = np.nanmin(corr_matrix.values)
vmax = 1
sns.heatmap(corr_matrix, cmap='coolwarm', annot=True, fmt=".2f", vmin=vmin, vmax=vmax)
plt.xlabel(f"Celltypes in {args.sample2}")
plt.ylabel(f"Celltypes in {', '.join(args.samples1)}")
plt.title("Pearson correlation between celltypes (filtered genes)")
plt.tight_layout()
plt.savefig(args.output, dpi=600)

print("✅ Done! Correlation heatmap saved.")

