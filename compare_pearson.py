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
                    help="Minimum average expression threshold for genes (default: 0.1)")
parser.add_argument("--output", default="correlation_celltypes.png", help="Output heatmap file")
args = parser.parse_args()

# -------------------------------
# LOAD DATA
# -------------------------------
print("🔹 Loading data...")
adata1 = sc.read_h5ad(args.h5ad1)
adata2 = sc.read_h5ad(args.h5ad2)

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
# GENE FILTERING
# -------------------------------
print("🔹 Filtering genes based on expression cutoff and presence...")
# Mean expression across all cell types in both datasets
mean_expr = pd.concat([expr1_ct.mean(axis=0), expr2_ct.mean(axis=0)], axis=1).mean(axis=1)

# Expression presence (gene expressed in at least one cell type)
expressed_in_ct = (
    ((expr1_ct > 0).sum(axis=0) > 0) |
    ((expr2_ct > 0).sum(axis=0) > 0)
)

# Keep genes that pass both filters
genes_to_keep = mean_expr[(mean_expr > args.gene_cutoff) & expressed_in_ct].index
print(f"✅ Retained {len(genes_to_keep)} genes out of {len(common_genes)} total.")

expr1_ct = expr1_ct.loc[:, genes_to_keep]
expr2_ct = expr2_ct.loc[:, genes_to_keep]

# -------------------------------
# STANDARDIZE ROWS (cell types)
# -------------------------------
print("🔹 Standardizing expression (z-score per celltype)...")
X1_std = expr1_ct.sub(expr1_ct.mean(axis=1), axis=0).div(expr1_ct.std(axis=1), axis=0)
X2_std = expr2_ct.sub(expr2_ct.mean(axis=1), axis=0).div(expr2_ct.std(axis=1), axis=0)

# Replace NaNs (if any celltype has 0 variance)
X1_std = X1_std.fillna(0)
X2_std = X2_std.fillna(0)

# -------------------------------
# VECTORISED PEARSON CORRELATION
# -------------------------------
print("🔹 Computing Pearson correlation between celltypes...")
corr_matrix = np.dot(X1_std.values, X2_std.values.T) / X1_std.shape[1]
corr_matrix = pd.DataFrame(corr_matrix, index=expr1_ct.index, columns=expr2_ct.index)

# -------------------------------
# PLOT HEATMAP
# -------------------------------
print(f"🔹 Plotting heatmap and saving to {args.output}...")
plt.figure(figsize=(8, 6))
sns.heatmap(corr_matrix, cmap='vlag', center=0, annot=True, fmt=".2f")
plt.xlabel(f"Celltypes in {args.sample2}")
plt.ylabel(f"Celltypes in {', '.join(args.samples1)}")
plt.title("Pearson correlation between celltypes (filtered genes)")
plt.tight_layout()
plt.savefig(args.output, dpi=600)

print("✅ Done! Correlation heatmap saved.")

