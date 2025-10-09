import scanpy as sc
import doubletdetection
import numpy as np
from scipy import sparse
import argparse
import os
import matplotlib.pyplot as plt

# -----------------------------
# Parse input arguments
# -----------------------------
parser = argparse.ArgumentParser(description="Detect doublets and plot UMAPs for all suggested thresholds")
parser.add_argument('myObject', help="Input AnnData h5ad file")
args = parser.parse_args()
myObject = args.myObject

# -----------------------------
# Output setup
# -----------------------------
fig_dir = "figures"
os.makedirs(fig_dir, exist_ok=True)
newObject = "markedDoublets.h5ad"

# -----------------------------
# Read AnnData
# -----------------------------
adata = sc.read_h5ad(myObject)

# Clean NaNs
if sparse.issparse(adata.X):
    adata.X.data[np.isnan(adata.X.data)] = 0
else:
    adata.X = np.nan_to_num(adata.X)

# -----------------------------
# Compute UMAP before doublets
# -----------------------------
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=None, show=False, save="_UMAP.png")

# -----------------------------
# Run doublet detection
# -----------------------------
clf = doubletdetection.BoostClassifier(n_iters=5, standard_scaling=True)
clf.fit(adata.X)

# Plot threshold curve to get suggested thresholds
fig_threshold = doubletdetection.plot.threshold(clf, show=False, p_step=6)
fig_threshold.savefig(os.path.join(fig_dir, "threshold_curve.png"), dpi=300)
plt.close(fig_threshold)

# doubletdetection does not return thresholds directly,
# but you can generate thresholds from the range shown in the figure.
# For example, sweep from 0.1 to 0.9 by 0.1 (adjust as needed):
suggested_thresholds = np.arange(0.1, 1.0, 0.1)

# -----------------------------
# Loop over all thresholds
# -----------------------------
for thresh in suggested_thresholds:
    doublets = clf.predict(p_thresh=1e-16, voter_thresh=thresh)
    adata.obs[f'predicted_doublet_{thresh:.2f}'] = doublets
    sc.pl.umap(adata, color=f'predicted_doublet_{thresh:.2f}', show=False,
               save=f"_UMAP_doublets_thresh{thresh:.2f}.png")

# -----------------------------
# Save doublet scores
# -----------------------------
adata.obs['doublet_score'] = clf.doublet_score()
adata.write(newObject, compression="gzip")
print(f"Doublet detection completed. Results saved to {newObject}")

