import scanpy as sc
import doubletdetection
import numpy as np
from scipy import sparse
import argparse
import os
import matplotlib.pyplot as plt

# Parse input file
parser = argparse.ArgumentParser()
parser.add_argument('myObject')
parser.add_argument('threshold')

args = parser.parse_args()
myObject = args.myObject
threshold = float(args.threshold) 
# Output file name
newObject = "noDoublets.h5ad" 

base_name = os.path.splitext(os.path.basename(newObject))[0]

# Read AnnData object
adata = sc.read_h5ad(myObject)

# Clean NaNs in adata.X
if sparse.issparse(adata.X):
    print("Cleaning sparse matrix (adata.X)...")
    adata.X.data[np.isnan(adata.X.data)] = 0  # Only fix non-zero values
    check_matrix = adata.X.toarray()
else:
    print("Cleaning dense matrix (adata.X)...")
    adata.X = np.nan_to_num(adata.X)
    check_matrix = adata.X

# Sanity check: ensure no NaNs remain
if np.isnan(check_matrix).any():
    raise ValueError("NaNs still present in adata.X after cleanup.")

# Run doublet detection
print("Running doublet detection...")
clf = doubletdetection.BoostClassifier(n_iters=5, standard_scaling=True)
doublets = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=threshold)



# Threshold plot
figname1 = "figures/"+args.threshold + "threshold_test.png"
print(figname1)
fig1 = doubletdetection.plot.threshold(clf, show=False, p_step=6)
fig1.savefig(figname1, dpi=300)
plt.close(fig1)

# Convergence plot
figname2 = "figures/"+args.threshold + "conversion_test.png"
print(figname2)
fig2 = doubletdetection.plot.convergence(clf, show=False, p_thresh=1e-16, voter_thresh=threshold)
fig2.savefig(figname2, dpi=300)
plt.close(fig2)


# Save results to AnnData
adata.obs['doublet_score'] = clf.doublet_score()
adata.obs['predicted_doublet'] = doublets

# Save updated AnnData object
adata.write(newObject, compression="gzip")
print("Doublet detection completed and results saved.")

