#!/usr/bin/env python3
import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import numpy as np
import os

# Fix importlib for older packages
sys.modules['importlib.metadata'] = importlib_metadata

# --- Parse input ---
parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()
myObject = args.myObject
newObject = "clustered.h5ad"  
basename ="clustered" 

# --- Load data ---
combined_adata = sc.read(myObject)

# --- Clean NaNs / Infs safely ---
if hasattr(combined_adata, "X") and combined_adata.X is not None:
    combined_adata.X = np.nan_to_num(combined_adata.X, nan=0, posinf=0, neginf=0)

# --- Leiden clustering ---
sc.tl.leiden(combined_adata, resolution=2.5, flavor="igraph", n_iterations=2, directed=False)


# --- Ensure figures are saved as PNG at 600 dpi ---
sc.settings.figdir = './figures'
sc.settings.file_format_figs = 'png'
plt.rcParams['savefig.dpi'] = 600

# --- UMAP plot for clusters ---
sc.pl.umap(combined_adata, color=["leiden"], save="_clusters.png", legend_loc="on data")

# --- Marker genes dictionary ---
marker_genes  = {
    "MG": ["Rlbp1","Gfap","Apoe","Notch1","Pax6","Slc1a3","Vim"],
    "Rod": ["Rho","Nrl","Crx","Rom1"],
    "Cones": ["Opn1mw","Opn1sw","Arr3","Thrb","Gnat2"],
    "BC": ["Vsx1", "Sebox","Bhlhe23","Cabp5","Vsx1","Pcp4","Isl1"] ,
    "AC": ["Gad1","Gad2","Slc6a9","Tfap2b","Prox1","Pax6","Calb2","Pcp4","Elavl3","Isl1"],
    "HC": ["Lhx1","Cbln4","Calb1","Nefl","Nefm", "Onecut1", "Onecut2"],
    "RGC": ["Nefl","Nefm","Sncg","Thy1","Ebf3","Rbfox3","Isl1","Isl2","Pou4f1","Pou4f3","Rbpms"],
    "Microglia": ["Ptprc","Csf2rb","Sall1"],
    "Astrocytes":["Pax2","Igf2", "Gfap"]
}

# --- Filter marker genes to keep only those present in data ---
filtered_marker_genes = {
    cluster: [gene for gene in genes if gene in combined_adata.var_names]
    for cluster, genes in marker_genes.items()
}

# Optional: print missing genes for awareness
for cluster, genes in marker_genes.items():
    missing = set(genes) - set(filtered_marker_genes[cluster])
    if missing:
        print(f"⚠️  Missing genes for cluster '{cluster}': {missing}")


# --- Scatter plots for each gene ---
for cluster, genes in filtered_marker_genes.items():
    for gene in genes:
        sc.pl.scatter(combined_adata, color=gene, title=f'{cluster} - {gene}', basis='umap', show=False)
        plt.savefig(os.path.join(sc.settings.figdir, f'{gene}.png'), dpi=600)
        plt.close()

# --- Dotplot ---

# --- Dotplot (normal) ---
sc.pl.dotplot(
    combined_adata,
    filtered_marker_genes,
    groupby="leiden",
    standard_scale="var",
    save=f'{basename}_markerGenes.png'
)

# --- Save processed AnnData ---
combined_adata.obs_names_make_unique()
combined_adata.write(newObject, compression="gzip")

