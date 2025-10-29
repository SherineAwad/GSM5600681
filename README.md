# Analysis of adult Retina WT sample 
# GSE184933: GSM5600681



##  Quality Control Plots

These violin plots visualize cell-level quality metrics **before and after filtering** low-quality cells.

---

###  `violin_QC.png`

This plot shows **raw quality control metrics** (e.g., number of genes per cell, mitochondrial percentage) **before filtering**. It helps ide.pngy thresholds for removing low-quality or dead cells.

![violin_QC.png](figures/violin_QC.png)

---

### Filtering 

We filter cells based on basic quality control metrics:  

- Keep cells with **800–8000 detected genes**  
- Keep cells with **1200–30000 total counts**  
- Remove cells with **high mitochondrial content (>25%)**  

```python
combined_adata = combined_adata[
    (combined_adata.obs['n_genes_by_counts'] > 800) &
    (combined_adata.obs['n_genes_by_counts'] < 8000) &
    (combined_adata.obs['total_counts'] > 1200) &
    (combined_adata.obs['total_counts'] < 30000) &
    (combined_adata.obs['pct_counts_mt'] < 25),
    :
]
``` 

###  `violin_AfterQC.png`

This plot shows the same QC metrics **after filtering**. It confirms that poor-quality cells were successfully removed based on chosen thresholds.

![violin_AfterQC.png](figures/violin_AfterQC.png)


## Sample UMAP 
![sample umap](figures/umap_WT.png)


## Clusters UMAP 
![clusters umap](figures/umap_clusters.png?v=3)


## Marker Genes Dotplot 
![markers dotplot](figures/dotplot_clustered_markerGenes.png)

## Marker Genes UMAP 

<img src="figures/Nrl.png" alt="Nrl" width="33%"><img src="figures/Crx.png" alt="Crx" width="33%"><img src="figures/Pcp4.png" alt="Pcp4" width="33%">
<img src="figures/Rom1.png" alt="Rom1" width="33%"><img src="figures/Rho.png" alt="Rho" width="33%"><img src="figures/Apoe.png" alt="Apoe" width="33%">
<img src="figures/Calb2.png" alt="Calb2" width="33%"><img src="figures/Arr3.png" alt="Arr3" width="33%"><img src="figures/Thrb.png" alt="Thrb" width="33%">
<img src="figures/Opn1sw.png" alt="Opn1sw" width="33%"><img src="figures/Prox1.png" alt="Prox1" width="33%"><img src="figures/Pax6.png" alt="Pax6" width="33%">
<img src="figures/Nefl.png" alt="Nefl" width="33%"><img src="figures/Sncg.png" alt="Sncg" width="33%"><img src="figures/Gnat2.png" alt="Gnat2" width="33%">
<img src="figures/Cabp5.png" alt="Cabp5" width="33%"><img src="figures/Rlbp1.png" alt="Rlbp1" width="33%"><img src="figures/Gad1.png" alt="Gad1" width="33%">
<img src="figures/Isl1.png" alt="Isl1" width="33%"><img src="figures/Opn1mw.png" alt="Opn1mw" width="33%"><img src="figures/Tfap2b.png" alt="Tfap2b" width="33%">
<img src="figures/Vim.png" alt="Vim" width="33%"><img src="figures/Elavl3.png" alt="Elavl3" width="33%"><img src="figures/Slc6a9.png" alt="Slc6a9" width="33%">
<img src="figures/Nefm.png" alt="Nefm" width="33%"><img src="figures/Calb1.png" alt="Calb1" width="33%"><img src="figures/Thy1.png" alt="Thy1" width="33%">
<img src="figures/Gad2.png" alt="Gad2" width="33%"><img src="figures/Sebox.png" alt="Sebox" width="33%"><img src="figures/Bhlhe23.png" alt="Bhlhe23" width="33%">
<img src="figures/Slc1a3.png" alt="Slc1a3" width="33%"><img src="figures/Vsx1.png" alt="Vsx1" width="33%"><img src="figures/Rbfox3.png" alt="Rbfox3" width="33%">
<img src="figures/Notch1.png" alt="Notch1" width="33%"><img src="figures/Igf2.png" alt="Igf2" width="33%"><img src="figures/Rbpms.png" alt="Rbpms" width="33%">
<img src="figures/Ebf3.png" alt="Ebf3" width="33%"><img src="figures/Onecut2.png" alt="Onecut2" width="33%"><img src="figures/Pou4f1.png" alt="Pou4f1" width="33%">
<img src="figures/Onecut1.png" alt="Onecut1" width="33%"><img src="figures/Isl2.png" alt="Isl2" width="33%"><img src="figures/Pou4f3.png" alt="Pou4f3" width="33%">
<img src="figures/Lhx1.png" alt="Lhx1" width="33%"><img src="figures/Gfap.png" alt="Gfap" width="33%"><img src="figures/Pax2.png" alt="Pax2" width="33%">
<img src="figures/Csf2rb.png" alt="Csf2rb" width="33%"><img src="figures/Cbln4.png" alt="Cbln4" width="33%"><img src="figures/Sall1.png" alt="Sall1" width="33%">
<img src="figures/Ptprc.png" alt="Ptprc" width="33%"><img src="figures/Chat.png" alt="Chat" width="33%"><img src="figures/Rpe65.png" alt="Rpe65"  width="33%">
<img src="figures/Acta2.png" alt="Acta2" width="33%"> <img src="figures/Otx2.png" alt="Otx2" width="33%">


## Annotations 


![annotations](figures/WT_annotations.png?v=4)

![annotationsON](figures/WT_annotationsON.png?v=4)


## 🔬 Comparing Neurog2 (Neurog2_9SA_5weeks, Neurog2_9SA_2mo) with WT

We used the script **`compare_pearson.py`** to evaluate transcriptional similarity between cell types in the Neurog2 samples and those in the WT sample.

The neurog2 samples are analysed in [Neurog2](https://github.com/SherineAwad/Neurog2)

### ⚙️  What `compare_pearson.py` does?

```
1. **Input and subsetting**  
   The analysis uses two `.h5ad` files — one for the Neurog2 samples and one for the WT sample.  
   The script selects only the specified samples and aligns both datasets to the same set of genes.

2. **Averaging gene expression per cell type**  
   For each dataset, the mean expression level of every gene is calculated within each cell type.  
   This step produces one representative expression profile per cell type.

3. **Filtering genes**  
   Genes are kept only if their average expression is above a user-defined cutoff and they are expressed in at least one cell.  
   This ensures that the analysis focuses on biologically meaningful and active genes.

4. **Pearson correlation between cell types**  
   The script then compares every Neurog2 cell type to every WT cell type by calculating the Pearson correlation between their average expression profiles.  
   This measures how similar the expression patterns are between the two cell types across all retained genes.
    
5. **Heatmap visualization**  
   The resulting correlation matrix is shown as a heatmap:  
   - The **rows** correspond to Neurog2 cell types.  
   - The **columns** correspond to WT cell types.  
   - The **color intensity** reflects the correlation strength — higher values indicate more similar transcriptional profiles.
```


#### Using cutoff 0.05 
![0.05](heatmap05.png)

```
python compare_pearson.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad annotated.h5ad --samples1 Neurog2_9SA_5weeks Neurog2_9SA_2mo --sample2 WT  --output heatmap05.png --gene_cutoff 0.05 
✅ Retained 17287 genes out of 19780 total.
``` 

#### Using cutoff 0.1 
![0.1](heatmap_1.png)

```
python compare_pearson.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad annotated.h5ad --samples1 Neurog2_9SA_5weeks Neurog2_9SA_2mo --sample2 WT  --output heatmap_1.png --gene_cutoff 0.1 
✅ Retained 14841 genes out of 19780 total.
``` 

#### Using cutoff 0.2 
![0.2](heatmap_2.png)

```
python compare_pearson.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad annotated.h5ad --samples1 Neurog2_9SA_5weeks Neurog2_9SA_2mo --sample2 WT  --output heatmap_1.png --gene_cutoff 0.2
✅ Retained 11131 genes out of 19780 total.
```

#### Using cutoff 0.5 
![0.5](heatmap_5.png)
```
python compare_pearson.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad annotated.h5ad --samples1 Neurog2_9SA_5weeks Neurog2_9SA_2mo --sample2 WT  --output heatmap_1.png --gene_cutoff 0.
✅ Retained 5437 genes out of 19780 total.
```


# 🧩 For Debugging Issues

```
# Name                     Version          Build            Channel
scanpy                     1.10.3           pypi_0           pypi
```

### 💾 Command & Output for Neuorg2
```bash
python print_h5ad.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad 
/nfs/turbo/umms-thahoang/sherine/miniconda/envs/scanpy_solo_env/lib/python3.9/site-packages/louvain/__init__.py:54: UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools<81.
  from pkg_resources import get_distribution, DistributionNotFound
Loading AnnData object from: annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad

==================================================
BASIC AnnData INFORMATION
==================================================
Overall shape: (39764, 39597) (cells: 39764, genes: 39597)
Observation names (cells): ['AAACCAAAGCCATACA-1', 'AAACCCGCAATCCGTC-1', 'AAACCCGCATCACTGC-1', 'AAACCCGCATCGTACC-1', 'AAACCCTGTTGCTGTG-1']...
Variable names (genes): ['Xkr4', 'Gm1992', 'Gm19938', 'ENSMUSG00000142494', 'Gm37381']...

=== adata.X (Primary Expression Matrix) ===
Shape: (39764, 39597) (cells: 39764, genes: 39597)
Type: dense
Data type: float32
Min: -6.156065
Max: 199.404129
Mean: 0.000000
Std: 0.999985
Non-zero values: 1574535108
Zero values: 0
Sparsity: 0.0000 (0.00%)

Top 5 genes (first 5 rows):
  Xkr4: min=-0.7485, max=4.9478, mean=0.0000
  Gm1992: min=-0.1218, max=10.4576, mean=0.0000
  Gm19938: min=-0.2123, max=9.5976, mean=0.0000
  ENSMUSG00000142494: min=-0.0188, max=53.2843, mean=-0.0000
  Gm37381: min=-0.0230, max=43.5026, mean=-0.0000

Top 5 cells (first 5 columns):
  AAACCAAAGCCATACA-1: min=-2.2320, max=45.7361, mean=-0.0210
  AAACCCGCAATCCGTC-1: min=-1.8717, max=31.5131, mean=0.0172
  AAACCCGCATCACTGC-1: min=-2.4021, max=115.1231, mean=0.0330
  AAACCCGCATCGTACC-1: min=-2.7288, max=99.6983, mean=-0.0177
  AAACCCTGTTGCTGTG-1: min=-1.8227, max=75.3620, mean=0.0293

adata.raw type: <class 'anndata._core.raw.Raw'>

=== adata.raw.X (Raw Expression Matrix) ===
Shape: (39764, 40007) (cells: 39764, genes: 40007)
Type: dense
Data type: float32
Min: -4.656951
Max: 10.000000
Mean: -0.007822
Std: 0.693484
Non-zero values: 1590838348
Zero values: 0
Sparsity: 0.0000 (0.00%)

Top 5 genes (first 5 rows):
  Xkr4: min=-0.7183, max=5.1481, mean=0.0526
  Gm1992: min=-0.1151, max=10.0000, mean=0.0014
  Gm19938: min=-0.2098, max=10.0000, mean=0.0111
  ENSMUSG00000142494: min=-0.0205, max=10.0000, mean=-0.0170
  Gm37381: min=-0.0204, max=10.0000, mean=-0.0151

Top 5 cells (first 5 columns):
  AAACCAAAGCCATACA-1: min=-2.1655, max=10.0000, mean=-0.0234
  AAACCCGCAATCCGTC-1: min=-1.9220, max=10.0000, mean=0.0131
  AAACCCGCATCACTGC-1: min=-2.2012, max=10.0000, mean=0.0207
  AAACCCGCATCGTACC-1: min=-2.9586, max=10.0000, mean=-0.0216
  AAACCCTGTTGCTGTG-1: min=-1.8224, max=10.0000, mean=0.0171

==================================================
COMPARISON: adata.X vs adata.raw.X
==================================================
Same shape: False

==================================================
CATEGORICAL OBSERVATIONS
==================================================
Found 3 categorical columns:

Column: sample (3 unique values)
Values: ['control_2mo', 'Neurog2_9SA_2mo', 'Neurog2_9SA_5weeks']
---
Column: leiden (51 unique values)
Values: ['0', '6', '3', '10', '5', '8', '15', '4', '1', '9', '13', '18', '2', '19', '20', '16', '21', '22', '12', '23', '24', '14', '17', '27', '28', '30', '34', '35', '36', '37', '38', '39', '40', '32', '7', '41', '29', '42', '43', '25', '44', '26', '45', '11', '31', '46', '47', '48', '49', '50', '33']
---
Column: celltype (6 unique values)
Values: ['MG', 'Rod', 'Cones', 'BC', 'MGPC', 'AC']
---

==================================================
ADDITIONAL METADATA
==================================================
adata.obs columns: ['sample', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'n_genes', 'doublet_score', 'predicted_doublet', 'leiden', 'celltype']
adata.var columns: ['mt', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'mean', 'std']

No additional layers found.
``` 

### 💾 Command & Output for WT
```bash
python print_h5ad.py annotated.h5ad 
/nfs/turbo/umms-thahoang/sherine/miniconda/envs/scanpy_solo_env/lib/python3.9/site-packages/louvain/__init__.py:54: UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools<81.
  from pkg_resources import get_distribution, DistributionNotFound
Loading AnnData object from: annotated.h5ad

==================================================
BASIC AnnData INFORMATION
==================================================
Overall shape: (21256, 21079) (cells: 21256, genes: 21079)
Observation names (cells): ['AAACCCAAGACTTAAG-1', 'AAACCCAAGCATTTCG-1', 'AAACCCAAGGCACGAT-1', 'AAACCCAAGGGTGGGA-1', 'AAACCCAAGGTGCGAT-1']...
Variable names (genes): ['Xkr4', 'Gm1992', 'Gm19938', 'Gm37381', 'Rp1']...

=== adata.X (Primary Expression Matrix) ===
Shape: (21256, 21079) (cells: 21256, genes: 21079)
Type: dense
Data type: float32
Min: -10.000000
Max: 10.000000
Mean: -0.010800
Std: 0.814151
Non-zero values: 448055224
Zero values: 0
Sparsity: 0.0000 (0.00%)

Top 5 genes (first 5 rows):
  Xkr4: min=-0.1643, max=10.0000, mean=0.0038
  Gm1992: min=-0.0211, max=10.0000, mean=-0.0178
  Gm19938: min=-0.1757, max=10.0000, mean=0.0054
  Gm37381: min=-0.0157, max=10.0000, mean=-0.0134
  Rp1: min=-1.6397, max=2.1227, mean=0.0249

Top 5 cells (first 5 columns):
  AAACCCAAGACTTAAG-1: min=-3.1532, max=10.0000, mean=-0.0704
  AAACCCAAGCATTTCG-1: min=-3.2957, max=10.0000, mean=-0.0572
  AAACCCAAGGCACGAT-1: min=-1.7961, max=10.0000, mean=0.0169
  AAACCCAAGGGTGGGA-1: min=-3.2656, max=10.0000, mean=-0.1094
  AAACCCAAGGTGCGAT-1: min=-2.2016, max=10.0000, mean=0.0160

adata.raw type: <class 'anndata._core.raw.Raw'>

=== adata.raw.X (Raw Expression Matrix) ===
Shape: (21256, 21079) (cells: 21256, genes: 21079)
Type: sparse
Data type: float32
Min: 0.000000
Max: 8.622109
Mean: 0.138170
Std: 0.470614
Non-zero values: 42633255
Zero values: 405421969
Sparsity: 0.9048 (90.48%)

Top 5 genes (first 5 rows):
  Xkr4: min=0.0000, max=2.8670, mean=0.0333
  Gm1992: min=0.0000, max=0.7869, mean=0.0003
  Gm19938: min=0.0000, max=3.1630, mean=0.0392
  Gm37381: min=0.0000, max=2.0586, mean=0.0003
  Rp1: min=0.0000, max=4.8305, mean=2.1371

Top 5 cells (first 5 columns):
  AAACCCAAGACTTAAG-1: min=0.0000, max=6.4226, mean=0.1057
  AAACCCAAGCATTTCG-1: min=0.0000, max=6.3472, mean=0.1179
  AAACCCAAGGCACGAT-1: min=0.0000, max=5.9292, mean=0.1594
  AAACCCAAGGGTGGGA-1: min=0.0000, max=5.8897, mean=0.0955
  AAACCCAAGGTGCGAT-1: min=0.0000, max=6.0853, mean=0.1577

==================================================
COMPARISON: adata.X vs adata.raw.X
==================================================
Same shape: True
Number of different elements: 448055221
Identical matrices: False

==================================================
CATEGORICAL OBSERVATIONS
==================================================
Found 3 categorical columns:

Column: sample (1 unique values)
Values: ['WT']
---
Column: leiden (55 unique values)
Values: ['0', '2', '1', '12', '14', '7', '15', '16', '17', '6', '21', '23', '4', '5', '18', '8', '25', '26', '27', '11', '29', '9', '30', '13', '31', '32', '3', '35', '36', '37', '38', '39', '40', '41', '42', '43', '24', '44', '45', '46', '48', '49', '51', '52', '53', '54', '57', '58', '28', '59', '60', '19', '61', '62', '63']
---
Column: celltype (7 unique values)
Values: ['BC', 'Rod', 'Cones', 'MG', 'AC', 'RGC', 'HC']
---

==================================================
ADDITIONAL METADATA
==================================================
adata.obs columns: ['sample', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'n_genes', 'leiden', 'celltype']
adata.var columns: ['mt', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'highly_variable', 'means',]()

