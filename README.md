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


#### Printing expressions for diagnoses 
```
python print_expression.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad neurog2.csv --samples Neurog2_9SA_5weeks Neurog2_9SA_2mo 
```


[📥 Download neurog2 Expressions](https://docs.google.com/spreadsheets/d/1OH2wrqN7HCiSpkCQfB68l-ik--XXCsDwTcvmpyZV-mQ/edit?usp=sharing)
