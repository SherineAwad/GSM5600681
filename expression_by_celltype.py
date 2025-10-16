import argparse
import scanpy as sc
import pandas as pd
import numpy as np

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Compute gene expression per cell type.")
    parser.add_argument("--input", "-i", required=True, help="Path to input .h5ad file")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file name")
    args = parser.parse_args()

    # Load the h5ad file
    print(f"🔍 Loading {args.input}...")
    adata = sc.read_h5ad(args.input)

    # Check if 'celltype' is present in .obs
    if 'celltype' not in adata.obs.columns:
        raise ValueError("The AnnData object does not contain a 'celltype' column in .obs")

    # Convert to dense if needed
    if not isinstance(adata.X, (pd.DataFrame, np.ndarray)):
        print("⚠️  Converting sparse matrix to dense (may use significant memory)...")
        adata.X = adata.X.toarray()

    # Create DataFrame with gene expression, indexed by celltype
    print("📊 Grouping cells by celltype and computing mean expression...")
    expr_df = pd.DataFrame(adata.X, index=adata.obs['celltype'], columns=adata.var_names)
    mean_expr_by_celltype = expr_df.groupby(level=0).mean()

    # Transpose so that genes are rows, cell types are columns
    result_df = mean_expr_by_celltype.T

    # Save to CSV
    print(f"💾 Saving output to {args.output}...")
    result_df.to_csv(args.output)
    print("✅ Done.")

if __name__ == "__main__":
    main()

