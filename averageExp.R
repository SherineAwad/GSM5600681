library(Seurat)
library(GenomeInfoDb)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)

set.seed(1234)

# ---- Input argument ----
args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]

# ---- Load analysed Seurat object ----
myRDS <- paste0(mysample, "_annotated.rds")
cat("Reading:", myRDS, "\n")
myObject <- readRDS(myRDS)

# ---- Set RNA assay ----
DefaultAssay(myObject) <- "RNA"

# ---- Calculate average expression per cluster ----
cat("Calculating average expression (RNA)...\n")
avg_exp <- AverageExpression(
  object = myObject,
  assays = "RNA",
  slot = "data"   # normalized data; use "counts" for raw counts if preferred
)

# ---- Extract and save average expression ----
avg_exp_df <- as.data.frame(avg_exp$RNA)
head(avg_exp_df)

file_name <- paste0(mysample, "_AverageExpression.csv")
write.csv(avg_exp_df, file = file_name)
cat("Average expression saved to:", file_name, "\n")

# ---- Save updated object ----
myRDS_out <- paste0(mysample, "_avgExpression.rds")
saveRDS(myObject, file = myRDS_out)
cat("Saved Seurat object to:", myRDS_out, "\n")

# ---- Optional summary message ----
cat("✅ Finished computing average expression per cluster for sample:", mysample, "\n")

