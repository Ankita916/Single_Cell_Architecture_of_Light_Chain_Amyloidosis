# 02_integration.R
# Batch integration of the merged object using RPCA (Seurat v5),
# followed by clustering and UMAP on the integrated embedding.
#
# Input:  results/objects/merged_seurat.rds
# Output: results/objects/merged_seurat_integrated.rds

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(1234)

# Parameters
obj_dir     <- "results/objects"
n_dims      <- 21
cluster_res <- 0.8

options(future.globals.maxSize = 16 * 1024^3)

merged_seurat <- readRDS(file.path(obj_dir, "merged_seurat.rds"))

# RPCA integration across samples
merged_seurat_int <- IntegrateLayers(
  object         = merged_seurat,
  method         = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.rpca",
  verbose        = FALSE
)

merged_seurat_int[["RNA"]] <- JoinLayers(merged_seurat_int[["RNA"]])

# Clustering and UMAP on the integrated embedding
merged_seurat_int <- FindNeighbors(merged_seurat_int, dims = 1:n_dims,
                                   reduction = "integrated.rpca")
merged_seurat_int <- FindClusters(merged_seurat_int,
                                  cluster.name = "rpca_clusters",
                                  resolution = cluster_res)
merged_seurat_int <- RunUMAP(merged_seurat_int, dims = 1:n_dims,
                             reduction = "integrated.rpca",
                             reduction.name = "umap.rpca")

saveRDS(merged_seurat_int, file.path(obj_dir, "merged_seurat_integrated.rds"))
