# 03_batch_mixing.R
# Evaluation of batch mixing before and after integration:
#   (1) UMAP coloured by batch
#   (2) Local inverse Simpson's index (LISI; Korsunsky et al., 2019)
#   (3) Seurat MixingMetric (Stuart et al., 2019)
# Lower MixingMetric values indicate better mixing; higher LISI values
# indicate better mixing.
#
# Input:  results/objects/merged_seurat.rds
#         results/objects/merged_seurat_integrated.rds
# Output: results/figures/, results/tables/

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(FNN)
})

set.seed(1234)

# Parameters
obj_dir <- "results/objects"
fig_dir <- "results/figures"
tab_dir <- "results/tables"
n_dims  <- 21

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

merged_seurat     <- readRDS(file.path(obj_dir, "merged_seurat.rds"))
merged_seurat_int <- readRDS(file.path(obj_dir, "merged_seurat_integrated.rds"))

# UMAP coloured by batch
p_unintegrated <- DimPlot(merged_seurat, reduction = "umap",
                          group.by = "Batch")
p_integrated   <- DimPlot(merged_seurat_int, reduction = "umap.rpca",
                          group.by = "Batch")

ggsave(file.path(fig_dir, "Batch_effect_removal.png"),
       p_unintegrated + p_integrated,
       width = 2000, height = 1300, units = "px", dpi = 300)

# LISI computed per cell on an embedding
# X: cells x dimensions matrix; labels: batch label per cell
compute_lisi <- function(X, labels, perplexity = 30, tol = 1e-5, max_iter = 50) {
  labels   <- as.character(labels)
  k        <- min(nrow(X) - 1, perplexity * 3)
  target_H <- log(perplexity)
  nn       <- FNN::get.knn(X, k = k)

  sapply(seq_len(nrow(X)), function(i) {
    D    <- nn$nn.dist[i, ]^2
    beta <- 1
    bmin <- -Inf
    bmax <- Inf

    # Binary search for the kernel width matching the target perplexity
    for (it in seq_len(max_iter)) {
      P    <- exp(-D * beta)
      sumP <- sum(P)
      H    <- if (sumP == 0) 0 else log(sumP) + beta * sum(D * P) / sumP
      if (abs(H - target_H) < tol) break
      if (H > target_H) {
        bmin <- beta
        beta <- if (is.infinite(bmax)) beta * 2 else (beta + bmax) / 2
      } else {
        bmax <- beta
        beta <- if (is.infinite(bmin)) beta / 2 else (beta + bmin) / 2
      }
    }

    # Inverse Simpson's index of the batch composition of the neighbourhood
    P   <- P / sum(P)
    p_b <- tapply(P, labels[nn$nn.index[i, ]], sum)
    1 / sum(p_b^2, na.rm = TRUE)
  })
}

lisi_unintegrated <- compute_lisi(Embeddings(merged_seurat, "umap"),
                                  merged_seurat$Batch)
lisi_integrated   <- compute_lisi(Embeddings(merged_seurat_int, "umap.rpca"),
                                  merged_seurat_int$Batch)

lisi_summary <- data.frame(
  stage        = c("unintegrated", "integrated_rpca"),
  median_iLISI = round(c(median(lisi_unintegrated), median(lisi_integrated)), 3),
  mean_iLISI   = round(c(mean(lisi_unintegrated), mean(lisi_integrated)), 3)
)
write.csv(lisi_summary, file.path(tab_dir, "LISI_batch_mixing_summary.csv"),
          row.names = FALSE)
print(lisi_summary)

# Seurat mixing metric
mix_unintegrated <- MixingMetric(merged_seurat, grouping.var = "Batch",
                                 reduction = "pca", dims = 1:n_dims,
                                 k = 5, max.k = 300)
mix_integrated   <- MixingMetric(merged_seurat_int, grouping.var = "Batch",
                                 reduction = "integrated.rpca", dims = 1:n_dims,
                                 k = 5, max.k = 300)

mixing_summary <- data.frame(
  stage         = c("unintegrated", "integrated_rpca"),
  median_metric = c(median(mix_unintegrated), median(mix_integrated))
)
write.csv(mixing_summary, file.path(tab_dir, "SeuratMixingMetric_summary.csv"),
          row.names = FALSE)
print(mixing_summary)
