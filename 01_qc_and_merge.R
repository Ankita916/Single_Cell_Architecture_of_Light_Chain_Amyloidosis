# 01_qc_and_merge.R
# Per-sample quality control, merging and unintegrated embedding of
# bone marrow scRNA-seq samples.
#
# Input:  config/samples.csv, data/<sample_id>/ (10x Genomics output)
# Output: results/objects/merged_seurat.rds

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(1234)

# Parameters
sample_sheet <- "config/samples.csv"
data_root    <- "data"
out_dir      <- "results/objects"
n_dims       <- 21
cluster_res  <- 0.8

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Quality control of a single sample
qc_sample <- function(sample_id, batch, diagnosis, lc_type) {
  counts <- Read10X(data.dir = file.path(data_root, sample_id))
  obj <- CreateSeuratObject(counts = counts, project = sample_id,
                            min.cells = 3, min.features = 200)

  # QC metrics
  obj[["percent.mt"]]       <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent_ribo"]]     <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  obj[["percent_hb"]]       <- PercentageFeatureSet(obj, pattern = "^HB[^(P)]")
  obj[["percent.ig"]]       <- PercentageFeatureSet(obj, pattern = "^IG[KLH]")
  obj$log10GenesPerUMI      <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)

  # Remove cells with low complexity or too few detected genes
  obj <- subset(obj, subset = log10GenesPerUMI > 0.8 & nFeature_RNA > 200)

  # Remove cells with high mitochondrial content
  # (cutoff: 3 x sample median, capped at 20%)
  mt_cutoff <- min(3 * median(obj$percent.mt), 20)
  obj <- subset(obj, cells = colnames(obj)[obj$percent.mt < mt_cutoff])

  # Remove cells in the top 5% of detected genes (potential multiplets)
  feature_cutoff <- sort(obj$nFeature_RNA)[round(ncol(obj) * 0.95)]
  obj <- subset(obj, cells = colnames(obj)[obj$nFeature_RNA < feature_cutoff])

  # Sample metadata
  obj$Batch     <- batch
  obj$new.ident <- sample_id
  obj$Diagnosis <- diagnosis
  obj$lc_type   <- lc_type
  obj
}

samples <- read.csv(sample_sheet, stringsAsFactors = FALSE)

seurat_list <- lapply(seq_len(nrow(samples)), function(i) {
  qc_sample(samples$sample_id[i], samples$batch[i],
            samples$diagnosis[i], samples$lc_type[i])
})
names(seurat_list) <- samples$sample_id

# Merge samples
merged_seurat <- merge(seurat_list[[1]], y = seurat_list[-1],
                       add.cell.ids = samples$sample_id,
                       project = "Amyloidosis")

# Normalisation, dimensionality reduction and clustering (unintegrated)
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:n_dims)
merged_seurat <- FindClusters(merged_seurat, resolution = cluster_res)
merged_seurat <- RunUMAP(merged_seurat, dims = 1:n_dims)

# Join layers, then split the RNA assay by sample for integration
merged_seurat[["RNA"]] <- JoinLayers(merged_seurat[["RNA"]])
merged_seurat[["RNA"]] <- split(merged_seurat[["RNA"]],
                                f = merged_seurat$orig.ident)

saveRDS(merged_seurat, file.path(out_dir, "merged_seurat.rds"))
