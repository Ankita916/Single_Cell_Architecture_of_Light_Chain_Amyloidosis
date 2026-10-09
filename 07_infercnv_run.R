# 07_infercnv_run.R
# Runs inferCNV on each sample, using non-malignant cell types as reference.
#
# Input:  results/objects/annotated_seurat.rds  (cell-type labels as identities,
#                                                sample ID in new.ident)
#         resources/gene_order.txt  (no header; tab-separated: gene, chromosome,
#                                    start, end)
# Output: results/infercnv/<sample_id>/

suppressPackageStartupMessages({
  library(Seurat)
  library(infercnv)
})

out_root        <- "results/infercnv"
gene_order_file <- "resources/gene_order.txt"
reference_cell_types <- c("CD8_T_cell", "NK", "CD4_T_cell", "B_cell", "CD14_Mono",
                          "CD16_Mono", "Late_Eryth", "GMP", "HSC", "Prolif_cells",
                          "pDC", "cDC2")

dir.create(file.path(out_root, "annotations"), recursive = TRUE, showWarnings = FALSE)

seurat_obj <- readRDS("results/objects/annotated_seurat.rds")
seurat_obj$cell_type <- as.character(Idents(seurat_obj))

for (sample_id in unique(seurat_obj$new.ident)) {
  out_dir <- file.path(out_root, sample_id)
  if (file.exists(file.path(out_dir, "run.final.infercnv_obj"))) next

  sample_obj <- seurat_obj[, seurat_obj$new.ident == sample_id]

  annotation <- data.frame(cell = colnames(sample_obj), cell_type = sample_obj$cell_type)
  annotation_file <- file.path(out_root, "annotations",
                               paste0("annotation_", sample_id, ".txt"))
  write.table(annotation, annotation_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)

  # Reference groups present in this sample
  ref_groups <- intersect(reference_cell_types, annotation$cell_type)

  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = LayerData(sample_obj, assay = "RNA", layer = "counts"),
    annotations_file  = annotation_file,
    gene_order_file   = gene_order_file,
    ref_group_names   = ref_groups)

  infercnv::run(infercnv_obj,
                cutoff = 0.1,   # recommended for 10x Genomics data
                out_dir = out_dir,
                cluster_by_groups = TRUE,
                denoise = TRUE,
                HMM = TRUE,
                HMM_type = "i3")
}
