# 05_differential_expression.R
# Differential gene expression and pathway analysis of plasma cells,
# AL_MM versus AL_MGUS:
#   (1) single-cell test (Wilcoxon) and pseudobulk test (DESeq2)
#   (2) genes significant in both tests, heatmaps of average expression
#   (3) pathway analysis with SCPA (MSigDB Hallmark gene sets)
#
# Input:  results/objects/merged_seurat_integrated.rds
#         config/cluster_annotation.csv   (columns: seurat_clusters, Final_label)
#         resources/xy_chromosome_genes.txt  (one gene symbol per line)
#         resources/gene_cytoband.csv     (columns: external_gene_name, final_name)
# Output: results/tables/, results/figures/differential_expression/

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(openxlsx)
  library(SCPA)
  library(msigdbr)
  library(ComplexHeatmap)
})

# Parameters
obj_dir     <- "results/objects"
fig_dir     <- "results/figures/differential_expression"
tab_dir     <- "results/tables"
p_cutoff    <- 0.05   # unadjusted p-value cutoff, applied to both tests
lfc_cutoff  <- 0.5      # absolute log2 fold-change cutoff, applied to both tests
group_1     <- "Plasma_AL_MM"
group_2     <- "Plasma_AL_MGUS"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Data preparation ------------------------------------------------------
seurat_full <- readRDS(file.path(obj_dir, "merged_seurat_integrated.rds"))

# Cell-type labels from the cluster annotation table
cluster_annotation <- read.csv("config/cluster_annotation.csv",
                               stringsAsFactors = FALSE)
label_lookup <- setNames(cluster_annotation$Final_label,
                         cluster_annotation$seurat_clusters)
seurat_full$customclassif <- unname(
  label_lookup[as.character(seurat_full$seurat_clusters)])

seurat_full$celltype.Diagnosis <- paste(seurat_full$customclassif,
                                        seurat_full$Diagnosis, sep = "_")
seurat_full$donor_id.Diagnosis <- paste0(seurat_full$Diagnosis, "-",
                                         seurat_full$new.ident)
Idents(seurat_full) <- "celltype.Diagnosis"


# Cytoband annotation of gene names for heatmap labels
gene_cytoband <- read.csv("resources/gene_cytoband.csv", stringsAsFactors = FALSE)
to_cytoband_name <- function(genes) {
  idx <- match(genes, gene_cytoband$external_gene_name)
  genes[!is.na(idx)] <- gene_cytoband$final_name[idx[!is.na(idx)]]
  genes
}

# ---- Differential expression -----------------------------------------------
# Single-cell level
de_sc <- FindMarkers(seurat_obj, ident.1 = group_1, ident.2 = group_2,
                     verbose = FALSE)

# Pseudobulk level: counts aggregated per diagnosis, sample and cell type
# (underscores in group names are converted to dashes by AggregateExpression)
pseudo_obj <- AggregateExpression(seurat_obj, assays = "RNA", return.seurat = TRUE,
                                  group.by = c("Diagnosis", "new.ident",
                                               "customclassif"))

pseudo_meta <- do.call(rbind, strsplit(Cells(pseudo_obj), "_"))
colnames(pseudo_meta) <- c("Diagnosis", "new_ident", "customclassif")
rownames(pseudo_meta) <- Cells(pseudo_obj)
pseudo_obj <- AddMetaData(pseudo_obj, as.data.frame(pseudo_meta))
pseudo_obj$celltype.Diagnosis <- paste(pseudo_obj$customclassif,
                                       pseudo_obj$Diagnosis, sep = "_")
Idents(pseudo_obj) <- "celltype.Diagnosis"

de_bulk <- FindMarkers(pseudo_obj,
                       ident.1 = gsub("_AL_", "_AL-", group_1),
                       ident.2 = gsub("_AL_", "_AL-", group_2),
                       test.use = "DESeq2")

# Combine the two sets of results
names(de_bulk) <- paste0(names(de_bulk), ".bulk")
de_bulk$gene   <- rownames(de_bulk)
names(de_sc)   <- paste0(names(de_sc), ".sc")
de_sc$gene     <- rownames(de_sc)

merge_dat <- merge(de_sc, de_bulk, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val.bulk), ]

# Genes significant in both tests
common_idx <- which(merge_dat$p_val.sc < p_cutoff & merge_dat$p_val.bulk < p_cutoff)
common     <- merge_dat$gene[common_idx]
message("Genes significant in both tests: ", length(common))

# Up- and downregulated in AL_MM with consistent large fold changes
up_df <- merge_dat[which(merge_dat$p_val.sc < p_cutoff &
                           merge_dat$p_val.bulk < p_cutoff &
                           merge_dat$avg_log2FC.sc > lfc_cutoff &
                           merge_dat$avg_log2FC.bulk > lfc_cutoff), ]
down_df <- merge_dat[which(merge_dat$p_val.sc < p_cutoff &
                             merge_dat$p_val.bulk < p_cutoff &
                             merge_dat$avg_log2FC.sc < -lfc_cutoff &
                             merge_dat$avg_log2FC.bulk < -lfc_cutoff), ]

up_df$direction   <- "UP_AL-MM"
down_df$direction <- "DOWN_AL-MM"
write.xlsx(rbind(up_df, down_df),
           file.path(tab_dir, "UP_Down_genes_MM_vs_MGUS_FC_cutoff_2.xlsx"))

# ---- Heatmaps of average expression per sample -----------------------------
plasma_de <- subset(seurat_obj, idents = c(group_1, group_2))

avg_expr <- AverageExpression(plasma_de, features = common,
                              group.by = "donor_id.Diagnosis",
                              layer = "data")$RNA
scaled_expr <- t(scale(t(as.matrix(avg_expr))))
rownames(scaled_expr) <- to_cytoband_name(rownames(scaled_expr))

heatmap_colors <- colorRampPalette(c("darkblue", "white", "darkred"))(100)
heatmap_breaks <- seq(-2, 2, length.out = 101)

# All genes significant in both tests (samples clustered)
png(file.path(fig_dir, "Average_expression_plasma_common_genes_MM_MGUS.png"),
    width = 2500, height = 1000, res = 300)
pheatmap(t(scaled_expr), cluster_cols = TRUE, cluster_rows = FALSE,
         cellwidth = 2.5, cellheight = 3, treeheight_row = 100,
         color = heatmap_colors, breaks = heatmap_breaks, fontsize = 2.5,
         legend = TRUE)
dev.off()

# Genes passing the fold-change cutoff
sig_genes  <- to_cytoband_name(union(up_df$gene, down_df$gene))
subset_mat <- scaled_expr[rownames(scaled_expr) %in% sig_genes, ]

# Sample and diagnosis labels for the columns (names use dashes after averaging)
sample_info <- unique(data.frame(
  donor     = gsub("_", "-", plasma_de$donor_id.Diagnosis),
  sample_id = plasma_de$new.ident,
  diagnosis = gsub("_", "-", plasma_de$Diagnosis)
))
sample_info <- sample_info[match(colnames(subset_mat), sample_info$donor), ]

colnames(subset_mat) <- sample_info$sample_id
annotation_col <- data.frame(diagnosis = sample_info$diagnosis,
                             row.names = sample_info$sample_id)
annotation_colors <- list(diagnosis = c("AL-MGUS" = "lightblue",
                                        "AL-MM"   = "lightcoral"))

tiff(file.path(fig_dir, "Average_expression_plasma_up_down_genes_MM_vs_MGUS.tiff"),
     width = 1500, height = 2500, res = 300)
pheatmap(subset_mat, cluster_cols = FALSE, cluster_rows = TRUE,
         cellwidth = 5, cellheight = 5, treeheight_row = 35, treeheight_col = 35,
         color = heatmap_colors, breaks = heatmap_breaks, legend = TRUE,
         annotation_col = annotation_col, annotation_colors = annotation_colors,
         fontsize = 5)
dev.off()

# ---- Pathway analysis with SCPA (all genes) --------------------------------
plasma_scpa <- subset(seurat_full, idents = c(group_1, group_2))

cells_mm <- seurat_extract(plasma_scpa, meta1 = "celltype.Diagnosis",
                           value_meta1 = group_1)
cells_mgus <- seurat_extract(plasma_scpa, meta1 = "celltype.Diagnosis",
                             value_meta1 = group_2)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  format_pathways()

scpa_out <- compare_pathways(samples = list(cells_mm, cells_mgus),
                             pathways = hallmark)
scpa_out$Pathway <- gsub("HALLMARK_", "", scpa_out$Pathway)
write.xlsx(scpa_out, file.path(tab_dir, "scpa_out_MM_MGUS.xlsx"))

# Heatmap of pathway q-values
tiff(file.path(fig_dir, "scpa_MM_vs_MGUS_pathway_Qval.tiff"),
     width = 1500, height = 2000, res = 300)
ht <- withVisible(plot_heatmap(scpa_out,
                               column_names = "PlasmaCells AL_MM vs AL_MGUS",
                               show_row_names = TRUE, row_fontsize = 7,
                               column_fontsize = 2))
if (ht$visible) print(ht$value)
dev.off()

# Pathway fold changes
p_fc <- ggplot(scpa_out, aes(reorder(Pathway, FC), FC)) +
  geom_col(aes(fill = adjPval < 0.05)) +
  coord_flip() +
  labs(x = "Pathway", y = "FC", title = "Hallmark pathways FC from SCPA") +
  theme_bw() +
  theme(axis.text  = element_text(size = 8, color = "black", face = "bold"),
        axis.title = element_text(size = 15, color = "black", face = "bold"))

ggsave(file.path(fig_dir, "scpa_MM_vs_MGUS_pathway_FC.tiff"), p_fc,
       width = 2000, height = 2000, units = "px", dpi = 300)
