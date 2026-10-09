# 09_plasma_cell_clonality.R
# Classifies plasma cells as malignant, residual normal or undetermined by
# immunoglobulin light-chain restriction, and cross-checks the calls with
# FISH-surrogate gene expression (CCND1, MAFB, CCND3).
#
# Each plasma cell expresses one light-chain class (kappa or lambda). A cell whose
# dominant class matches the clinical light-chain type of the patient is called
# malignant; a cell with the opposite class is called residual normal.
#
# Input:  results/objects/annotated_seurat.rds  (cell-type labels as identities,
#                                                sample ID in new.ident)
#         config/samples.csv  (sample_id, lc_type)
# Output: results/figures/clonality/, results/tables/

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

fig_dir <- "results/figures/clonality"
tab_dir <- "results/tables"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# Genes whose expression stands in for a FISH-detected translocation
translocation_surrogates <- c(P03 = "MAFB",    # t(14;20)
                              P04 = "CCND1",   # t(11;14)
                              P06 = "CCND1",
                              P07 = "CCND1",
                              P09 = "CCND1",
                              P16 = "CCND1",
                              P18 = "CCND3")   # t(6;14)

clonal_colors <- c(malignant = "firebrick", residual_normal = "steelblue",
                   undetermined = "grey70")

# ---- Plasma cells and light-chain genes ------------------------------------
seurat_obj <- readRDS("results/objects/annotated_seurat.rds")
samples    <- read.csv("config/samples.csv")
plasma_obj <- subset(seurat_obj, idents = "Plasma")

igk_genes <- grep("^IGKC$", rownames(plasma_obj), value = TRUE)
igl_genes <- grep("^IGLC[0-9]+$", rownames(plasma_obj), value = TRUE)
message("Lambda constant-region genes: ", paste(igl_genes, collapse = ", "))

# Lambda gene expression per sample
p_igl <- VlnPlot(plasma_obj, features = igl_genes, group.by = "new.ident",
                 pt.size = 0.1, ncol = 3)
ggsave(file.path(fig_dir, "VlnPlot_IGLC_genes_per_sample.png"), p_igl,
       width = 16, height = 10, dpi = 300)

igl_summary <- FetchData(plasma_obj, vars = c(igl_genes, "new.ident")) %>%
  group_by(new.ident) %>%
  summarise(across(all_of(igl_genes), list(mean = mean, pct_pos = ~ mean(. > 0) * 100)))
write.csv(igl_summary, file.path(tab_dir, "IGLC_expression_per_sample.csv"),
          row.names = FALSE)

# ---- Light-chain restriction -----------------------------------------------
expr <- GetAssayData(plasma_obj, assay = "RNA", layer = "data")
igk_sum <- Matrix::colSums(expr[igk_genes, , drop = FALSE])
igl_sum <- Matrix::colSums(expr[igl_genes, , drop = FALSE])

plasma_lc <- data.frame(
  cell      = colnames(plasma_obj),
  sample_id = plasma_obj$new.ident,
  clinical_lc = samples$lc_type[match(plasma_obj$new.ident, samples$sample_id)],
  dominant_lc = case_when(igk_sum > igl_sum ~ "kappa",
                          igl_sum > igk_sum ~ "lambda",
                          TRUE ~ "undetermined")) %>%
  mutate(clonal_call = case_when(dominant_lc == "undetermined" ~ "undetermined",
                                 dominant_lc == clinical_lc ~ "malignant",
                                 TRUE ~ "residual_normal"))
print(table(plasma_lc$sample_id, plasma_lc$clonal_call))

# Percentages among cells with a defined light chain
sample_summary <- plasma_lc %>%
  filter(dominant_lc != "undetermined") %>%
  group_by(sample_id, clinical_lc) %>%
  summarise(n_cells = n(),
            pct_malignant = round(100 * mean(clonal_call == "malignant"), 1),
            pct_residual_normal = round(100 * mean(clonal_call == "residual_normal"), 1),
            .groups = "drop")
write.csv(sample_summary, file.path(tab_dir, "PlasmaCell_ClonalityClassification.csv"),
          row.names = FALSE)

# Bar plot, all plasma cells per sample
plot_data <- plasma_lc %>%
  count(sample_id, clonal_call) %>%
  group_by(sample_id) %>%
  mutate(pct = 100 * n / sum(n))

p_bar <- ggplot(plot_data, aes(sample_id, pct, fill = clonal_call)) +
  geom_col() +
  scale_fill_manual(values = clonal_colors) +
  labs(x = "Sample", y = "% of plasma cells", fill = "Clonal classification") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "PlasmaCell_LightChainRestriction_per_sample.png"), p_bar,
       width = 12, height = 6, dpi = 300)

# ---- Cross-check with FISH-surrogate genes ---------------------------------
cross_validation <- lapply(names(translocation_surrogates), function(sample_id) {
  gene  <- translocation_surrogates[[sample_id]]
  calls <- plasma_lc %>%
    filter(sample_id == !!sample_id, clonal_call %in% c("malignant", "residual_normal"))
  if (!gene %in% rownames(expr) || nrow(calls) == 0) return(NULL)

  gene_expr <- expr[gene, calls$cell]
  is_malignant <- calls$clonal_call == "malignant"
  data.frame(sample_id = sample_id, surrogate_gene = gene,
             mean_expr_malignant = mean(gene_expr[is_malignant]),
             mean_expr_residual  = mean(gene_expr[!is_malignant]),
             n_malignant = sum(is_malignant),
             n_residual  = sum(!is_malignant),
             n_concordant_malignant = sum(is_malignant & gene_expr > 0))
})
cross_validation <- do.call(rbind, cross_validation)
write.csv(cross_validation, file.path(tab_dir, "FISH_surrogate_crosscheck.csv"),
          row.names = FALSE)

# Share of malignant-called cells that express the surrogate gene
message("Overall concordance with FISH surrogates: ",
        round(100 * sum(cross_validation$n_concordant_malignant) /
                sum(cross_validation$n_malignant), 1), "%")
