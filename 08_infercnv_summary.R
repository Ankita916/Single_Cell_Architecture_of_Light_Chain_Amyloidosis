# 08_infercnv_summary.R
# Summarises inferCNV results for malignant plasma cells across samples:
#   (1) mean inferCNV signal per gene and per cytoband in each sample
#   (2) heatmaps with clinical annotations
#   (3) Kaplan-Meier comparison of two sample subgroups
#   (4) Cox regression of overall survival on the copy-number signal
#
# Input:  results/infercnv/<sample_id>/run.final.infercnv_obj  (from 07)
#         resources/gene_order.txt, resources/gene_cytoband.csv
#         config/survival.csv  (sample_id, os_status [1 = event], overall_survival_days)
#         config/patient_annotation.csv  (sample_id plus clinical and FISH columns;
#                                         0/1 columns are drawn as binary tracks)
# Output: results/figures/infercnv/, results/tables/

suppressPackageStartupMessages({
  library(dplyr)
  library(infercnv)
  library(pheatmap)
  library(Polychrome)
  library(survival)
  library(survminer)
})

infercnv_dir <- "results/infercnv"
fig_dir      <- "results/figures/infercnv"
tab_dir      <- "results/tables"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

malignant_cell_type <- "Plasma"
subgroup_1 <- c("P01", "P17", "P02", "P06", "P03", "P04", "P10", "P11", "P18", "P13_D")
subgroup_2 <- c("P05", "P07", "P08", "P09", "P12", "P14", "P15", "P16", "P19")
cox_level  <- "cytoband"                        # "cytoband" or "gene"
excluded_from_cox <- "P13_R"                    # second sample of the same patient

# ---- Mean inferCNV signal of malignant plasma cells per sample -------------
sample_files <- list.files(infercnv_dir, pattern = "run.final.infercnv_obj",
                           recursive = TRUE, full.names = TRUE)

mean_signal <- lapply(sample_files, function(f) {
  obj   <- readRDS(f)
  cells <- obj@observation_grouped_cell_indices[[malignant_cell_type]]
  rowMeans(obj@expr.data[, cells, drop = FALSE], na.rm = TRUE)
})
names(mean_signal) <- basename(dirname(sample_files))

# Samples (rows) x genes in genomic order
gene_order <- read.table("resources/gene_order.txt", sep = "\t", header = FALSE,
                         col.names = c("gene", "chromosome", "start", "end"))
all_genes  <- Reduce(union, lapply(mean_signal, names))
gene_mat   <- do.call(rbind, lapply(mean_signal, function(v) v[all_genes]))
genes      <- gene_order$gene[gene_order$gene %in% colnames(gene_mat)]
gene_mat   <- gene_mat[sort(rownames(gene_mat)), genes]

# Mean signal per cytoband (autosomes), in genomic order
gene_cytoband <- read.csv("resources/gene_cytoband.csv", stringsAsFactors = FALSE) %>%
  filter(grepl("^[0-9]+[pq]", cytoband)) %>%
  distinct(external_gene_name, .keep_all = TRUE)
gene_map <- data.frame(gene = genes) %>%
  inner_join(gene_cytoband, by = c("gene" = "external_gene_name"))
cytoband_order <- unique(gene_map$cytoband)
cytoband_mat <- sapply(cytoband_order, function(cb) {
  rowMeans(gene_mat[, gene_map$gene[gene_map$cytoband == cb], drop = FALSE],
           na.rm = TRUE)
})

# ---- Annotations -----------------------------------------------------------
patient_annotation <- read.csv("config/patient_annotation.csv", check.names = FALSE,
                               stringsAsFactors = FALSE)
rownames(patient_annotation) <- patient_annotation$sample_id
patient_annotation$sample_id <- NULL
patient_annotation <- patient_annotation[rownames(gene_mat), , drop = FALSE]

is_binary <- sapply(patient_annotation, function(x) all(na.omit(x) %in% c(0, 1)))
patient_annotation[is_binary] <- lapply(patient_annotation[is_binary], as.character)

clinical_colors <- c(
  list(Mayo_2004_Staging = c(I = "lightseagreen", II = "#FFD700", III = "indianred"),
       Light_chain_type = c(kappa = "aquamarine2", lambda = "mediumpurple3"),
       Disease_status = c(Diagnosis = "#1B9E77", Relapse = "#E69F00"),
       Gender = c(Male = "#1F77B4", Female = "#DB7093"),
       Diagnosis_sub_category = c(AL_MGUS = "lightblue", AL_MM = "lightcoral",
                                  AL_SMM = "#D8BFD8")),
  setNames(rep(list(c("0" = "white", "1" = "black")), sum(is_binary)),
           names(patient_annotation)[is_binary]))

glasbey <- glasbey.colors(32)
chromosome_colors <- function(chromosomes, drop) {
  chromosomes <- unique(chromosomes)
  setNames(glasbey[-drop][seq_along(chromosomes)], chromosomes)
}

# ---- Heatmaps --------------------------------------------------------------
heatmap_colors <- colorRampPalette(c("darkblue", "white", "darkred"))(25)
heatmap_colors[13] <- "white"

plot_cnv_heatmap <- function(mat, filename, col_annotation, chromosome_drop,
                             breaks, cluster_rows = FALSE, gaps_col = NULL,
                             width = 3500, height = 2000) {
  annotation_colors <- c(
    list(chromosome = chromosome_colors(col_annotation$chromosome, chromosome_drop)),
    clinical_colors)
  tiff(filename, width = width, height = height, res = 300)
  pheatmap(mat, cluster_rows = cluster_rows, cluster_cols = FALSE,
           show_colnames = FALSE, show_rownames = TRUE,
           annotation_col = col_annotation, annotation_row = patient_annotation,
           annotation_colors = annotation_colors, color = heatmap_colors,
           breaks = breaks, na_col = "grey", gaps_col = gaps_col)
  dev.off()
}

# Per gene
gene_annotation <- data.frame(chromosome = gene_order$chromosome[match(genes, gene_order$gene)],
                              row.names = genes)
plot_cnv_heatmap(gene_mat, file.path(fig_dir, "InferCNV_heatmap_Plasma_cells_gene.tiff"),
                 gene_annotation, chromosome_drop = 1,
                 breaks = seq(0.9, 1.1, length.out = 25))

# Per cytoband
cytoband_annotation <- data.frame(
  chromosome = paste0("chr", sub("[pq].*", "", cytoband_order)),
  row.names = cytoband_order)
plot_cnv_heatmap(cytoband_mat,
                 file.path(fig_dir, "InferCNV_heatmap_Plasma_cells_cytoband.tiff"),
                 cytoband_annotation, chromosome_drop = c(1, 5, 16, 19),
                 breaks = seq(0.94, 1.06, length.out = 25))

# Per cytoband, samples clustered, chromosomes separated
plot_cnv_heatmap(cytoband_mat,
                 file.path(fig_dir, "InferCNV_heatmap_Plasma_cells_cytoband_clustered.tiff"),
                 cytoband_annotation, chromosome_drop = c(1, 5, 16, 19),
                 breaks = seq(0.94, 1.06, length.out = 25), cluster_rows = TRUE,
                 gaps_col = cumsum(rle(cytoband_annotation$chromosome)$lengths),
                 width = 5000, height = 4000)

# ---- Survival --------------------------------------------------------------
survival_df <- read.csv("config/survival.csv") %>%
  rename(time = overall_survival_days, status = os_status) %>%
  filter(!is.na(time), !is.na(status), time > 0)

# Kaplan-Meier: subgroup 1 versus subgroup 2
km_df <- survival_df %>%
  mutate(group = factor(case_when(sample_id %in% subgroup_1 ~ "Subgroup1",
                                  sample_id %in% subgroup_2 ~ "Subgroup2"),
                        levels = c("Subgroup1", "Subgroup2"))) %>%
  filter(!is.na(group))
group_n <- table(km_df$group)

km_fit <- survfit(Surv(time, status) ~ group, data = km_df)
km_plot <- ggsurvplot(km_fit, data = km_df, pval = TRUE, pval.method = TRUE,
                      conf.int = TRUE, risk.table = TRUE, risk.table.col = "group",
                      risk.table.height = 0.25, tables.theme = theme_cleantable(),
                      legend.title = "Patient group",
                      legend.labs = paste0(names(group_n), " (n=", group_n, ")"),
                      palette = c("#E7B800", "#2E9FDF"), ggtheme = theme_bw())
tiff(file.path(fig_dir, "KM_inferCNV_subgroups.tiff"), width = 2000, height = 2000, res = 300)
print(km_plot)
dev.off()

# Cox regression of overall survival on the inferCNV signal
cox_mat <- if (cox_level == "cytoband") cytoband_mat else gene_mat
cox_mat <- cox_mat[!rownames(cox_mat) %in% excluded_from_cox, , drop = FALSE]
cox_df  <- inner_join(survival_df, data.frame(sample_id = rownames(cox_mat), cox_mat,
                                              check.names = FALSE),
                      by = "sample_id")

cox_results <- lapply(colnames(cox_mat), function(feature) {
  x <- cox_df[[feature]]
  if (sum(!is.na(x)) < 5) return(NULL)
  fit <- tryCatch(coxph(Surv(cox_df$time, cox_df$status) ~ x), error = function(e) NULL)
  if (is.null(fit) || fit$nevent == 0) return(NULL)
  s <- summary(fit)
  data.frame(feature = feature, n_obs = s$n,
             hazard_ratio = s$coefficients[1, "exp(coef)"],
             lower_CI = s$conf.int[1, "lower .95"],
             upper_CI = s$conf.int[1, "upper .95"],
             p_value = s$coefficients[1, "Pr(>|z|)"])
})
cox_results <- do.call(rbind, cox_results)
cox_results$padj <- p.adjust(cox_results$p_value, method = "fdr")
write.csv(cox_results[order(cox_results$p_value), ],
          file.path(tab_dir, paste0("infercnv_cox_", cox_level, ".csv")),
          row.names = FALSE)
