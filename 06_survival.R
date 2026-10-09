# 06_survival.R
# Survival analysis of immune and plasma cell expression programmes:
#   (1) ligand-receptor genes: per-sample mean expression in a cell type,
#       univariate Cox regression and Kaplan-Meier curves (median split)
#   (2) T/NK cell exhaustion score: per-sample mean module score,
#       Kaplan-Meier curves (median and upper-quartile split) and Cox models
#
# Input:  results/objects/annotated_seurat.rds  (cell-type labels as identities,
#                                                sample ID in new.ident)
#         config/survival.csv  (sample_id, os_status [1 = event], overall_survival_days)
# Output: results/figures/survival/, results/tables/

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(survminer)
})

fig_dir <- "results/figures/survival"
tab_dir <- "results/tables"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Inputs and parameters -------------------------------------------------
seurat_obj  <- readRDS("results/objects/annotated_seurat.rds")
survival_df <- read.csv("config/survival.csv") %>%
  rename(time = overall_survival_days, status = os_status) %>%
  filter(!is.na(time), !is.na(status))

lr_genes <- unique(c(
  "APP", "BSG", "CD38", "CD47", "CD74", "COL9A2", "CXCR4", "LGALS9", "MDK", "MIF",
  "NCL", "P4HB", "PECAM1", "PPIA", "SDC1", "THBS1", "TNFRSF13B", "TNFRSF17",
  "TNFSF13B", "ADGRE5", "CD4", "CD44", "CD55", "CD8A", "CD8B", "CD94", "HLA-A",
  "HLA-B", "HLA-C", "HLA-DOB", "HLA-E", "HLA-F", "ICAM2", "ITGA4", "ITGAL",
  "ITGAM", "ITGB1", "ITGB2", "KLRC1", "LILRB1", "LRP1", "NKG2A", "PTGES3",
  "PTGER2", "PTGER4", "KLRD1", "PTGES"))
lr_cell_types <- c("Plasma", "NK")

exhaustion_markers <- list(
  NK         = c("KLRC1", "KLRD1", "TIGIT", "LAG3", "CD96", "ZNF683", "HAVCR2",
                 "TOX", "CD244", "CD160"),
  CD8_T_cell = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "TOX", "CD244",
                 "CD160", "ENTPD1", "EOMES", "ZNF683", "IKZF2"),
  CD4_T_cell = c("PDCD1", "LAG3", "TIGIT", "CTLA4", "TOX", "ENTPD1", "IKZF2",
                 "HAVCR2", "CD244"))

group_colors <- c(Low = "#1F78B4", High = "#E31A1C")

# ---- Helper functions ------------------------------------------------------
# Mean log-normalised expression per sample (rows: samples, columns: genes)
sample_means <- function(cells, genes) {
  expr <- GetAssayData(cells, assay = "RNA", layer = "data")[genes, , drop = FALSE]
  by_sample <- split(seq_len(ncol(expr)), cells$new.ident)
  means <- sapply(by_sample, function(i) Matrix::rowMeans(expr[, i, drop = FALSE]))
  data.frame(sample_id = colnames(means), t(means), check.names = FALSE)
}

# High / Low groups at a quantile cutoff
make_groups <- function(x, prob = 0.5, inclusive = FALSE) {
  cutoff <- quantile(x, prob, na.rm = TRUE)
  high <- if (inclusive) x >= cutoff else x > cutoff
  factor(ifelse(high, "High", "Low"), levels = c("Low", "High"))
}

# Hazard ratio, 95% CI and p-values of a single-term Cox model
cox_summary <- function(fit) {
  s <- summary(fit)
  data.frame(hazard_ratio = s$conf.int[1, "exp(coef)"],
             lower_CI     = s$conf.int[1, "lower .95"],
             upper_CI     = s$conf.int[1, "upper .95"],
             p_value      = s$coefficients[1, "Pr(>|z|)"],
             p_logrank_test = s$logtest[["pvalue"]])
}

# Univariate Cox regression per gene (expression scaled to SD units)
run_cox_genes <- function(df, genes) {
  out <- lapply(genes, function(g) {
    if (sd(df[[g]]) == 0) return(NULL)
    df$z <- as.numeric(scale(df[[g]]))
    cbind(gene = g, cox_summary(coxph(Surv(time, status) ~ z, data = df)))
  })
  res <- do.call(rbind, out)
  res$significant <- res$p_value < 0.05
  res
}

plot_hazard <- function(res, filename, width = 2000, height = 1800) {
  p <- ggplot(res, aes(reorder(gene, hazard_ratio), hazard_ratio, color = significant)) +
    geom_point(size = 3, shape = 15) +
    geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.8, linewidth = 1.2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_text(aes(y = upper_CI, label = sprintf("p=%.3g", p_value)),
              hjust = -0.2, size = 3) +
    coord_flip() +
    scale_color_manual(values = c("black", "red")) +
    labs(title = "Univariate Cox hazard ratios with 95% CI", x = "Gene",
         y = "Hazard ratio (per SD increase)") +
    theme_bw()
  ggsave(filename, p, width = width, height = height, units = "px", dpi = 300)
}

plot_km <- function(df, filename, title) {
  fit <- survfit(Surv(time, status) ~ group, data = df)
  p <- ggsurvplot(fit, data = df, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
                  legend.labs = levels(df$group), palette = unname(group_colors),
                  title = title)
  tiff(filename, width = 2000, height = 2000, res = 300)
  print(p)
  dev.off()
}

# ---- (1) Ligand-receptor genes ---------------------------------------------
for (ct in lr_cell_types) {
  cells <- subset(seurat_obj, idents = ct)
  genes <- intersect(lr_genes, rownames(cells))
  message(ct, ": genes not found: ", paste(setdiff(lr_genes, genes), collapse = ", "))

  df  <- merge(survival_df, sample_means(cells, genes), by = "sample_id")
  out <- file.path(fig_dir, ct)
  dir.create(out, showWarnings = FALSE)

  cox_res <- run_cox_genes(df, genes)
  write.csv(cox_res, file.path(tab_dir, paste0("cox_LR_genes_", ct, ".csv")),
            row.names = FALSE)
  plot_hazard(cox_res, file.path(out, "LR_hazard_ratios.tiff"))

  for (g in genes) {
    df$group <- make_groups(df[[g]])
    plot_km(df, file.path(out, paste0(g, "_KM.tiff")), paste("Survival by", g, "in", ct))
  }
}

# ---- (2) Exhaustion score --------------------------------------------------
pastel <- c("#FADADD", "#FACBAA", "#FAF1D6", "#DFF7E1", "#B4E7D9", "#B1D6E6",
            "#C0B9E2", "#E2C3E8", "#F2CDD2", "#F9B7AC", "#FFDAC1", "#E2F0CB",
            "#BDE0D0", "#A9DEF9", "#CBBEE7", "#E4C1F9", "#FBE4FF", "#FCCDE2",
            "#F8B195", "#F67280")
cutoffs <- c(median = 0.5, upper_quartile = 0.75)

cox_scores <- list()
for (ct in names(exhaustion_markers)) {
  cells <- subset(seurat_obj, idents = ct)
  cells <- AddModuleScore(cells, features = list(exhaustion_markers[[ct]]),
                          name = "Exhaustion")

  p_vln <- VlnPlot(cells, features = "Exhaustion1", group.by = "new.ident",
                   pt.size = 0.1, cols = rep_len(pastel, n_distinct(cells$new.ident))) +
    geom_jitter(size = 0.4, width = 0.2, alpha = 0.6) +
    labs(title = paste(ct, "exhaustion score"), y = "Exhaustion score") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  ggsave(file.path(fig_dir, paste0("Exhaustion_score_", ct, ".png")), p_vln,
         width = 2000, height = 1500, units = "px", dpi = 300)

  score_df <- cells@meta.data %>%
    group_by(sample_id = new.ident) %>%
    summarise(score = mean(Exhaustion1), .groups = "drop") %>%
    inner_join(survival_df, by = "sample_id")

  for (cut in names(cutoffs)) {
    score_df$group <- make_groups(score_df$score, cutoffs[[cut]], inclusive = TRUE)
    plot_km(score_df, file.path(fig_dir, paste0("Exhaustion_KM_", ct, "_", cut, ".tiff")),
            paste("Survival by", ct, "exhaustion score"))
    cox_scores[[paste(ct, cut)]] <- cbind(
      cell_type = ct, model = paste("High vs Low,", cut),
      cox_summary(coxph(Surv(time, status) ~ group, data = score_df)))
  }
  cox_scores[[paste(ct, "continuous")]] <- cbind(
    cell_type = ct, model = "Continuous score",
    cox_summary(coxph(Surv(time, status) ~ score, data = score_df)))
}

write.csv(do.call(rbind, cox_scores), file.path(tab_dir, "cox_exhaustion_score.csv"),
          row.names = FALSE)
