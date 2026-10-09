# 00_prepare_resources.R
# Builds the gene annotation files used by 05_differential_expression.R from
# Ensembl (biomaRt):
#   resources/xy_chromosome_genes.txt  genes on chromosomes X and Y
#   resources/gene_cytoband.csv        gene symbol and "GENE (cytoband)" label

suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
})

dir.create("resources", showWarnings = FALSE)

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
anno_gene <- getBM(attributes = c("external_gene_name", "chromosome_name", "band"),
                   mart = ensembl)

anno_gene <- anno_gene %>%
  filter(chromosome_name %in% c(1:22, "X", "Y"), external_gene_name != "")

xy_genes <- unique(anno_gene$external_gene_name[anno_gene$chromosome_name %in% c("X", "Y")])
writeLines(xy_genes, "resources/xy_chromosome_genes.txt")

gene_cytoband <- anno_gene %>%
  mutate(cytoband = paste0(chromosome_name, band),
         final_name = paste0(external_gene_name, " (", cytoband, ")")) %>%
  distinct(external_gene_name, .keep_all = TRUE) %>%
  select(external_gene_name, final_name, cytoband)
write.csv(gene_cytoband, "resources/gene_cytoband.csv", row.names = FALSE)
