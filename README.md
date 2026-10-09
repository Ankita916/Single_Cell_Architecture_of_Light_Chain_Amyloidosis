# Single_Cell_Architecture_of_Light_Chain_Amyloidosis

Code for quality control, RPCA integration, batch-mixing evaluation,
cell-cell communication and differential expression analysis of bone marrow
scRNA-seq samples from patients with AL amyloidosis.

## Requirements

- R (>= 4.0)
- Seurat v5, ggplot2, patchwork, FNN
- CellChat, openxlsx, NMF, ggalluvial (cell-cell communication analysis)
- infercnv, Polychrome
- survival, survminer (survival analysis); biomaRt (annotation files)
- DESeq2, SCPA, msigdbr, pheatmap, ComplexHeatmap, dplyr (differential expression and pathway analysis)

## Repository layout

```
scripts/00_prepare_resources.R  builds resources/ gene annotation files from Ensembl
config/survival.csv          survival status (1 = event) and time in days per sample
config/samples.csv           sample sheet (sample ID, batch, diagnosis, light-chain type)
scripts/01_qc_and_merge.R    per-sample QC, merging, unintegrated PCA/UMAP
scripts/02_integration.R     RPCA integration, clustering, UMAP
scripts/03_batch_mixing.R    UMAP by batch, LISI and Seurat MixingMetric
scripts/04_cellchat.R        cell-cell communication analysis with CellChat
scripts/05_differential_expression.R  plasma cell DE (single-cell and pseudobulk) and SCPA
scripts/06_survival.R        ligand-receptor gene and exhaustion-score survival analysis
scripts/07_infercnv_run.R    inferCNV per sample
scripts/08_infercnv_summary.R  inferCNV heatmaps, subgroup survival and Cox regression
scripts/09_plasma_cell_clonality.R  light-chain restriction of plasma cells and FISH-surrogate check
config/patient_annotation.csv  clinical and FISH annotations per sample (heatmap tracks)
resources/gene_order.txt     inferCNV gene order file (gene, chromosome, start, end)
config/cluster_annotation.csv  cluster-to-cell-type table (seurat_clusters, Final_label)
resources/xy_chromosome_genes.txt  sex-chromosome genes excluded from DE (one symbol per line)
resources/gene_cytoband.csv  heatmap gene labels (external_gene_name, final_name)
```

`04_cellchat.R` expects an annotated Seurat object at
`results/objects/annotated_seurat.rds`, with cell-type labels as the active
identities. `05_differential_expression.R` uses
`results/objects/merged_seurat_integrated.rds` and the three files above.

## Data

Raw data are not included. Place the 10x Genomics output of each sample in
`data/<sample_id>/` (matrix, features and barcodes files), where
`<sample_id>` matches the `sample_id` column in `config/samples.csv`.

## Usage

Run from the repository root:

```r
source("scripts/00_prepare_resources.R")  # once, needs internet
source("scripts/01_qc_and_merge.R")
source("scripts/02_integration.R")
source("scripts/03_batch_mixing.R")
source("scripts/04_cellchat.R")
source("scripts/05_differential_expression.R")
source("scripts/06_survival.R")
source("scripts/07_infercnv_run.R")
source("scripts/08_infercnv_summary.R")
source("scripts/09_plasma_cell_clonality.R")
```

Outputs are written to `results/objects`, `results/figures` and
`results/tables`.

## QC filters (per sample)

- Genes detected in at least 3 cells; cells with at least 200 genes
- log10(genes) / log10(UMIs) > 0.8
- Mitochondrial content below 3 x the sample median (maximum 20%)
- Number of detected genes below the 95th percentile of the sample
