# 04_cellchat.R
# Cell-cell communication analysis with CellChat on the annotated bone marrow
# object: network inference, pathway visualisation, signalling roles,
# communication patterns and summary tables.
#
# Input:  results/objects/annotated_seurat.rds
#         (Seurat object with cell-type labels as the active identities)
# Output: results/objects/cellchat.rds, results/figures/cellchat/,
#         results/tables/

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(openxlsx)
  library(NMF)
  library(ggalluvial)
})

options(future.globals.maxSize = Inf)

# Parameters
obj_dir             <- "results/objects"
fig_dir             <- "results/figures/cellchat"
tab_dir             <- "results/tables"
min_cells           <- 10   # minimum cells per group for a communication to be kept
n_patterns_outgoing <- 3
n_patterns_incoming <- 4
vertex_receiver     <- c(1, 2, 3, 4, 8)   # groups shown on the left of hierarchy plots

pathways_focus <- c("MHC-I", "MHC-II", "CCL", "CXCL", "BAFF", "COLLAGEN",
                    "THBS", "ICAM", "PECAM1", "IFN-II", "GRN", "NCAM", "MIF",
                    "CD45", "SELPLG", "MK")

dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# Open a PNG or TIFF device depending on the file extension
open_device <- function(filename, width, height) {
  path <- file.path(fig_dir, filename)
  if (grepl("\\.tiff$", filename)) {
    tiff(path, width = width, height = height, units = "in", res = 300,
         compression = "lzw")
  } else {
    png(path, width = width, height = height, units = "in", res = 300)
  }
}

# Draw a plot into a file; visible return values (ggplot, heatmap objects) are
# printed, as at the console
save_plot <- function(filename, width, height, plot_fun) {
  open_device(filename, width, height)
  on.exit(dev.off())
  res <- withVisible(plot_fun())
  if (res$visible) print(res$value)
  invisible(NULL)
}

# ---- CellChat object -------------------------------------------------------
seurat_obj <- readRDS(file.path(obj_dir, "annotated_seurat.rds"))

meta <- data.frame(labels = Idents(seurat_obj))
cellchat <- createCellChat(object = seurat_obj, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels")

cell_types <- levels(cellchat@idents)
group_size <- as.numeric(table(cellchat@idents))

# Ligand-receptor database (human)
cellchat@DB <- CellChatDB.human

# ---- Communication inference ----------------------------------------------
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = min_cells)

# Ligand-receptor pairs between plasma cells and all other cell types
plasma_idx <- which(cell_types == "Plasma")
other_idx  <- setdiff(seq_along(cell_types), plasma_idx)

df_net_plasma     <- subsetCommunication(cellchat, sources.use = plasma_idx,
                                         targets.use = other_idx)
df_net_plasma_rev <- subsetCommunication(cellchat, sources.use = other_idx,
                                         targets.use = plasma_idx)
write.xlsx(list(Net_Plasma = df_net_plasma, Net_Plasma_Rev = df_net_plasma_rev),
           file.path(tab_dir, "cellchat_net_plasma.xlsx"))

# Pathway-level communication and aggregated networks
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file.path(obj_dir, "cellchat.rds"))

net_count  <- cellchat@net$count
net_weight <- cellchat@net$weight

pathways_all   <- unique(cellchat@netP$pathways)
pathways_focus <- intersect(pathways_focus, pathways_all)

# ---- Aggregated networks ---------------------------------------------------
save_plot("netVisual_circle_all_count.png", 6, 6, function() {
  netVisual_circle(net_count, vertex.weight = group_size, weight.scale = TRUE,
                   label.edge = FALSE, title.name = "Number of interactions")
})
save_plot("netVisual_circle_all.tiff", 8, 8, function() {
  netVisual_circle(net_weight, vertex.weight = group_size, weight.scale = TRUE,
                   label.edge = FALSE, title.name = "Interaction weights/strength")
})
save_plot("netVisual_circle_all_count.tiff", 8, 8, function() {
  netVisual_circle(net_count, vertex.weight = group_size, weight.scale = TRUE,
                   label.edge = FALSE, title.name = "Number of interactions")
})

# Outgoing signalling of each cell type (other groups masked)
plot_sender_network <- function(net_matrix, cell_type, title) {
  sender <- matrix(0, nrow = nrow(net_matrix), ncol = ncol(net_matrix),
                   dimnames = dimnames(net_matrix))
  sender[cell_type, ] <- net_matrix[cell_type, ]
  netVisual_circle(sender, vertex.weight = group_size, weight.scale = TRUE,
                   edge.weight.max = max(net_matrix), title.name = title,
                   vertex.label.cex = 2.5)
}

for (ct in rownames(net_weight)) {
  save_plot(paste0("net_circle_", ct, ".png"), 10, 10,
            function() plot_sender_network(net_weight, ct, ct))
}
for (ct in rownames(net_count)) {
  save_plot(paste0("net_circle_count_", ct, ".tiff"), 10, 10,
            function() plot_sender_network(net_count, ct,
                                           paste(ct, "(# interactions)")))
}

# ---- Signalling pathways ---------------------------------------------------
for (pw in pathways_focus) {
  save_plot(paste0(pw, "_netVisual_aggregate.tiff"), 8, 8, function() {
    netVisual_aggregate(cellchat, signaling = pw, layout = "chord")
  })
  save_plot(paste0(pw, "_netVisual_heatmap.tiff"), 10, 6, function() {
    netVisual_heatmap(cellchat, signaling = pw, color.heatmap = "Reds")
  })
}

# ---- Ligand-receptor pairs between plasma cells and other cell types -------
for (ext in c("tiff", "png")) {
  save_plot(paste0("netVisual_bubble_pathway_source_plasma.", ext), 10, 10,
            function() netVisual_bubble(cellchat, sources.use = plasma_idx,
                                        targets.use = other_idx,
                                        remove.isolate = FALSE, font.size = 14))
  save_plot(paste0("netVisual_bubble_pathway_source_other_cells_to_plasma.", ext),
            10, 10,
            function() netVisual_bubble(cellchat, sources.use = other_idx,
                                        targets.use = plasma_idx,
                                        remove.isolate = FALSE, font.size = 14))
}

save_plot("netVisual_chord_gene_source_plasma.tiff", 12, 12, function() {
  netVisual_chord_gene(cellchat, sources.use = plasma_idx,
                       targets.use = other_idx, lab.cex = 0.7, legend.pos.y = 30)
})
save_plot("netVisual_chord_gene_source_other_cells_to_plasma.tiff", 12, 12,
          function() {
  netVisual_chord_gene(cellchat, sources.use = other_idx,
                       targets.use = plasma_idx, legend.pos.x = 15,
                       lab.cex = 0.8, small.gap = 5, big.gap = 25,
                       annotationTrackHeight = 0.08)
})

# ---- Signalling roles ------------------------------------------------------
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

for (pw in pathways_focus) {
  save_plot(paste0("netAnalysis_signalingRole_network_", pw, ".tiff"), 8, 8,
            function() netAnalysis_signalingRole_network(
              cellchat, signaling = pw, width = 8, height = 2.5, font.size = 10))
}

save_plot("netAnalysis_signalingRole_heatmap.tiff", 12, 12, function() {
  ht_out <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",
                                              width = 12, height = 12,
                                              color.heatmap = "Purples")
  ht_in  <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",
                                              width = 12, height = 12,
                                              color.heatmap = "Oranges")
  ht_out + ht_in
})

# Hierarchy plot and ligand-receptor contribution for every pathway
for (pw in pathways_all) {
  save_plot(paste0(pw, "_netVisual.png"), 8, 8, function() {
    netVisual(cellchat, signaling = pw, vertex.receiver = vertex_receiver,
              layout = "hierarchy")
  })
  gg <- netAnalysis_contribution(cellchat, signaling = pw)
  ggsave(file.path(fig_dir, paste0(pw, "_L-R_contribution.pdf")),
         plot = gg, width = 5, height = 5, units = "in")
}

# ---- Communication patterns ------------------------------------------------
# Inspect selectK output to choose the number of patterns
selectK(cellchat, pattern = "outgoing")

open_device("identifyCommunicationPatterns_heatmap_outgoing_global.tiff", 12, 12)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing",
                                          k = n_patterns_outgoing,
                                          width = 5, height = 10, font.size = 8)
dev.off()

save_plot("identifyCommunicationPatterns_river_outgoing_global.tiff", 10, 8,
          function() netAnalysis_river(cellchat, pattern = "outgoing"))
save_plot("identifyCommunicationPatterns_dotplot_outgoing_global.tiff", 10, 6,
          function() netAnalysis_dot(cellchat, pattern = "outgoing"))

selectK(cellchat, pattern = "incoming")

open_device("identifyCommunicationPatterns_heatmap_incoming_global.tiff", 12, 12)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming",
                                          k = n_patterns_incoming,
                                          width = 5, height = 10, font.size = 8)
dev.off()

save_plot("identifyCommunicationPatterns_river_incoming_global.tiff", 10, 8,
          function() netAnalysis_river(cellchat, pattern = "incoming"))
save_plot("identifyCommunicationPatterns_dotplot_incoming_global.tiff", 10, 6,
          function() netAnalysis_dot(cellchat, pattern = "incoming"))

# ---- Functional similarity of signalling pathways --------------------------
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

save_plot("netVisual_embedding_functional.tiff", 8, 6, function() {
  netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
})

# ---- Summary tables --------------------------------------------------------
total_incoming <- colSums(net_count)
total_outgoing <- rowSums(net_count)

header_style <- createStyle(textDecoration = "bold")

# Interaction counts and weights between cell types
cell_type_legend <- data.frame(
  Cell_Type = c("Plasma", "CD8_T_cell", "NK", "CD4_T_cell", "CD14_Mono",
                "CD16_Mono", "Late_Eryth", "B_cell", "GMP", "HSC",
                "Prolif_cells", "pDC", "cDC2"),
  Description = c("Plasma cells", "CD8+ cytotoxic T cells",
                  "Natural killer cells", "CD4+ helper T cells",
                  "Classical monocytes", "Non-classical monocytes",
                  "Late erythroblasts", "B cells",
                  "Granulocyte-monocyte progenitors", "Hematopoietic stem cells",
                  "Proliferating cells", "Plasmacytoid dendritic cells",
                  "Conventional dendritic cells type 2")
)

wb <- createWorkbook()

addWorksheet(wb, "CellType_Interaction_Counts_AL")
writeData(wb, "CellType_Interaction_Counts_AL", net_count,
          rowNames = TRUE, colNames = TRUE, headerStyle = header_style)
writeData(wb, "CellType_Interaction_Counts_AL",
          "Number of significant ligand-receptor pairs between cell types.",
          startRow = nrow(net_count) + 3, startCol = 1)

addWorksheet(wb, "CellType_CommProb_AL")
writeData(wb, "CellType_CommProb_AL", net_weight,
          rowNames = TRUE, colNames = TRUE, headerStyle = header_style)
writeData(wb, "CellType_CommProb_AL",
          paste("Interaction weight between cell types (sum of communication",
                "probabilities of significant ligand-receptor pairs). Higher",
                "values indicate stronger signalling."),
          startRow = nrow(net_weight) + 3, startCol = 1)

addWorksheet(wb, "Cell_Type_Legend_AL")
writeData(wb, "Cell_Type_Legend_AL", cell_type_legend,
          headerStyle = createStyle(textDecoration = "bold", fgFill = "#D3D3D3"))

saveWorkbook(wb, file.path(tab_dir, "CellChat_AL_Analysis.xlsx"), overwrite = TRUE)

# Total incoming and outgoing interactions per cell type
comm_table <- rbind(
  data.frame(Metric = "Total_Incoming", t(total_incoming)),
  data.frame(Metric = "Total_Outgoing", t(total_outgoing))
)

wb <- createWorkbook()
addWorksheet(wb, "Comm_Summary_AL")
writeData(wb, "Comm_Summary_AL",
          "Total cell-cell communication (AL amyloidosis)",
          startRow = 1, startCol = 1)
writeData(wb, "Comm_Summary_AL",
          "Total number of significant ligand-receptor pairs per cell type (CellChat count matrix).",
          startRow = 3, startCol = 1)
writeData(wb, "Comm_Summary_AL", comm_table, startRow = 6,
          headerStyle = createStyle(textDecoration = "bold", fgFill = "#D9D9D9"))

saveWorkbook(wb, file.path(tab_dir, "Supplementary_Comm_Summary_AL.xlsx"),
             overwrite = TRUE)
