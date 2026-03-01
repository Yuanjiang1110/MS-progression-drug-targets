#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

# =========================
# Inputs (edit here)
# =========================
seurat_rds <- "data/scRNA/ms_single_cell.rds"
outdir     <- "results/scRNA_DE"

# cell types to keep (edit to your study)
cell_types_to_keep <- c(
  "microglial cell",
  "oligodendrocyte",
  "oligodendrocyte precursor cell",
  "astrocyte",
  "endothelial cell",
  "neuron",
  "pericyte",
  "progenitor cell"
)

# target genes (Ensembl IDs in your example)
target_genes <- c(
  "ENSG00000197959", # DNM3 (example)
  "ENSG00000075292", # ZNF638
  "ENSG00000174227", # PIGG
  "ENSG00000135636"  # DYSF
)

# disease labels in metadata
group_var <- "disease"
case_level <- "multiple sclerosis"
ctrl_level <- "normal"

# FindMarkers settings (minimal & robust)
min_pct <- 0.1
test_use <- "wilcox"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Load Seurat object
# =========================
obj <- readRDS(seurat_rds)

# Required metadata checks
meta_cols <- colnames(obj@meta.data)
stopifnot("cell_type" %in% meta_cols)
stopifnot(group_var %in% meta_cols)

# Subset cell types of interest
obj <- subset(obj, subset = cell_type %in% cell_types_to_keep)

# Keep only target genes present in object
target_present <- intersect(target_genes, rownames(obj))
if (length(target_present) == 0) {
  stop("None of target_genes found in Seurat object rownames().")
}

# =========================
# Differential expression per cell type
# =========================
res_list <- list()

for (ct in cell_types_to_keep) {
  
  sub <- subset(obj, subset = cell_type == ct)
  
  # check both groups exist
  groups <- unique(as.character(sub@meta.data[[group_var]]))
  if (!(case_level %in% groups && ctrl_level %in% groups)) {
    message("[SKIP] ", ct, ": missing one of groups (", ctrl_level, ", ", case_level, ")")
    next
  }
  
  # run DE only on target genes (fast + focused)
  mk <- FindMarkers(
    object = sub,
    group.by = group_var,
    ident.1  = case_level,
    ident.2  = ctrl_level,
    features = target_present,
    logfc.threshold = 0,
    min.pct = min_pct,
    test.use = test_use
  )
  
  if (nrow(mk) == 0) next
  
  mk$gene <- rownames(mk)
  mk$cell_type <- ct
  
  # FDR
  mk$p_val_adj_fdr <- p.adjust(mk$p_val, method = "fdr")
  
  res_list[[ct]] <- mk
}

res <- data.table::rbindlist(res_list, fill = TRUE)

# Save
data.table::fwrite(res, file.path(outdir, "MS_vs_normal_DE_target_genes.tsv"), sep = "\t")

# Also save significant results only
res_sig <- res[p_val_adj_fdr < 0.05]
data.table::fwrite(res_sig, file.path(outdir, "MS_vs_normal_DE_target_genes_sig.tsv"), sep = "\t")

message("[DONE] Wrote results to: ", outdir)