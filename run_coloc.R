#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

# =========================
# User inputs (edit here)
# =========================
GWAS_N <- 12584
QTL_N  <- 376          # brain pQTL example; change to 7213 for plasma
TYPE1  <- "quant"      # your GWAS is treated as quantitative in your original script
TYPE2  <- "quant"

# Input file must contain BOTH GWAS + QTL columns (see required columns below)
infile  <- "data/coloc_inputs/brain_msseverity_coloc_input.tsv"
outfile <- "results/coloc/brain_msseverity_coloc_summary.tsv"

# Optional: run only selected proteins (set NULL to run all proteins found)
proteins_keep <- NULL
# proteins_keep <- c("DNM3","CAB39L")

# =========================
# Helpers
# =========================
check_cols <- function(dt, required) {
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
}

safe_p <- function(p) {
  p <- as.numeric(p)
  p[p == 0] <- 1e-300
  p
}

run_one_protein <- function(d, gwas_n, qtl_n, type1="quant", type2="quant") {
  # coloc wants one row per SNP
  d <- d[!duplicated(SNP)]
  d <- d[complete.cases(d[, .(SNP, pos, beta_gwas, se_gwas, p_gwas, maf, beta_qtl, se_qtl, p_qtl)])]
  
  if (nrow(d) < 50) return(NULL)  # too few SNPs to be meaningful; adjust threshold if needed
  
  d[, `:=`(
    varbeta_gwas = se_gwas^2,
    varbeta_qtl  = se_qtl^2,
    p_gwas = safe_p(p_gwas),
    p_qtl  = safe_p(p_qtl)
  )]
  
  fit <- coloc.abf(
    dataset1 = list(
      beta     = d$beta_gwas,
      varbeta  = d$varbeta_gwas,
      snp      = d$SNP,
      position = d$pos,
      pvalues  = d$p_gwas,
      type     = type1,
      N        = gwas_n,
      MAF      = d$maf
    ),
    dataset2 = list(
      beta     = d$beta_qtl,
      varbeta  = d$varbeta_qtl,
      snp      = d$SNP,
      position = d$pos,
      pvalues  = d$p_qtl,
      type     = type2,
      N        = qtl_n,
      MAF      = d$maf
    )
  )
  
  s <- as.data.table(as.list(fit$summary))
  s
}

# =========================
# Main
# =========================
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

dt <- fread(infile)

# Required columns in infile:
# SNP, protein, pos, maf,
# beta_gwas, se_gwas, p_gwas,
# beta_qtl,  se_qtl,  p_qtl
required_cols <- c(
  "SNP","protein","pos","maf",
  "beta_gwas","se_gwas","p_gwas",
  "beta_qtl","se_qtl","p_qtl"
)
check_cols(dt, required_cols)

# optional subset
if (!is.null(proteins_keep)) {
  dt <- dt[protein %in% proteins_keep]
}

# split by protein and run coloc
res_list <- lapply(split(dt, dt$protein), function(d) {
  out <- run_one_protein(d, gwas_n = GWAS_N, qtl_n = QTL_N, type1 = TYPE1, type2 = TYPE2)
  if (is.null(out)) return(NULL)
  out[, protein := unique(d$protein)]
  out
})

res <- rbindlist(res_list, fill = TRUE)
setcolorder(res, c("protein", setdiff(names(res), "protein")))

fwrite(res, outfile, sep = "\t")

message("[DONE] Wrote: ", outfile)