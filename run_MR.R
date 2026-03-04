#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

# =========================
# Inputs (edit here)
# =========================
rdata_in <- "data/derived/pQTL_MS_data.Rdata"   # contains exposures + outcomes objects
outdir   <- "results/mr_basic"

# sample sizes (edit to your real Ns)
N_gwas <- 12584
N_exp  <- list(
  plasma = 20013,
  brain  = 376
)

# Optional: clumping parameters
DO_CLUMP  <- TRUE
CLUMP_R2  <- 0.001
CLUMP_KB  <- 10000

# Filter outcome SNPs 
KEEP_OUTCOME_P_GT <- 5e-8

# MR methods (keep minimal)
MR_METHODS <- c("mr_ivw_mre", "mr_wald_ratio")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Load data
# =========================
load(rdata_in)

# Expect the Rdata to contain (names can be adjusted here):
# exposures (lists): plasma_pQTL,  brain_pQTL_rosmap
# outcomes (data.frames): out_ms_plasma,  out_ms_brain

# If your objects are named differently, map them here:
exposure_map <- list(
  plasma = plasma_pQTL,
  brain  = brain_pQTL_rosmap
)

outcome_map <- list(
  plasma = out_ms_plasma,
  brain  = out_ms_brain
)

# =========================
# Helper: run MR for one tissue
# =========================
run_mr_tissue <- function(exposure_list, outcome_df, n_exp, n_gwas,
                          tissue_name,
                          do_clump=TRUE, clump_r2=0.001, clump_kb=10000,
                          keep_outcome_p_gt=5e-8,
                          methods=c("mr_ivw_mre","mr_wald_ratio")) {
  
  # outcome filter
  outcome_df <- outcome_df[outcome_df$pval.outcome > keep_outcome_p_gt, ]
  
  res_list <- vector("list", length(exposure_list))
  harm_list <- vector("list", length(exposure_list))
  steiger_list <- vector("list", length(exposure_list))
  
  # iterate exposures (list names are protein ids; if unnamed, will use index)
  for (i in seq_along(exposure_list)) {
    exp_dat <- exposure_list[[i]]
    
    # optional clump (applied consistently)
    if (do_clump) {
      exp_dat <- tryCatch(
        clump_data(exp_dat, clump_r2 = clump_r2, clump_kb = clump_kb),
        error = function(e) exp_dat  # if clumping fails, continue with unclumped
      )
    }
    
    dat <- harmonise_data(
      exposure_dat = exp_dat,
      outcome_dat  = outcome_df,
      action = 1
    )
    
    if (nrow(dat) == 0) next
    
    dat$samplesize.exposure <- n_exp
    dat$samplesize.outcome  <- n_gwas
    
    mr_res <- mr(dat, method_list = methods)
    mr_or  <- generate_odds_ratios(mr_res)
    
    st <- tryCatch(directionality_test(dat), error = function(e) NULL)
    
    res_list[[i]]    <- mr_or
    harm_list[[i]]   <- dat
    steiger_list[[i]]<- st
  }
  
  res <- rbindlist(res_list, fill = TRUE)
  harm <- rbindlist(harm_list, fill = TRUE)
  stg <- rbindlist(steiger_list, fill = TRUE)
  
  # minimal column selection (similar to your original)
  if (nrow(res) > 0) {
    keep_cols <- intersect(names(res), c(
      "id.exposure","id.outcome","exposure","outcome",
      "method","nsnp","b","se","pval",
      "or","or_lci95","or_uci95"
    ))
    res <- res[, ..keep_cols]
    res[, FDR := p.adjust(pval, method = "BH")]
  }
  
  list(result = res, harmonised = harm, steiger = stg)
}

# =========================
# Run all tissues
# =========================
for (tissue in names(exposure_map)) {
  
  out <- run_mr_tissue(
    exposure_list = exposure_map[[tissue]],
    outcome_df    = outcome_map[[tissue]],
    n_exp         = N_exp[[tissue]],
    n_gwas        = N_gwas,
    tissue_name   = tissue,
    do_clump      = DO_CLUMP,
    clump_r2      = CLUMP_R2,
    clump_kb      = CLUMP_KB,
    keep_outcome_p_gt = KEEP_OUTCOME_P_GT,
    methods       = MR_METHODS
  )
  
  fwrite(out$result,     file.path(outdir, paste0("mr_", tissue, "_result.tsv")), sep = "\t")
  fwrite(out$harmonised, file.path(outdir, paste0("mr_", tissue, "_harmonised.tsv")), sep = "\t")
  fwrite(out$steiger,    file.path(outdir, paste0("mr_", tissue, "_steiger.tsv")), sep = "\t")
  
  message("[DONE] ", tissue, " -> ", outdir)

}
