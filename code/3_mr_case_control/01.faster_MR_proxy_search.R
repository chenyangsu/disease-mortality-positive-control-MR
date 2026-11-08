#######################################################################
# 01.faster_MR_proxy_search.R — BMI -> 6 disease outcomes
#
# Repurposed from `1-myanalysis/01.faster_MR_proxy_search.R`.
# Outer loop is the disease outcome (selected by the array index `i`),
# inner loop is over exposure files in `path_to_exposures`. Here the
# exposure dir contains only `BMI.30124842.catalog.tsv.gz` (the LD-
# clumped lead SNPs), so a single exposure is iterated per outcome.
#
# Args:
#   1. i (integer, 1..6) — row in outcome.all.tsv to run
#   2. path_to_exposures (dir) — dir of LD-clumped exposure tables
#   3. outdir (dir)            — top-level results dir
#######################################################################

start_time <- Sys.time()

suppressPackageStartupMessages({
  library(tidyverse)
  library(TwoSampleMR)
  library(stringr)
  library(gwasvcf)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 01.faster_MR_proxy_search.R <i> <path_to_exposures> <outdir>")
}
i                 <- as.integer(args[1])
path_to_exposures <- args[2]
outdir            <- args[3]

OUTCOME_TSV  <- "<PROJECT_ROOT>/02.BMI_MR_6diseases/outcome.all.tsv"
PLINK_PATH   <- "<SOFTWARE_ROOT>/plink"
# UKB 50k MAF>0.01 per-chromosome bfile prefix used by the original
# BMI LD-clumping script. get_ld_proxies() expects `<prefix><chr>`.
BFILE_PREFIX <- "<DATA_ROOT>/data/02.ukb50k_maf0.01/"

cat("Outcome table   :", OUTCOME_TSV, "\n")
cat("Exposure dir    :", path_to_exposures, "\n")
cat("Output dir      :", outdir, "\n")
cat("PLINK           :", PLINK_PATH, "\n")
cat("Bfile prefix    :", BFILE_PREFIX, "\n\n")

# -- 1. Pick the i-th outcome and load it ----------------------------------
outcome_path     <- fread(OUTCOME_TSV, fill = TRUE, sep = "\t")
ith_outcome_path <- outcome_path[i, ]
ith_outcome_name <- ith_outcome_path$outcome_name
ith_outcome_df   <- as.data.frame(fread(ith_outcome_path$gwas_path, fill = TRUE))

# Defensive cleaning: strip stray CR on rs_id (CRC file), drop NA/blank rsids
# (HF has ~189k NAs that otherwise leak NA rows into base-R logical subsetting),
# and dedup on rs_id+alleles so downstream `&&` checks always see length-1 rows.
n_raw <- nrow(ith_outcome_df)
ith_outcome_df$rs_id <- sub("\r$", "", ith_outcome_df$rs_id)
ith_outcome_df <- ith_outcome_df[!is.na(ith_outcome_df$rs_id) &
                                 ith_outcome_df$rs_id != "", ]
ith_outcome_df <- ith_outcome_df[
  !duplicated(ith_outcome_df[, c("rs_id", "effect_allele", "other_allele")]), ]
cat("Outcome cleaning: ", n_raw, " -> ", nrow(ith_outcome_df),
    " rows (dropped NA/blank rsids and rsid+allele duplicates)\n", sep = "")

cat("Job Array Number:", i, "- Outcome:", ith_outcome_name, "\n")
cat("Outcome rows    :", nrow(ith_outcome_df), "\n")
if ("p_value" %in% names(ith_outcome_df)) {
  cat("Min outcome p   :", suppressWarnings(min(ith_outcome_df$p_value, na.rm = TRUE)), "\n")
}

# -- 2. Discover exposure files --------------------------------------------
exposure_paths <- list.files(path_to_exposures, full.names = TRUE,
                             pattern = "\\.tsv(\\.gz)?$|\\.txt(\\.gz)?$")
cat("---- Running loop for :", length(exposure_paths), "exposures ----\n\n")

# -- 3. PLINK once for the whole run ---------------------------------------
set_plink(path = PLINK_PATH)

for (j in seq_along(exposure_paths)) {
  tryCatch({

    exp_path <- exposure_paths[j]
    # Use the first two dot-separated chunks: "BMI.30124842.catalog.tsv.gz" -> "BMI"
    # (matches the convention in the original 1-myanalysis script.)
    protname <- strsplit(basename(exp_path), "\\.")[[1]][1]
    cat("Analyzing exposure", j, ":", protname, "----\n")

    Ins_df <- as.data.frame(fread(exp_path, fill = TRUE))

    Ins <- format_data(
      Ins_df, type = "exposure", header = TRUE,
      snp_col = "rs_id", beta_col = "beta",
      se_col = "standard_error", eaf_col = "effect_allele_frequency",
      effect_allele_col = "effect_allele", other_allele_col = "other_allele",
      pval_col = "p_value", min_pval = 1e-300, samplesize_col = "n"
    )

    ith_outcome <- ith_outcome_df %>% filter(rs_id %in% Ins$SNP)

    # ---- Proxy search ------------------------------------------------------
    proxy_list <- Ins$SNP[!Ins$SNP %in% ith_outcome$rs_id]
    cat("Need proxies for ", length(proxy_list), " of ", nrow(Ins),
        " exposure SNPs missing from outcome\n", sep = "")

    for (rsid in proxy_list) {

      chr   <- Ins_df[Ins_df$rs_id == rsid, ]$chromosome
      bfile <- paste0(BFILE_PREFIX, chr)

      raw_proxies <- tryCatch(
        get_ld_proxies(
          rsid, bfile,
          tag_kb = 5000, tag_nsnp = 5000, tag_r2 = 0.8,
          threads = 4, out = tempfile()
        ),
        error = function(e) {
          cat("  proxy search error for ", rsid, ": ", conditionMessage(e), "\n", sep = "")
          NULL
        }
      )
      if (is.null(raw_proxies) || nrow(raw_proxies) == 0) next

      raw_proxies$R_sq  <- raw_proxies$R^2
      raw_proxies       <- raw_proxies[raw_proxies$MAF_A <= 0.42 &
                                       raw_proxies$MAF_B <= 0.42, ]
      if (nrow(raw_proxies) == 0) next

      # Walk proxies in descending R^2 (which is the get_ld_proxies output order).
      for (idx in seq_len(nrow(raw_proxies))) {
        proxy_row <- raw_proxies[idx, ]
        if (!(proxy_row$SNP_B %in% ith_outcome_df$rs_id)) next

        Ins_row <- Ins[Ins$SNP == rsid, ]

        outcome_row_a <- ith_outcome_df[
          ith_outcome_df$rs_id == proxy_row$SNP_B &
          ith_outcome_df$effect_allele == proxy_row$B1 &
          ith_outcome_df$other_allele  == proxy_row$B2, ]
        outcome_row_b <- ith_outcome_df[
          ith_outcome_df$rs_id == proxy_row$SNP_B &
          ith_outcome_df$effect_allele == proxy_row$B2 &
          ith_outcome_df$other_allele  == proxy_row$B1, ]

        new_proxy_row <- NULL
        if (Ins_row$effect_allele.exposure == proxy_row$A1) {
          if (!is_empty(outcome_row_a$effect_allele) &&
              proxy_row$B1 == outcome_row_a$effect_allele) {              # Case 1
            new_proxy_row <- outcome_row_a %>%
              mutate(rs_id = proxy_row$SNP_A,
                     effect_allele = Ins_row$effect_allele.exposure,
                     other_allele  = Ins_row$other_allele.exposure,
                     effect_allele_frequency = Ins_row$eaf.exposure)
          } else if (!is_empty(outcome_row_b$effect_allele) &&
                     proxy_row$B2 == outcome_row_b$effect_allele) {       # Case 2
            new_proxy_row <- outcome_row_b %>%
              mutate(rs_id = proxy_row$SNP_A,
                     effect_allele = Ins_row$other_allele.exposure,
                     other_allele  = Ins_row$effect_allele.exposure,
                     effect_allele_frequency = 1 - Ins_row$eaf.exposure)
          }
        } else if (Ins_row$effect_allele.exposure == proxy_row$A2) {
          if (!is_empty(outcome_row_a$effect_allele) &&
              proxy_row$B1 == outcome_row_a$effect_allele) {              # Case 3
            new_proxy_row <- outcome_row_a %>%
              mutate(rs_id = proxy_row$SNP_A,
                     effect_allele = Ins_row$other_allele.exposure,
                     other_allele  = Ins_row$effect_allele.exposure,
                     effect_allele_frequency = 1 - Ins_row$eaf.exposure)
          } else if (!is_empty(outcome_row_b$effect_allele) &&
                     proxy_row$B2 == outcome_row_b$effect_allele) {       # Case 4
            new_proxy_row <- outcome_row_b %>%
              mutate(rs_id = proxy_row$SNP_A,
                     effect_allele = Ins_row$effect_allele.exposure,
                     other_allele  = Ins_row$other_allele.exposure,
                     effect_allele_frequency = Ins_row$eaf.exposure)
          }
        }
        if (!is.null(new_proxy_row)) {
          ith_outcome <- rbind(ith_outcome, new_proxy_row)
          break
        }
      }
    }

    if (nrow(ith_outcome) != nrow(Ins)) {
      cat("After proxy search: ", nrow(ith_outcome), " / ", nrow(Ins),
          " exposure SNPs resolved\n", sep = "")
    }

    # ---- MR ----------------------------------------------------------------
    if (nrow(ith_outcome) == 0) {
      cat("Instrument not found in outcome, skipping MR...\n\n")
    } else {

      outcome_dat <- format_data(
        dat = ith_outcome, type = "outcome",
        snp_col = "rs_id",
        beta_col = "beta", se_col = "standard_error",
        eaf_col = "effect_allele_frequency",
        effect_allele_col = "effect_allele", other_allele_col = "other_allele",
        pval_col = "p_value", min_pval = 1e-300,
        chr_col = "chromosome", pos_col = "base_pair_location",
        log_pval = FALSE
      )
      outcome_dat$samplesize.outcome <- ith_outcome_path[1, ]$sample_size

      exp_dat_outcome <- harmonise_data(exposure_dat = Ins, outcome_dat = outcome_dat)

      if (nrow(exp_dat_outcome[exp_dat_outcome$mr_keep == "TRUE", ]) == 0) {
        cat("EMPTY:", protname, "__", ith_outcome_name, "\n")
      } else {
        output_dir <- paste0(outdir, "/", protname, "/")
        filename   <- paste0(protname, "--", ith_outcome_name)
        for (sub in c("harmonized", "or", "hetero", "pleio", "steiger")) {
          dir.create(paste0(output_dir, sub), recursive = TRUE, showWarnings = FALSE)
        }
        setwd(output_dir)

        write_tsv(exp_dat_outcome,
                  file = paste0(output_dir, "harmonized/", filename, ".harmonized.txt"))

        mr_results <- mr(exp_dat_outcome)
        OR <- generate_odds_ratios(mr_results)
        write_tsv(OR, file = paste0(output_dir, "or/", filename, ".or.txt"))

        pleio_res <- mr_pleiotropy_test(exp_dat_outcome)
        if (nrow(pleio_res) > 0) {
          write_tsv(pleio_res, file = paste0(output_dir, "pleio/", filename, ".pleio.txt"))
        }

        hetero_res <- mr_heterogeneity(exp_dat_outcome)
        if (nrow(hetero_res) > 0) {
          hetero_res$i_squared <- (hetero_res$Q - hetero_res$Q_df) / hetero_res$Q
          hetero_res$i_squared <- ifelse(hetero_res$i_squared < 0, 0, hetero_res$i_squared)
          write_tsv(hetero_res, file = paste0(output_dir, "hetero/", filename, ".hetero.txt"))
        }

        steiger <- directionality_test(exp_dat_outcome)
        write_tsv(steiger, file = paste0(output_dir, "steiger/", filename, ".steiger.txt"))

        cat(filename, ": done\n\n")
      }
    }
  }, error = function(e) cat("ERROR :", conditionMessage(e), "\n\n"))
}

cat("Time taken:", format(Sys.time() - start_time), "\n")
