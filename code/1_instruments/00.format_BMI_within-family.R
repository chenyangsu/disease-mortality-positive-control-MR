# Format to GWAS catalog columns, Label cis/trans, LD clumping

start_time <- Sys.time()

# .libPaths(c("<R_LIBS_USER>",
#             "<R_LIBS_SITE>",
#             "<R_HOME>/library"))

library(magrittr)
library(dplyr)
library(tidyverse)
library(readxl)
library(data.table)
library(tidyr)
library(stringr)
library(ieugwasr)
library(gwasvcf)             # vcf_to_granges()
library(VariantAnnotation)   # readVcf()

# args = commandArgs(trailingOnly=TRUE)
# 
# # test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }
# 
# path <- args[1]  


path <- "<PROJECT_ROOT>/claude_output/ieu-b-4815.vcf.gz"


# # 2. Retain variants on autosomes only ----


# GWAS-VCF (ieu-b-4815, Howe et al. 2022 within-sibship BMI) → GWAS-Catalog schema
# FORMAT = ES:SE:LP:SS:ID  (effect size, SE, -log10 P, sample size, study SNP ID)
# INFO   = AF              (effect-allele frequency, where present)
vcf <- VariantAnnotation::readVcf(path)
df  <- gwasvcf::vcf_to_granges(vcf) %>% as_tibble()
head(df)

std_df <- data.frame(

  # 0  — coerce factor/Rle to integer so downstream `chromosome %in% 1:22` and
  #      `chromosome == 6` work cleanly (NA for X/Y/MT, dropped by autosome filter)
  chromosome = suppressWarnings(as.integer(as.character(df$seqnames))),

  # 1
  base_pair_location = df$start,

  # 2
  effect_allele = df$ALT,

  # 3
  other_allele = df$REF,

  # 4
  beta = df$ES,

  # 5
  standard_error = df$SE,

  # 6
  effect_allele_frequency = df$AF,

  # 7
  p_value = 10^(-df$LP),

  # Encouraged
  variant_id = NA,
  rs_id = df$ID,
  info = NA,
  ci_lower = NA,
  ci_upper = NA,
  ref_allele = NA,
  n = df$SS,

  stringsAsFactors = FALSE
)
df <- as.data.table(std_df)
rm(vcf, std_df)
gc()
head(df)


##################### speed up everything ######################
str(df)
df <- df[df$p_value < 5e-8, ]
################################################################

####################

# Remove HLA region
# chr6:28,510,120-33,480,577

df <- df[!(df$chromosome == 6 & df$base_pair_location < 33480577 & df$base_pair_location > 28510120), ]

#####################

unique(df$chromosome)

write_tsv(df, "<PROJECT_ROOT>/00.preprocess-within-family-BMI/within-family-BMI-no-eaf.tsv")
