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


df <- fread("<PROJECT_ROOT>/00.preprocess-within-family-BMI/gwas_catalog/BMI_within_family.txt.gz")

# > unique(df$CHROM)
# [1]  1 10 11 12 13 14 15 16 17 18 19  2 20 21 22  3  4  5  6  7  8  9 23

autosomes <- seq(1,22)

# > autosomes
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22

autosomes_df <- df[df$chromosome %in% autosomes, ]
unique(autosomes_df$chromosome)

# 3. Filter MAF > 0.01  ----
autosomes_df <- autosomes_df[autosomes_df$effect_allele_frequency > 0.01 & autosomes_df$effect_allele_frequency < 0.99, ]


  
  
# LD clumping ----


# Only clump variants with p < 5e-8
autosomes_df <- autosomes_df[autosomes_df$p_value < 5e-8, ]

# add new pval column (needs to be called pval for )
autosomes_df$pval <- as.numeric(autosomes_df$p_value)
typeof(autosomes_df$pval)

# max(autosomes_df$minus_log10p)  # here you see minus log 10 p actually shows min p value is not 0
min(autosomes_df$pval)  # here you see min pval is 0

# Rename for ld_clump_local function
autosomes_df %<>% dplyr::rename("rsid" = "rs_id")


# Initialize empty dataframe to store clumping results for all chromosomes
header <- colnames(autosomes_df)
new_df <- data.frame(matrix(ncol = length(header), nrow = 0))
colnames(new_df) <- header
head(new_df)

for (i in unique(autosomes_df$chromosome)) {
  tryCatch({
    df_clumped <- ld_clump_local(
      dat = autosomes_df,
      bfile = paste('<DATA_ROOT>/data/02.ukb50k_maf0.01/', i, sep = ""),  # UKB 50k with MAF > 0.01
      plink_bin = '<SOFTWARE_ROOT>/plink',
      clump_kb = 1000,
      clump_p = 5e-8,
      clump_r2 = 0.001
    )
    new_df <- rbind(new_df, df_clumped)
    print(new_df)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  
}

# drop pval column since don't need it anymore
new_df$pval <- NULL

# rename back to original column names
new_df <- new_df %>% dplyr::rename("rs_id" = "rsid")






#################
# prepare an output file name
#################

# output_path
directory_path <- "<PROJECT_ROOT>/ld_clumped/"
system(paste0("mkdir -p ", directory_path))
output_basename <- "BMI_within_family.tsv"  # input is a VCF; write a TSV with a sensible name
output_path <- paste0(directory_path, output_basename)

# Write the merged file to a compressed TSV file
write_tsv(new_df, output_path)
    

end_time <- Sys.time()


print(paste0("Time taken: ", end_time - start_time))


