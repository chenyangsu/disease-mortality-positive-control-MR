# Format to GWAS catalog columns, Label cis/trans, LD clumping

start_time <- Sys.time()

.libPaths(c("<R_LIBS_USER>",
            "<R_LIBS_SITE>",
            "<R_HOME>/library"))

library(magrittr)
library(dplyr)
library(tidyverse)
library(readxl)
library(data.table)
library(tidyr)
library(stringr)
library(ieugwasr)

# args = commandArgs(trailingOnly=TRUE)
# 
# # test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }
# 
# path <- args[1]  


path <- "<DATA_ROOT>/data/04.prep_MR_outcome/socioeconomic/gwas_catalog/Household.catalog.tsv.gz"


# # 2. Retain variants on autosomes only ----


df <- fread(path)
head(df)

# > head(df)
# chromosome base_pair_location effect_allele other_allele     beta standard_error effect_allele_frequency p_value variant_id       rs_id   info ci_lower ci_upper ref_allele       n
# <int>              <int>        <char>       <char>    <num>          <num>                   <num>   <num>     <lgcl>      <char> <lgcl>   <lgcl>   <lgcl>     <char>   <int>
#   1:          5           29439168             T            C -0.00032        0.00179               0.3765480 0.86000         NA    rs667647     NA       NA       NA        #NA 3037499
# 2:          5           86633075             T            C -0.00105        0.00359               0.0746997 0.76960         NA rs113534962     NA       NA       NA        #NA 3037499
# 3:          2          170110443             T            C  0.00724        0.00939               0.9893748 0.44090         NA rs559397866     NA       NA       NA        #NA 3037499
# 4:         10          126652663             T            C  0.00336        0.00174               0.4646600 0.05336         NA   rs2366866     NA       NA       NA        #NA 3037499
# 5:          5           46390943             A            C  0.01236        0.02326               0.9988451 0.59520         NA rs540077909     NA       NA       NA        #NA 3037499
# 6:          3           62721844             T            C -0.00748        0.00372               0.0610560 0.04447         NA    rs472303     NA       NA       NA        #NA 3037499


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

for (i in 1:22) {
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
output_path <- paste0(directory_path, basename(path))

# Write the merged file to a compressed TSV file
write_tsv(new_df, output_path)
    

end_time <- Sys.time()

print(paste0("Time taken: ", end_time - start_time))

