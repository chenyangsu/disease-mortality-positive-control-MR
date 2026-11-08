####
######### Fast because outer loop is GWAS disease outcome (big file) and in this script, we read in the exposure file (strict v2g cis)
####

start_time <- Sys.time()

.libPaths(c("<R_LIBS_USER>",
            "<R_LIBS_SITE>",
            "<R_HOME>/library"))

###############
# For compute canada
##############

library(tidyverse)
#library(vroom)
library(TwoSampleMR)
library(stringr)
library(gwasvcf)
library(data.table)


args <- commandArgs(trailingOnly=TRUE) #input UNZIPPED file
i <- as.integer(args[1])  # outcome

outcome_path <- fread('<DATA_ROOT>/data/04.prep_MR_outcome/disease_mortality/outcome.all.tsv', fill = T, sep = '\t')
# for (i in 1:nrow(outcome_path)) {
  ith_outcome_path <- outcome_path[i,] 
  ith_outcome_name <- paste0(str_split(basename(ith_outcome_path$gwas_path),'[.]')[[1]][1]) # e.g., "Age_Abdomen.35418184"
  ith_outcome_df <- fread(ith_outcome_path$gwas_path, fill = T)
  cat("Job Array Number:", i, "- Outcome:", ith_outcome_name, "\n")
  print(min(ith_outcome_df$p_value))
# }



# 2. Read in file with all variants to use as exposure.
path_to_exposures <- "<PROJECT_ROOT>/ld_clumped/"
exposure_paths = list.files(path_to_exposures, full.names = T)

# 3. Give directory to store output of MR results
# outdir <- args[3]
outdir <- '<PROJECT_ROOT>/output/'
cat("Output directory:", outdir, '\n')

# 4. Give directory where the proxies for each exposure's variant is
# proxy_dir <- args[4]
#proxy_dir <- '<DATA_ROOT>/01.pQTL/01.aric/05.MR_strict_v2g/0.find_proxies/found_proxies/'
# cat("Proxy directory:", proxy_dir, '\n\n')


cat("---- Running loop for :", length(exposure_paths), "exposures ---- \n\n")


for(j in 1:length(exposure_paths)){
  tryCatch({ # try catch outcome error
    
    
    exp_path <- exposure_paths[j]
    protname <- strsplit(basename(exp_path), "\\.")[[1]][1]  # get first part of string: "EducationalAttainment.35361970.catalog.tsv.gz"
    # protname <- "EducationalAttainment"
    

    Ins_df <- fread(exp_path, fill  = T)  # vroom may read T as TRUE so use fread
    
    Ins<-format_data(Ins_df, type = "exposure", header = TRUE,
                     #phenotype_col = "", 
                     snp_col = "rs_id", beta_col = "beta",
                     se_col = "standard_error", eaf_col = "effect_allele_frequency", effect_allele_col = "effect_allele",
                     other_allele_col = "other_allele", pval_col = "p_value", min_pval = 1e-300, samplesize_col = 'n')

    # retain variants in exposure only
    ith_outcome <- ith_outcome_df %>% filter(rs_id %in% Ins$SNP)
    
    
    ################### Perform proxy search if input GWAS cis pQTLs are not in outcome GWAS ##################################
    # MR proxy search
    # Step 1. Check if variant is in outcome GWAS. If yes, then good and don't do anything, if not then do 2.
    # Step 2. Perform proxy search with gwasvcf get_ld_proxies() function and filter out rows where the variant OR proxy of variant has MAF > 0.42
    # Step 3. If proxy is in the outcome GWAS, relabel effect allele frequency (EAF) such that the minor allele is the minor allele
    
    
    
    # Step 1. ----
    proxy_list <- Ins$SNP[!Ins$SNP %in% ith_outcome$rs_id]  # list of cis-pQTLs that we need to find proxies for
    cat("List of variants that we need to find proxies for: ", proxy_list, "\n")
    
    # Step 2. ----
    # Get LD proxy
    # https://mrcieu.github.io/gwasvcf/reference/get_ld_proxies.html
    set_plink(path = "<SOFTWARE_ROOT>/plink")  # else get_ld_proxy() may give error
    
    for (rsid in proxy_list) {  # get proxy for each rsid that is present in cis-pQTLs but missing from outcome gwas rsid
      
      ################### We get to skip this rate limiting step ##############################
      chr <- Ins_df[Ins_df$rs_id == rsid, ]$chromosome  # extract chromosome
      bfile <- paste0("<DATA_ROOT>/data/02.ukb50k_maf0.01/", chr)

      raw_proxies <- get_ld_proxies(
        rsid,  # single rsid
        bfile,  # plink file corresponding to chromosome of the variant
        tag_kb = 5000,
        tag_nsnp = 5000,
        tag_r2 = 0.8,
        threads = 4,
        out = tempfile()
      )
      raw_proxies$R_sq <- raw_proxies$R^2  # compute rsq ourselves since function only returns R

      # Example of first 2 rows of output (found 30 proxies) when querying rs564382789 NOTE: *_A are columns of original variant, *_B are proxy variant columns

      # CHR_A     BP_A SNP_A       MAF_A CHR_B     BP_B SNP_B      PHASE MAF_B     R A1    B1    A2    B2     R_sq
      # <int>    <int> <chr>       <dbl> <int>    <int> <chr>      <chr> <dbl> <dbl> <chr> <chr> <chr> <chr> <dbl>
      #     19 39419059 rs564382789 0.401    19 39417652 rs10425217 CTTC  0.390 0.985 C     T     T     C     0.970
      #     19 39419059 rs564382789 0.401    19 39417440 rs7508539  CGTA  0.390 0.985 C     G     T     A     0.970


      # Filter out proxies where MAF > 0.42 (equivalently, retain proxies with MAF <= 0.42)
      raw_proxies <- raw_proxies[raw_proxies$MAF_A <= 0.42 & raw_proxies$MAF_B <= 0.42, ]
      #########################################################################
      
      # # just load in the raw proxies file with the proxies for this rsid that we found already
      # raw_proxies <- fread(paste0(proxy_dir, "/", rsid, "_raw_proxies.txt"))
      
      
      # If no proxies are found, go to next rsid that we need to find a proxy for. This is a GOTCHA!
      # Without this check, if nrow(proxies) == 0, the statement "if(proxy_row$SNP_B %in% ith_outcome_df$rs_id) {"  will throw the following error:
      #### Error in if (proxy_row$SNP_B %in% ith_outcome_df$rs_id) { : 
      #### argument is of length zero
      # This will cause proxy search for the list of rsids in "proxy_list" to terminate early and stop all the code below including MR from running!
      if (nrow(raw_proxies) == 0) {  
        # NOTE (added for public release): `next` was intended here. As run, this `break` ends the
        # proxy search at the first instrument with no usable proxy, so no proxies were substituted
        # for the mortality outcomes. Left unchanged so the code reproduces the published results;
        # see Methods (Instrument selection and harmonisation) and Discussion of the paper.
        break
      }
      
      # Step 3. ----
      for (idx in 1:nrow(raw_proxies)) {  # for each proxy, starting from highest R2
        proxy_row <- raw_proxies[idx, ]
        
        if(proxy_row$SNP_B %in% ith_outcome_df$rs_id) {  # check if the proxy is in the outcome GWAS variants. If not, do nothing and go to next proxy
          
          Ins_row <- Ins[Ins$SNP == rsid, ]  # get row for the queried rsid
          
          
          # Extract row of the proxy from the outcome gwas
          outcome_row_a <- ith_outcome_df[ith_outcome_df$rs_id == proxy_row$SNP_B &  # match on proxy rsid from outcome gwas
                                            ith_outcome_df$effect_allele == proxy_row$B1 &  # need to match on effect/other allele since the proxy rsid may have duplications in outcome gwas
                                            ith_outcome_df$other_allele == proxy_row$B2,]
          
          outcome_row_b <- ith_outcome_df[ith_outcome_df$rs_id == proxy_row$SNP_B &  # match on proxy rsid from outcome gwas
                                            ith_outcome_df$effect_allele == proxy_row$B2 &  # need to match on effect/other allele since the proxy rsid may have duplications in outcome gwas
                                            ith_outcome_df$other_allele == proxy_row$B1,]
          
          
          # Need to check 4 cases
          # Case 1: Input GWAS EA matches Proxy queried A1 AND Proxy B1 matches Outcome EA
          # Case 2: Input GWAS EA matches Proxy queried A1 AND Proxy B2 matches Outcome EA
          # Case 3: Input GWAS EA matches Proxy queried A2 AND Proxy B1 matches Outcome EA
          # Case 4: Input GWAS EA matches Proxy queried A2 AND Proxy B2 matches Outcome EA
          
          
          # Case 1 + 2 ----
          if (Ins_row$effect_allele.exposure == proxy_row$A1) {
            if (!is_empty(outcome_row_a$effect_allele) && proxy_row$B1 == outcome_row_a$effect_allele) {  # Case 1: effect allele is the minor allele ("A1" column is always the minor allele of the queried rsid for get_ld_proxies() function)
              
              # Replace 4 columns: rsid, effect allele, other allele, effect allele frequency
              new_proxy_row <- outcome_row_a %>%
                mutate(
                  rs_id = proxy_row$SNP_A,  # replace outcome proxy rsid with the queried variant's rsid
                  effect_allele = Ins_row$effect_allele.exposure,         # replace outcome proxy EA with the queried variant's EA
                  other_allele = Ins_row$other_allele.exposure,  # replace outcome proxy OA with the queried variant's OA
                  effect_allele_frequency = Ins_row$eaf.exposure
                )
              
            } else if (!is_empty(outcome_row_b$effect_allele) && proxy_row$B2 == outcome_row_b$effect_allele) {  # Case 2
              
              # Replace 4 columns: rsid, effect allele, other allele, effect allele frequency
              new_proxy_row <- outcome_row_b %>%
                mutate(
                  rs_id = proxy_row$SNP_A,  # replace outcome proxy rsid with the queried variant's rsid
                  effect_allele = Ins_row$other_allele.exposure,         # replace outcome proxy EA with the queried variant's OA
                  other_allele = Ins_row$effect_allele.exposure,  # replace outcome proxy OA with the queried variant's EA
                  effect_allele_frequency = 1 - Ins_row$eaf.exposure
                )
            }
            
            # For case 3 + 4 ----
          } else if (Ins_row$effect_allele.exposure == proxy_row$A2) {
            if (!is_empty(outcome_row_a$effect_allele) && proxy_row$B1 == outcome_row_a$effect_allele) {   # Case 3
              
              # Replace 4 columns: rsid, effect allele, other allele, effect allele frequency
              new_proxy_row <- outcome_row_a %>%
                mutate(
                  rs_id = proxy_row$SNP_A,  # replace outcome proxy rsid with the queried variant's rsid
                  effect_allele = Ins_row$other_allele.exposure,         # replace outcome proxy EA with the queried variant's OA
                  other_allele = Ins_row$effect_allele.exposure,  # replace outcome proxy OA with the queried variant's EA
                  effect_allele_frequency = 1 - Ins_row$eaf.exposure
                )
              
            } else if (!is_empty(outcome_row_b$effect_allele) && proxy_row$B2 == outcome_row_b$effect_allele) {  # Case 4
              
              # Replace 4 columns: rsid, effect allele, other allele, effect allele frequency
              new_proxy_row <- outcome_row_b %>%
                mutate(
                  rs_id = proxy_row$SNP_A,  # replace outcome proxy rsid with the queried variant's rsid
                  effect_allele = Ins_row$effect_allele.exposure,         # replace outcome proxy EA with the queried variant's EA
                  other_allele = Ins_row$other_allele.exposure,  # replace outcome proxy OA with the queried variant's OA
                  effect_allele_frequency = Ins_row$eaf.exposure
                )
            }
            
          }
          
          ith_outcome <- rbind(ith_outcome, new_proxy_row)  # add the new found proxy row to ith_outcome
          
          break  # found proxy so break out of inner for loop (proxy loop) and go to next rsid that we need to find proxy for
        }  # close check if proxy is present in outcome GWAS: if(proxy_row$SNP_B %in% ith_outcome_df$rs_id) {
      }  # close loop for proxies: for (idx in 1:nrow(raw_proxies)) {
    }  # close loop for rsd:   for (rsid in proxy_list) {
    
    
    # if none of the found proxies for these queried variants were in the outcome gwas
    if (nrow(ith_outcome) == 0) {
      print(paste0("None of the ", nrow(raw_proxies), " proxies for rsid ", rsid, " were found in the outcome GWAS ..."))
      
    } else if (nrow(ith_outcome) != nrow(Ins)) {  # One of the queried rsid did not get replaced by a proxy since no proxies were present in outcome
      print(paste0("For ", length(proxy_list), " queried rsids, only ", nrow(Ins) - nrow(ith_outcome), " proxies were found in the outcome GWAS ..."))
    }
    #################### Finish proxy search
    
    
    
  
    if (nrow(ith_outcome) == 0) {  # e.g. if exposure has single IV which isn't in outcome GWAS and proxy search of this IV also returns nothing then no instrument available for MR
      cat("Instrument not found in outcome, skipping MR...\n\n")  
    } else {  # Run MR
      # format outcome data
      outcome_dat <- format_data(dat = ith_outcome,
                                 type = 'outcome', #phenotype_col = 'id',
                                 snp_col = 'rs_id',
                                 beta_col = 'beta', se_col = 'standard_error', eaf_col = 'effect_allele_frequency', 
                                 effect_allele_col = 'effect_allele', other_allele_col = 'other_allele', pval_col = 'p_value',
                                 #samplesize_col = 'SS',
                                 #id_col = 'id',
                                 min_pval = 1e-300, 
                                 chr_col = 'chromosome', pos_col = 'base_pair_location', log_pval = FALSE)
      # note exposure's smaple size is present but outcome sample size is not
      outcome_dat$samplesize.outcome <- ith_outcome_path[1,]$sample_size
      
      ############
      ##harmonise the exposure and outcome data
      # note this can be zero for some protein-outcome combination (that's why we need if else to escape skip combination)
      # if exp_dat_outcome is not empty, we'll go forward. Otherwise (else) do nothing. 
      ############
      exp_dat_outcome <- NULL
      exp_dat_outcome <- harmonise_data( 
        exposure_dat = Ins, 
        outcome_dat = outcome_dat)
      if(nrow(exp_dat_outcome[exp_dat_outcome$mr_keep=='TRUE',])==0){ # no variant to be tested in MR
        print(paste0('EMPTY:', protname, '__', ith_outcome_name))
        # do nothing
      } else{ # only if exp_dat_outcome is not emply, proceed to MR
        
        # create an output_dir
        output_dir <- paste0(outdir, '/', protname, '/')
        filename <- paste0(protname, '--', ith_outcome_name) # to make naming files easy in the downstream analyses
        system(paste0('mkdir -p ', output_dir, 'harmonized/'))
        system(paste0('mkdir -p ', output_dir, 'or/'))
        system(paste0('mkdir -p ', output_dir, 'hetero/'))
        system(paste0('mkdir -p ', output_dir, 'pleio/'))
        system(paste0('mkdir -p ', output_dir, 'steiger/'))
        setwd(output_dir) # setwd 
        
        exp_data_outcome_name <- paste0(output_dir,"harmonized/", filename, ".harmonized.txt")

        write_tsv(exp_dat_outcome, file=exp_data_outcome_name)

        # mr result
        mr_results <- mr(exp_dat_outcome)
        #mr_name <- paste0(output_dir, "mr/", protname, ".mr.txt")
        #write.table(mr_results, file=mr_name, sep = '\t', quote = F)

        # odds ratio
        OR <- generate_odds_ratios(mr_results)
        OR_name <- paste0(output_dir,"or/", filename, ".or.txt")
        write_tsv(OR, file=OR_name)

        # # scatter plot
        # pdf_name <- paste0(output_dir, "pdf/", protname, ".scatter.pdf")
        # pdf(pdf_name)
        # mr_scatter_plot(mr_results, exp_dat_outcome)[[1]] + xlim(0,1)
        # dev.off()

        # horizontal pleiotropy
        pleio_res <- mr_pleiotropy_test(exp_dat_outcome)
        if (nrow(pleio_res) > 0) {
          pleio_name <- paste0(output_dir,"pleio/", filename, ".pleio.txt")
          write_tsv(pleio_res, file=pleio_name)
        }
        
        # hetero test
        hetero_res <- mr_heterogeneity(exp_dat_outcome)
        if (nrow(hetero_res) > 0) {
          
          # res_single <- mr_singlesnp(exp_dat_outcome)   # res_single includes IVW and MR-Egger in the last two rows
          # res_single_beta <- res_single$b[1:nrow(exp_dat_outcome)]
          # res_single_se <- res_single$se[1:nrow(exp_dat_outcome)]
          # res_Isq <- Isq(res_single_beta, res_single_se)
          # then add to hetero_res$i_squared
          
          # Only compute I^2 when Q statistic is > Q_df and p value is significant (p<0.05) (indicating heterogeneity)
          # if Q < Q_df then I^2 will be below 0 which does not makes sense. Also Q<1 indicates no heterogeneity so computing I^2 is nonsensical
          hetero_res$i_squared <- (hetero_res$Q - hetero_res$Q_df)/hetero_res$Q  # I2 = 100%×(Q - df)/Q     (here show between 0 and 1 so no multiplication by 100%) 
          hetero_res$i_squared <- ifelse(hetero_res$i_squared < 0, 0, hetero_res$i_squared)  # negative values set to 0
          hetero_name <- paste0(output_dir,"hetero/", filename, ".hetero.txt")
          write_tsv(hetero_res, file=hetero_name)
        }
        
        # steiger
        # note exposure's smaple size is present but outcome sample size is not
        # exp_dat_outcome$samplesize.outcome <- ith_outcome_path[1,]$sample_size

        steiger <- directionality_test(exp_dat_outcome)
        steiger_name <- paste0(output_dir, "steiger/", filename, ".steiger.txt")
        write_tsv(steiger, file=steiger_name)
        
        cat(filename, ': done\n\n')
        
      } # close if-else
    }
  }, error=function(e){cat("ERROR :", conditionMessage(e), "\n\n")}) # tryCatch({ # try catch outcome error
}  # close if

end_time <- Sys.time()

print(paste0("Time taken: ", end_time - start_time))