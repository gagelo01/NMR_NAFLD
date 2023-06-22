#!/usr/bin/env Rscript
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(GagnonMR)
library(writexl)
library(furrr)

setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Triglycerides_dis/")
all_inst_mvmr <- fread( "Data/Modified/all_inst_mvmr.txt")
all_outcome_mvmr <- fread( "Data/Modified/all_outcome_mvmr.txt" )
inst_all_sign_clump <- fread( "Data/Modified/inst_all_sign_clump.txt")

##############
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao_small <- ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/"), ]
ao_small[, tomirge := tolower(id)]
options(future.globals.maxSize= 5e9)
plan(multicore, workers = 20, gc = TRUE)#plan(sequential)
#######################
#change only this section
#univariable
idserver_exposure <- c("trait-14-8") #paste0("trait-16-", c(1,2,4))
idserver_outcome <- c("trait-14-8", "dis-2-1")
idmrbase_exposure <- c("ukb-b-9405", "ieu-a-61", ao_small[pmid ==32203549,]$id, ao_small[author == "Borges CM", id])
idmrbase_outcome <- NULL
#multivariable
u<-c("logOR", "log odds", "log orr")
k <- c(df_index[unit %in% u,trait],  ao_small[unit %in% u, id]) #remove binary as exposure
k <-c(df_index[id %in% idserver_exposure & !(trait %in% k), trait],
      all_inst_mvmr[tolower(exposure) %in% tolower(idmrbase_exposure) & !(exposure %in% k), unique(exposure)])
exposurefor<-NULL #The phenotype you want to correct for , but will be added as "exposure" which means there instruments will be selected
correctfor <- list("UKB-b-9405")#The phenotype you want to correct for , but will be added as "correctfor" which means there instruments won't selected
if(!is.null(exposurefor)){k <- purrr::cross2(k,exposurefor)}
k<-lapply(k, unlist)
mvmr_exposure <- k[sapply(k, function(x) (length(unique(x))==length(x)))]


mvmr_exposure <- unique(mvmr_exposure)

mvmr_outcome <- c("NAFLD", "Fat_Liver")
pval <- c(5e-8)
if(is.null(correctfor)){correctfor<-"NULL"}
arguments_mvmr <- purrr::cross(list(mvmr_exposure, mvmr_outcome, correctfor, pval))
test<- sapply(arguments_mvmr, function(x) !any(grepl(paste(x[[3]], collapse = "|"), x[[1]])))
arguments_mvmr <- arguments_mvmr[test]
arguments_mvmr <- c(arguments_mvmr)
###############
k <- c(df_index[unit %in% u,trait],  ao_small[unit %in% u, id])
exp_vec<-c(all_inst_mvmr[tolower(exposure) %in% tolower(idmrbase_exposure), unique(exposure)],
           df_index[id %in% idserver_exposure, trait])
out_vec<-c(all_inst_mvmr[tolower(exposure) %in% tolower(idmrbase_outcome), unique(exposure)],
           df_index[id %in% idserver_outcome, trait])
arguments_uni <- rbind(tidyr::crossing(exposure = exp_vec[!(exp_vec%in%k)], outcome = out_vec),
                       tidyr::crossing(exposure = out_vec[!(out_vec%in%k)], outcome = exp_vec))

arguments_uni <- distinct(arguments_uni)
setDT(arguments_uni)
arguments_uni <- arguments_uni[!(exposure == outcome),]
harm_univariate <- map(split(arguments_uni, 1:nrow(arguments_uni)), function(x) {
  harm <- TwoSampleMR::harmonise_data(exposure_dat = inst_all_sign_clump[exposure == x[,exposure], ],
                                      outcome_dat = all_outcome_mvmr[outcome == x[,outcome], ],
                                      action = 1)
  return(harm)}) %>% rbindlist(.,fill = TRUE)


setDT(harm_univariate)
fwrite(harm_univariate, "Data/Modified/harm_univariate.txt")
harm_univariate[, exposure_outcome :=paste0(exposure, "_", outcome)]

harm_univariate <- TwoSampleMR::steiger_filtering(harm_univariate) %>% as.data.table(.)
harm_univariate <- harm_univariate[mr_keep==TRUE,]
fwrite(harm_univariate, "Data/Modified/harm_univariate.txt")
harm_univariate <- harm_univariate[steiger_dir == TRUE,]
egger_intercept <- mr_pleiotropy_test(harm_univariate)
harm_univariate <- harm_univariate[!(exposure == "Fasting_Insulin" & SNP == "rs1121980"),] #in the FTO region
list_harm_univariate <- split(harm_univariate, harm_univariate$exposure_outcome)

list_res_univariate <- future_map(list_harm_univariate,
                                  function(x) {
                                    GagnonMR::all_mr_methods(x,
                                                             skip_presso = FALSE)},
                                  .options = furrr_options(seed = TRUE))

res_univariate <- rbindlist(list_res_univariate, fill = TRUE)
res_univariate[, c("id.exposure", "id.outcome") := NULL]

FandQ <- lapply(list_harm_univariate, function(x) {
  res <- TwoSampleMR::mr_heterogeneity(x)
  if(nrow(res)==0){res<-data.frame(exposure = x$exposure, outcome = x$outcome)}
  x <- TwoSampleMR::add_rsq(x)
  res$fstat<-GagnonMR::fstat_fromdat(list(x))
  res$rsq <- sum(x$rsq.exposure)
  return(res)
}) %>% rbindlist(.,fill = TRUE)


res_ivw_nopleiotropicregion <- TwoSampleMR::mr(harm_univariate[pleiotropic_region != "GCKR",], method_list = "mr_ivw")
fwrite(res_ivw_nopleiotropicregion, "Data/Modified/res_ivw_nopleiotropicregion")
fwrite(FandQ, "Data/Modified/FandQ_univariate.txt")
fwrite(res_univariate, "Data/Modified/res_univariate.txt")
setDT(egger_intercept)
fwrite(egger_intercept, "Data/Modified/egger_intercept.txt")

#mvmr
performmvmr <- function(exposure_vec, outcome_vec, correctfor=NULL,pval_threshold = 1,
                        clump_exp_arg = "none", pval_inst = 5e-8,
                        clump_r2 = 0.01, clump_kb = 1000) { #clump_exp_arg either none, first, or second
  message(paste0("MVMR for the effect of ", paste(exposure_vec, collapse = " + "), " on ", outcome_vec, " while correcting for ",  paste(correctfor, collapse = " + ")))
  exposure_dat <- inst_all_sign_clump[inst_all_sign_clump$exposure %in% exposure_vec & inst_all_sign_clump$pval.exposure < pval_inst,]
  k<-all_inst_mvmr[(all_inst_mvmr$exposure %in% correctfor) & (all_inst_mvmr$SNP %in% unique(exposure_dat$SNP)),]
  exposure_dat<- rbind(exposure_dat,k, fill = TRUE)
  d1 <- all_inst_mvmr[(all_inst_mvmr$exposure %in% c(exposure_vec,correctfor)) & (all_inst_mvmr$SNP %in% unique(exposure_dat$SNP)),]

  if(clump_exp_arg == "none") {clump_exp<-NULL} else if(clump_exp_arg == "first"){clump_exp<-exposure_vec[1]} else if(clump_exp_arg == "second"){clump_exp<-exposure_vec[2]}
  inst_mvmr <- prepare_for_mvmr(exposure_dat = exposure_dat, d1 =d1,clump_r2 = clump_r2, clump_kb = clump_kb, pval_threshold = pval_threshold, clump_exp = clump_exp, harmonise_strictness = 1)

  exposure_outcome_harmonized <- TwoSampleMR::mv_harmonise_data(exposure_dat = inst_mvmr,
                                                                outcome_dat = all_outcome_mvmr[all_outcome_mvmr$outcome == outcome_vec,],
                                                                harmonise_strictness = 1)
  mvmr_results <- GagnonMR::mv_multiple_MendelianRandomization(exposure_outcome_harmonized = exposure_outcome_harmonized)
  mvmr_results$clump_exposure <- mvmr_results$clump_exp %>% ifelse(is.null(.), "none", .)
  mvmr_results <- mvmr_results[!(mvmr_results$exposure %in% correctfor),]
  mvmr_results$correctfor <- paste(correctfor, collapse = " + ")
  return(mvmr_results)
}

resmvmr <- future_map(arguments_mvmr, function(x) {
  resnone<-performmvmr(exposure_vec = x[[1]],
                       outcome_vec = x[[2]],
                       correctfor=if(x[[3]]=="NULL"){NULL}else{x[[3]]},
                       clump_exp_arg = "none",
                       pval_inst = x[[4]])
  return(resnone)
}, .options = furrr_options(seed = TRUE))

names(resmvmr) <- sapply(arguments_mvmr, function(x) paste0(paste(x[[1]],collapse = " + "), " ~ ", x[[2]], " correctfor = ", x[[3]]," (pval=", x[[4]],")"))

saveRDS(resmvmr, "Data/Modified/res_mvmr.rds")

message("This script finished without errors")

