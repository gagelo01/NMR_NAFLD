#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(writexl)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/Triglycerides_dis")

all_inst_mvmr <- fread( "Data/Modified/all_inst_mvmr.txt")
all_outcome_mvmr <- fread( "Data/Modified/all_outcome_mvmr.txt")
inst_all_sign_clump <- fread( "Data/Modified/inst_all_sign_clump.txt")
res_univariate <- fread( "Data/Modified/res_univariate.txt")
corrmat<-readRDS("Data/Modified/correlationmatrix.rds")
resmvmr <- readRDS( "Data/Modified/res_mvmr.rds")
dat<-fread("Data/Modified/observationalfulldata.txt")
final <- fread("Data/Modified/sample_summary_table.txt")
FandQ <- fread( "Data/Modified/FandQ_univariate.txt")
egger_intercept <- fread("Data/Modified/egger_intercept.txt")
coxHR<-fread( "Data/Modified/coxHR.txt")
trad <- fread( "Data/Modified/trad")
dictionnary<-tribble(
  ~id,  ~trait,  ~exposure,
  "Fat_Liver", "Liver Fat", "Fat_Liver",
  "NAFLD",   "Non-alcoholic fatty liver disease", "NAFLD",
  "UKB-b-9405",   "Waist circumference", "UKB-b-9405",)
trad<-rbind(trad, dictionnary)

###################
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao_small <- ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/"), ]

#####################
#Table 1
writexl::write_xlsx(x = final,
                    col_names = FALSE,
                    path = "Results/Table1.xlsx")

##Supplementary Table
#supptable 1 description of cohorts
idserver_exposure <- c("dis-2-1", "trait-14-8")
idmrbase_exposure <- c("ukb-b-9405", ao_small[author == "Borges CM", id][1])
dataset <- rbindlist(list(ao[id %in% idmrbase_exposure,], df_index[id %in% idserver_exposure,]), fill = TRUE)
dataset <- dataset[,.(trait,group_name, year, author,  consortium, sex, population, unit, nsnp, sample_size,ncase, ncontrol,pmid)]
dataset[author=="Borges CM", pmid := 35692035]
dataset[author=="Borges CM", consortium := "UKBiobank"]
dataset[author=="Borges CM", unit :="SD"]
dataset[, trait := trait %>% gsub("Glutamine", "NMR_Metabolites", . ) %>% sub("Fat_Liver", "Liver_Fat", .)]

#
coxHR<-coxHR[exposure!="GLOBAL",][cov_inc== "+ age_enrollment + eth + sex + med + smoking + alcohol + townsend",]
coxHR[,IV:=NULL]
coxHR[,outcome := gsub("nafld", "NAFLD", outcome)]
#
inst <- inst_all_sign_clump
inst[,c("id.outcome", "pval_origin.outcome", "action", "pval_origin.exposure", "id.exposure", "exposure_outcome", "mr_keep") := NULL]
# inst <- merge(ao_small[,.(id, trait)], inst, by.x = "id", by.y = "exposure", all.y = TRUE)

egger_intercept[,c("id.exposure","id.outcome") := NULL] 

#
mvmr_results <- lapply(resmvmr, function(x)
  x[, otherexposure := apply(.SD, 1,function(x) paste(setdiff(unique(exposure), x), collapse = "+")), .SDcols = "exposure"])
mvmr_results<-rbindlist(mvmr_results)
mvmr_results <- mvmr_results[outcome == "NAFLD"][correctfor == "UKB-b-9405",]
mvmr_results[, c("clump_exposure", "otherexposure") := NULL]
##########
cleanify <- function(dat, argu = c("correctfor", "outcome", "exposure"), trad) {
  for(i in 1:length(argu)) {
    if(argu[i] %in% colnames(dat)) {
      dat[, (argu[i]) :=gsub("^met-d-", "", get(argu[i]))] # 
      dat <- merge(trad[,.(id, trait)], dat, by.x = "id", by.y = argu[i])
      setnames(dat, c("id", "trait"), c(argu[i], paste0(argu[i], "_long")))
    }
  }
  return(dat)
}


obj<-c("coxHR", "inst", "FandQ", "res_univariate", "egger_intercept", "mvmr_results")
for(i in 1:length(obj)) {
  assign(obj[i], cleanify(dat = get(obj[i]), trad = trad))
}

#
dt_title <- data.table(title = paste0("Supplementary Table ", 1:7),
                       caption = c( "cox Hazard ratio results",
                                    "Description of the datasets used.",
                                   "Instruments and relevant statistics",
                                   "Instrument strength and heterogeneity statistics for univariable MR",
                                   "Univariable Mendelian Randomization results",
                                   "Univariable Mendelian Randomization Egger's intercept",
                                   "Multivariable Mendelian randomization results."))
                       
#                       
list_supdat <- list("Tables captions and titles" = dt_title,
                    "Supplementary Table 1" = coxHR,
                    "Supplementary Table 2" = dataset,
                    "Supplementary Table 3" = inst,
                    "Supplementary Table 4" = FandQ,
                    "Supplementary Table 5" = res_univariate, 
                    "Supplementary Table 6" = egger_intercept, 
                    "Supplementary Table 7" = mvmr_results)
for(i in 1:length(list_supdat)) {
  writexl::write_xlsx(x = list_supdat[[i]],
                      path = paste0("Results/", gsub(" ", "_", names(list_supdat)[[i]]), ".xlsx"))
}
writexl::write_xlsx(x = list_supdat,
                    path = "Results/supplementary_tables_clean.xlsx")

message("This script finished without errors")



