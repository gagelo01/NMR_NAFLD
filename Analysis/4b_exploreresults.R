#!/usr/bin/env Rscript
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(GagnonMR)
library(VennDiagram)


setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Triglycerides_dis/")
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao_small <- ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/"), ]

#########
resmvmr <- readRDS( "Data/Modified/res_mvmr.rds")
res_univariate <- fread("Data/Modified/res_univariate.txt")
res_univariate <- res_univariable[!(method  %in% c("Robust adjusted profile score (RAPS)", "Weighted mode")), ]
rescox <- fread( "Data/Modified/coxHR.txt")
corrmat<-readRDS("Data/Modified/correlationmatrix.rds")
egger_intercept <- fread( "Data/Modified/egger_intercept.txt")

#####
res_univariate<-merge(res_univariate, ao_small[, .(id, trait)], by.x = "exposure", by.y = "id", all.x = TRUE)
res_univariate <- res_univariate[!(method %in% c("MR Egger", "Weighted mode","Robust adjusted profile score (RAPS)")),]
mvmr_results <- lapply(resmvmr, function(x)
  x[, correctedfor := apply(.SD, 1,function(x) paste(setdiff(unique(exposure), x), collapse = "+")), .SDcols = "exposure"])
mvmr_results<-rbindlist(mvmr_results)
mvmr_results<-merge(mvmr_results, ao_small[, .(id, trait)], by.x = "exposure", by.y = "id", all.x = TRUE)

#######
trad <- ao_small[author == "Borges CM", .(id,trait)]
trad<- trad[!grepl("ratio ", tolower(trait)) | id %in% c("met-d-PUFA_pct", "met-d-Omega_6_by_Omega_3"),]
trad[,exposure:=id]
trad[,id := id %>% gsub("met-d-", "", .)]
fwrite(trad, "Data/Modified/trad")
nonratioid<-trad$id 

#####
expunisign <- res_univariate[exposure %in% c(ao_small[author == "Borges CM", id]) & outcome == "NAFLD", ][, all(pval<0.05) & (all(b>0)|all(b<0)), by = "exposure"][V1==TRUE]$exposure
expunisign <- expunisign[expunisign %in% paste0("met-d-",nonratioid)]
expunisign <- expunisign[expunisign %in%  egger_intercept[outcome == "NAFLD", ][pval > 0.05,]$exposure]

#####
corect <- c("logTG_GLGC_2022", "UKB-b-9405", "HDL_GLGC_2022") 
list_exp_multi <- vector(mode = "list", length = length(corect))
names(list_exp_multi)<- corect
for(i in 1:length(corect)) {
expmultisign <- mvmr_results[correctfor %in% corect[i] & outcome == "NAFLD" & !(method %in% c("Multivariable Egger Intercept", "Multivariable Egger"))][, all(pval<0.05) & (all(b>0)|all(b<0)), by = "exposure"][V1==TRUE]$exposure
k <- mvmr_results[correctfor %in% corect[i] & outcome == "NAFLD" & method == "Multivariable Egger Intercept" & pval > 0.05,]$exposure
expmultisign <- expmultisign[expmultisign%in%k]
expmultisign <- expmultisign[expmultisign %in% paste0("met-d-",nonratioid)]
list_exp_multi[[i]] <- expmultisign
}

causals <- intersect(expunisign, list_exp_multi$`UKB-b-9405`)
saveRDS(causals, "Data/Modified/causals.rds")
saveRDS(list_exp_multi, "Data/Modified/list_exp_multi.rds")
saveRDS(expunisign, "Data/Modified/expunisign.rds")
#####
hdl <- intersect(causals, list_exp_multi$logTG_GLGC_2022)
trad[id %in% gsub("met-d-","", hdl),]
res_univariate[exposure %in% hdl,]

tg <- intersect(causals, list_exp_multi$HDL_GLGC_2022)



###############
mvmr_results[correctfor %in% "logTG_GLGC_2022" & outcome == "NAFLD" & exposure %in% causals, ]
res_univariate[outcome == "NAFLD" & exposure == "met-d-L_VLDL_TG", ]

