#!/usr/bin/env Rscript
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(GagnonMR)
library(dendextend)

setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Triglycerides_dis/")
ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao_small <- ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/"), ]
trad <- fread( "Data/Modified/trad")
resmvmr <- readRDS( "Data/Modified/res_mvmr.rds")
res_univariate <- fread("Data/Modified/res_univariate.txt")
res_univariate <- res_univariate[!(method  %in% c("MR Egger", "Robust adjusted profile score (RAPS)", "Weighted mode")), ]
rescox <- fread( "Data/Modified/coxHR.txt")
rescoxtg <- fread("Data/Modified/rescoxtg.txt")

corrmat<-readRDS("Data/Modified/correlationmatrix.rds")
causals<- readRDS("Data/Modified/causals.rds")
list_exp_multi <- readRDS("Data/Modified/list_exp_multi.rds")
expunisign <- readRDS("Data/Modified/expunisign.rds")
datobs <- fread("Data/Modified/observationalfulldata.txt")
table1<- fread("Data/Modified/sample_summary_table.txt")
nmrcoxmr <- fread("Data/Modified/nmrcoxmr.txt")
list_exp_multi <- readRDS( "Data/Modified/list_exp_multi.rds")
corremat <- readRDS( "Data/Modified/correlationmatrix.rds")
FandQ<-fread("Data/Modified/FandQ_univariate.txt")
harm_univariate <- fread( "Data/Modified/harm_univariate.txt")
nsnpfstat <- merge(FandQ[method == "Inverse variance weighted" & outcome == "NAFLD",],
                   harm_univariate[outcome == "NAFLD", .N, by = "exposure"],
                   by = "exposure")
egger_intercept <- fread( "Data/Modified/egger_intercept.txt")
caseevent <- fread( "Data/Modified/tablecombineddyslipidemia.txt")


#######
k <- gsub("met-d-", "", causals)
otter_dendro <- as.dendrogram(hclust(d = dist(x = corremat[k, k])))
tree <- dendextend::cutree(otter_dendro,k=3)
tree<-factor(tree, labels = c("hdl", "triglyceride", "pufa"))

######
return_format_data<-function(data) {
  return(data[, paste0(format(round(exp(b), digits = 2), nsmall = 2), " 95% CI=", format(round(exp(lci), digits = 2), nsmall = 2), "-",  format(round(exp(uci), digits = 2), nsmall = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))])
}
return_format_data_noexp <-function(data) {
  return(data[, paste0(format(round(HR, digits = 2), nsmall = 2) , " 95% CI=", format(round(lci, digits = 2), nsmall = 2), "-",  format(round(uci, digits = 2), nsmall = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))])
}
return_format_fstat <-function(data) {
  return(data[, paste0(N, " SNPs (r2 = ", round(rsq*100, digits =2), "%; F-statistics = ",  round(fstat, digits = 0), ")")])
}

##Abstract
nmrcoxmr[, length(unique(exposure)), by = typefig]
table1[V1 == "N",as.numeric(V4)+as.numeric(V5)]
datobs[toinclude==1,median(nafld_time)/365.25]
table1[V1 == "N",as.numeric(V4)]
length(causals)
nmrcoxmr[panel %in% c("cox", "UVMR") & exposure %in% gsub("met-d-", "", causals), all(pval<0.05) & (all(b>0)|all(b<0)), by = "exposure"][V1==TRUE, .N]
table(tree)
k<-paste0("met-d-", c("Acetate", "L_HDL_P", "Total_TG", "PUFA_pct"))
return_format_data(res_univariate[exposure %in% k & method == "Inverse variance weighted" & outcome == "NAFLD",])
# list_mean_HR<-vector(mode = "list", length = length(unique(tree)))
# for(i in unique(tree)) {
# list_mean_HR[[i]] <- nmrcoxmr[panel %in% "UVMR", ][exposure %in% names(tree)[tree==i], paste(round(c(mean(HR),mean(exp(lci)), mean(exp(uci))), digits = 2), collapse = " ")]
# }
# list_mean_HR


#######Intro
nmrcoxmr[, length(unique(exposure)), by = typefig]

########Results
nmrcoxmr[, length(unique(exposure)), by = typefig]

#########cox
##Para 1
datobs[NMR==1,.N]
datobs[NMR==1][toinclude==0,.N]
datobs[toinclude==1, sum(nafld_censored)]
table1[V1== "N", V4]
table1[V1== "Age at enrollment, years (SD)", V4]
table1[V1== "Male, n (%)", V4]
table1[V1== "N", V5]
table1[V1== "Age at enrollment, years (SD)", V5]
table1[V1== "Male, n (%)", V5]
format(round(datobs[toinclude==1,median(nafld_time)/365.25], digits = 1), nsmall = 1)
datobs[nafld_censored==1 & toinclude ==1, .N]
#para 2
trad[,.N]
k<-rescox[cov_inc == "+ age_enrollment + eth + sex + med + smoking + alcohol + townsend + WC"  &
         exposure %in% trad$id,]
format(round(k[, c(min(HR), max(HR))], digits = 2), nsmall = 2)
k[pval<0.05/ntest, .N]
k[pval<0.05/ntest, sum(HR<1)]
k[pval<0.05/ntest, sum(HR>1)]

#######UVMR
#para 3
FandQ[outcome == "NAFLD", min(fstat)]
GagnonMR::from_genecard_to_generegion("GCKR", window = 1e6)
#para 4
nmrcoxmr[typefig == "Metabolites",  length(unique(exposure))]
return_format_fstat(nsnpfstat[exposure == "met-d-Acetate", ])
return_format_data(res_univariate[exposure %in% "met-d-Acetate" & method == "Inverse variance weighted" & outcome == "NAFLD",])
egger_intercept[outcome == "NAFLD" & exposure == "met-d-Acetate"]
res_univariate[outcome == "NAFLD" & exposure == "met-d-Acetate",]

#para 5
nmrcoxmr[typefig == "Lipids",  length(unique(exposure))]
return_format_data_noexp(rescox[cov_inc == "+ age_enrollment + eth + sex + med + smoking + alcohol + townsend + WC" & exposure == "PUFA_pct", ])
return_format_fstat(nsnpfstat[exposure == "met-d-PUFA_pct", ])
return_format_data(res_univariate[exposure %in% "met-d-PUFA_pct" & method == "Inverse variance weighted" & outcome == "NAFLD",])
return_format_data_noexp(rescox[exposure %in% "Omega_6_by_Omega_3",])

#para 6-7
dattg <- dat[!(!is.na(nafld_date) & f.53.0.0 > nafld_date), ]
dattg[!(is.na(hdl) & is.na(tg)),table(nafld_censored)]
return_format_data_noexp(rescoxtg[exposure != "GLOBAL" &exposure == "tg",])
return_format_data_noexp(rescoxtg[exposure != "GLOBAL" &exposure == "hdl",])

return_format_data_noexp(rescoxtg[exposure != "GLOBAL" &exposure == "tg_all_quantile5",])
return_format_data_noexp(rescoxtg[exposure != "GLOBAL" &exposure == "hdl_all_quantile5",])

#para 7
caseevent[sex=="All" & variable == "Hazard ratios for incident NAFLD",FALSE_FALSE]

#####Multivariable MR
#para 8
list_exp_multi$`UKB-b-9405` %>% length(.)
k <- dcast(nmrcoxmr[panel %in% c("Univariable MR", "Multivariable MR with WC") & outcome == "nafld", ], exposure ~ panel, value.var = c("b", "pval"))
k[, ratio := `b_Multivariable MR with WC` / `b_Univariable MR` ]
1-2*(1-pnorm(q=1,mean = 0, sd = 1))
k[`pval_Univariable MR`<0.05, paste(round(c(mean(ratio), sd(ratio)), digits = 2))]
k <- merge(k, trad[,.(id, trait)], by.x = "exposure", by.y = "id")
k[`pval_Univariable MR`<0.05,][order(ratio),][c(1,.N),]
k[exposure=="Acetate"]

#para 8
length(causals)
harm_univariate[steiger_dir == FALSE & outcome == "NAFLD" & exposure %in% causals,]

clust2<-names(tree)[tree=="hdl"] #HDL-C
length(clust2)
format(round(min(corremat[clust2,clust2]), digits = 2), nsmall = 2)
clust3<-names(tree)[tree=="triglyceride"] #Triglycerides
length(clust3)
format(round(min(corremat[clust3,clust3]), digits = 2), nsmall = 2)

#########Discussion   ###########
nmrcoxmr[, length(unique(exposure)), by = typefig]
length(causals)
table1[V1== "N", V2]
table1[V1== "N", V3]

#limit
{nrows<-Inf
ukbpath<-"/home/couchr02/Mendel_UKB/Source/Phenotype/February_2023_Update/ukb671338.tab" #this datafield is not available in the release of march 2022
columnname <- fread(ukbpath, nrows = 2)
columnname <- names(columnname)
list_id<-NULL
fieldid <- c("f.40061")
data_disease <- fread(input = ukbpath,
                      select = c(1, which(grepl(paste(c(unlist(list_id), fieldid ), collapse = "|"), columnname))),
                      nrows = nrows)
data_disease <- data_disease[2:.N,]
data_disease[!is.na(f.40061.2.0), paste0(sum(f.40061.2.0>5), "/", length(f.40061.2.0), 
                                         " or ", round(sum(f.40061.2.0>5) / length(f.40061.2.0), digits = 2)*100, "%")]}

##########Method  ############
nmrcoxmr[, length(unique(exposure)), by = typefig]
short_func<-function(x) {
  x %>% str_split_fixed(.,pattern = fixed("("), n =2) %>% .[1,1] %>% as.numeric
}
table1[V1 == "Male, n (%)", (short_func(V4) + short_func(V5))] / table1[V1 == "N",as.numeric(V4) + as.numeric(V5)] 
table1[V1 == "Ethnicity (White), n (%)", (short_func(V4) + short_func(V5))] / table1[V1 == "N",as.numeric(V4) + as.numeric(V5)] 
datobs[toinclude==1,median(nafld_time)/365.25]
table1[V1 == "N",as.numeric(V4)]
######Figures legends#####
caseevent[sex=="All" & variable == "Cases/total (event rate, %)", lapply(.SD, function(x) {
  x %>% gsub("^.*/", "", .) %>% gsub(" .*$", "", .) %>% as.numeric(.)
}),.SDcols = c("TRUE_TRUE", "TRUE_FALSE", "FALSE_TRUE", "FALSE_FALSE")]
