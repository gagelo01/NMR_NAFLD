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


####### 
k <- gsub("met-d-", "", causals)
otter_dendro <- as.dendrogram(hclust(d = dist(x = corremat[k, k])))
tree <- dendextend::cutree(otter_dendro,k=4)
tree<-factor(tree, labels = c("acetate", "hdl", "triglyceride", "pufa"))

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
format(round(datobs[toinclude==1,median(nafld_time)/365.25], digits = 2), nsmall = 2)

#para 2
trad[,.N]
k<-rescox[cov_inc == "+ age_enrollment + eth + sex + med + smoking + alcohol + townsend + WC"  & 
         exposure %in% trad$id,]
format(round(k[, c(min(HR), max(HR))], digits = 2), nsmall = 2)
k[pval<0.05, .N]
k[pval<0.05, sum(HR<1)]
k[pval<0.05, sum(HR>1)]

#######UVMR
#para 3
FandQ[outcome == "NAFLD", min(fstat)]

#para 4
nmrcoxmr[typefig == "Metabolites",  length(unique(exposure))]
return_format_fstat(nsnpfstat[exposure == "met-d-Acetate", ])
return_format_data(res_univariate[exposure %in% "met-d-Acetate" & method == "Inverse variance weighted" & outcome == "NAFLD",])
egger_intercept[outcome == "NAFLD" & exposure == "met-d-Acetate"]

#para 5
nmrcoxmr[typefig == "Lipids",  length(unique(exposure))]
return_format_data_noexp(rescox[cov_inc == "+ age_enrollment + eth + sex + med + smoking + alcohol + townsend + WC" & exposure == "PUFA_pct", ])
return_format_fstat(nsnpfstat[exposure == "met-d-PUFA_pct", ])
return_format_data(res_univariate[exposure %in% "met-d-PUFA_pct" & method == "Inverse variance weighted" & outcome == "NAFLD",])
return_format_data_noexp(rescox[exposure %in% "Omega_6_by_Omega_3",])

#para 6
dattg <- dat[!(!is.na(nafld_date) & f.53.0.0 > nafld_date | is.na(tg)), ]
dattg[,table(nafld_censored)]
return_format_data_noexp(rescoxtg[exposure != "GLOBAL" &exposure == "tg",])
return_format_data_noexp(rescoxtg[exposure != "GLOBAL" &exposure == "tg_quantile5",])

rescox[cov_inc %in% "+ age_enrollment + eth + sex + med + smoking + alcohol + townsend + WC",]
nmrcoxmr[typefig %in% c("Lipoproteins") & category != "",  length(unique(exposure))]
return_format_data(res_univariate[exposure %in% "met-d-Total_TG" & method == "Inverse variance weighted" & outcome == "NAFLD",])

k<-res_univariable[exposure=="Fat_Liver" & outcome %in% causals,][method == "Inverse variance weighted", ]
k[pval < 0.05, ]

#####Multivariable MR
#para 7
list_exp_multi$`UKB-b-9405` %>% length(.)
k <- dcast(nmrcoxmr[panel %in% c("UVMR", "MVMR with WC"), ], exposure ~ panel, value.var = c("b", "pval"))
k[, ratio := `b_MVMR with WC` / b_UVMR ] 
1-2*(1-pnorm(q=1,mean = 0, sd = 1))
k[pval_UVMR<0.05, paste(round(c(mean(ratio), sd(ratio)), digits = 2))]
trad[,exposure:=gsub("met-d-", "", id)]
k <- merge(k, trad, by = "exposure") 
k[pval_UVMR<0.05,][order(ratio),][c(1,.N),]
k[exposure=="Acetate"]      

#para 8
length(causals)
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

##########Method  ############
nmrcoxmr[, length(unique(exposure)), by = typefig]

