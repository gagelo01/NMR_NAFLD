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
res_ivw_nopleiotropicregion <- fread("Data/Modified/res_ivw_nopleiotropicregion")
list_table <- readRDS("Data/Modified/list_tablesexspecific.rds")
list_table <- list(rbindlist(list_table[1:3]), rbindlist(list_table[4:6]))
list_table <- lapply(list_table, function(table) {
  table[, variable :=  gsub("model", "Hazard ratios for incident NAFLD",variable)]
  table[is.na(Sex), Sex := ""]
  setNames(table, c("Sex", "variable", paste0("quintile", 1:5)))
})

table4 <- fread("Data/Modified/tablecombineddyslipidemia.txt")
top <- tribble(
  ~colA, ~colB, ~colC, ~colD, ~colE, ~colF,
  "Sex",   "HDL cholesterol levels", "High", "High", "Low", "Low",
  "",   "Triglyceride levels", "Low", "High", "Low", "High",
)
table4 <- rbind(top, setNames(table4, names(top)))

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
coxHR<-coxHR[exposure!="GLOBAL",][cov_inc== "+ age_enrollment + eth + sex + med + smoking + alcohol + townsend + WC",]
coxHR[,IV:=NULL]
coxHR[,outcome := gsub("nafld", "NAFLD", outcome)]
toremove<-c("ph_chisq", "ph_df", "ph_p")
coxHR[,(toremove):=NULL]
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


obj<-c("coxHR", "inst", "FandQ", "res_univariate", "egger_intercept", "mvmr_results","res_ivw_nopleiotropicregion")
for(i in 1:length(obj)) {
  assign(obj[i], cleanify(dat = get(obj[i]), trad = trad))
}

##### Supllementary table description
dt_title <- data.table(title = paste0("ST", c(1:10)),
                       caption = c( "Cox Hazard ratio results for all 170 metabolic factors on NAFLD adjusted for age enrollment, ethnicity, sex, medicine, smoking, alcohol, Townsend deprivation index, waist circumference",
                                    "Description of the datasets used for Two-Sample Mendelian randomization.",
                                   "Instruments and relevant statistics. Rsid, chromosome, position, beta, standard error, for each exposure and each SNPs",
                                   "Instrument strength (F-statistics) and heterogeneity statistics (Cochran's Q) for univariable MR",
                                   "Univariable Mendelian Randomization results using different pleiotropy robust MR analyses",
                                   "Univariable Mendelian Randomization Egger's intercept",
                                   "Event rates and hazard ratios for non-alcoholic fatty liver disease in UK Biobank participants for the entire cohort are presented, as well as in men and women separately. Hazard ratios for incident NAFLD are adjusted for age, ethnicity, smoking status (ever smoked/never smoked), alcohol intake, Townsend's deprivation index, waist circumference as well as medication for cholesterol, type 1 diabetes, and hypertension.",
                                   "Event rates and hazard ratios for non-alcoholic fatty liver disease in UK Biobank participants for the entire cohort are presented, as well as in men and women separately. Hazard ratios for incident NAFLD are adjusted for age, ethnicity, smoking status (ever smoked/never smoked), alcohol intake, Townsend's deprivation index, waist circumference as well as medication for cholesterol, type 1 diabetes, and hypertension.",
                                   "Event rates and hazard ratios for non-alcoholic fatty liver disease in UK Biobank participants for the entire cohort are presented, as well as in men and women separately. Hazard ratios for incident NAFLD are adjusted for age, ethnicity, smoking status (ever smoked/never smoked), alcohol intake, Townsend's deprivation index, waist circumference as well as medication for cholesterol, type 1 diabetes, and hypertension.",
                                   "Multivariable Mendelian randomization results using different pleiotropy robust MVMR analyses",
                                   "IVW for all metabolic measures with removing GCKR region"),
                       legends = c("exposure = unique exposure acronym; exposure_long = exposure name without acronym; outcome = unique outcome acronym; outcome_long = outcome name without acronym; HR = Hazard ratio; lci = lower confidence interval; uci = upper confidence interval; pval = P-value; cov_inc = covariate included in the model.",
                                   "trait = Unique trait identifier; group_name = the name of the study group;  year = year data was published; author = author of the data; consortium = consortium for the sample; sex = sex of the sample (in this study always Males and Females); population = Ancestry of the sample (in this study always European); unit = The unit either standard deviation (SD) or log Odds ratio (log(orr)); nsnp = number of SNPs in the summary statistics; sample_size = The maximum sample size; ncase = The maximum number of cases; ncontrol = the maximum number of controls; pmid = Pubmed ID",
                                   "SNP = the rsid;  effect_allele.exposure = effect allele of the exposure; other_allele.exposure = non effect allele of the exposure effect_allele.outcome = effect allele of the outcome; other_allele.outcome = non effect allele of the outcome; beta.exposure =  SNP effect on exposure; beta.outcome = SNP effect on outcome; eaf.exposure = effect allele frequency in the exposure;    eaf.outcome = effect allele frequency in the outcome; remove = should you remove the SNP; palindromic = is the SNP palindromic; ambiguous = is the SNP A/T or C/G; outcome = Unique identifier for the outcome; chr.outcome = chromosome of the SNP integer from 1 to 22; pos.outcome = position on the chromosome (GRCH37); se.outcome = standard error of the SNP outcome association; pval.outcome = p-value of the SNP outcome association; samplesize.outcome = Maximum sample size of the outcome GWAS; ncase.outcome = Number of case of the outcome GWAS;        ncontrol.outcome = Number of controls of the outcome GWAS; mr_keep.outcome = Should you remove this genetic instruments because the allele frequency does not match;",
                                   "outcome = Unique name for the outcome; exposure = Unique name for the exposure; method = the method to compute the cochran's Q test; Q = cochran's Q value; Q_df = The cochran's Q statistics number of degree of freedom; Q_pval = The Cochran's Q  p-value; fstat = the F-statistics of the instrumental variables; rsq = the variance explained of the exposure by the instrumental variables.",
                                   "outcome = unique name of the outcome; exposure = unique name of the exposure; method = name of the method to estimate the causal effect of the exposure on the outcome; nsnp = The number of SNPs genetic instruments used to compute the estimate; b = the estimate scaled per SD or log(orr);  se =  standard error of the estimate; pval = the p-value of the estimate; lci = the 95% confidence intervall lower bound of the effect; uci = the 95% confidence intervall upper bound of the effect; type_of_test = categorisation of the method based on Slob and Burgess pmid = 32249995.",
                                   "outcome = unique name for the outcome; exposure = unique name for the exposure; egger_intercept = the egger intecept estimate; se = the egger intecept standard error; pval = the egger intecept p-value.",
                                   "Sex = sex; variable = different value; quintile = 5 quintiles",
                                   "Sex = sex; variable = different value; quintile = 5 quintiles",
                                   "Different value of triglyceride and cholesterol and different hazard ratio.",
                                   "exposure = unique name for the exposure; outcome = unique name for the outcome; correctfor = A unique id for the variable that was included in MVMR to adust for the exposure; correctfor_long = The variable correctfor, but with a longer more comprehensive name  b = the direct effect (adjusted for the other exposures); se = the standard error of the direct effect; lci = The 95% confidence intervall lower bound of the direct effect; uci = The 95% confidence intervall upper bound of the direct effect; pval = The p-value of the direct effect; cochranQ = The cochran's Q statitics; cochranQpval =The cochran's Q statitics  p-value; nsnp = the number of SNPs (genetic instruments); method = The name of the multivarialble MR method used;   F_stastistics = The conditionnal F-statistics.",
                                   "outcome = unique name of the outcome; exposure = unique name of the exposure; method = name of the method to estimate the causal effect of the exposure on the outcome; nsnp = The number of SNPs genetic instruments used to compute the estimate; b = the estimate scaled per SD or log(orr);  se =  standard error of the estimate; pval = the p-value of the estimate."
                       ))

#
list_supdat <- list("Tables captions and titles" = dt_title,
                    "ST1" = coxHR,
                    "ST2" = dataset,
                    "ST3" = inst,
                    "ST4" = FandQ,
                    "ST5" = res_univariate,
                    "ST6" = egger_intercept,
                    "ST7" = list_table[[1]],
                    "ST8" = list_table[[2]],
                    "ST9" = table4,
                    "ST10" = mvmr_results,
                    "ST11" = res_ivw_nopleiotropicregion)
for(i in 1:length(list_supdat)) {
  writexl::write_xlsx(x = list_supdat[[i]],
                      path = paste0("Results/", gsub(" ", "_", names(list_supdat)[[i]]), ".xlsx"))
}
writexl::write_xlsx(x = list_supdat,
                    path = "Results/supplementary_tables_clean.xlsx")

message("This script finished without errors")


