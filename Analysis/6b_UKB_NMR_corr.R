  #!/usr/bin/env Rscript
  #import instrument
  library(data.table)
  library(tidyverse)
  library(furrr)
  library(survival)
  library(lubridate)
  library(survminer)
  setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Triglycerides_dis")
  
  nrows<-Inf#Inf
  ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
  ao_small <- ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/"), ]
  dat <- fread("Data/Modified/UKB_NMR_Data_QCed_biomarkers.txt", nrows = nrows)
  
  regname <- colnames(dat)[!(colnames(dat) %in% c("eid", "visit_index")) ]
  dat[,M_LDL_FC_by_CE:=gsub("Inf", "NA", M_LDL_FC_by_CE)%>%as.numeric(.)]
  dat[,c(regname) := lapply(.SD, function(x) (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)), .SDcols = regname] #standardise
  dat[,NMR:=1]
  dat<- dat[, .SD[1], by = "eid"] #I choose the first measurment for every participants. 1426 participants with two measurments. 3712 with only measurment at time 1. 116595 with only measurment at time 0
  ukbpath<-"/home/couchr02/Mendel_UKB/Source/Phenotype/March_2022_Update/ukb51070.tab"
  columnname <- fread(ukbpath, nrows = 2)
  columnname <- names(columnname)
  
  ##########Find exclusion field
  #Unique thing to change
  fieldidok<-c( "K740", "K742", "K758","K760","K769" )
  #######################
  tolocate<-paste(unique(substr(fieldidok, 1, 3)), collapse = "|")
  fieldid <- c("f.41280.0", "f.41270.0") #date, icd10
  data_disease <- fread(input = ukbpath, 
                        select = c(1, which(grepl(paste(fieldid, collapse = "|"), columnname))),
                        nrows = nrows)
  data_disease <- data_disease[2:.N,]
  
  measure.vars<-colnames(data_disease)[grepl("f.41280.0", colnames(data_disease))]
  data_disease[, (measure.vars):=lapply(.SD,as.IDate), .SDcols = measure.vars]
  
  index<-data_disease[,apply(.SD, 1, function(x) any(grepl(tolocate, x))), .SDcols = colnames(data_disease)[grepl("41270", colnames(data_disease))]]
  k <- data_disease[index, ]
  list_res<-vector(mode = "list", length = length(fieldid))
  for(i in 1:length(fieldid)) {
    measure.vars<-colnames(data_disease)[grepl(fieldid[i], colnames(data_disease))]
    long <- melt(k, id.vars = "f.eid", measure.vars = measure.vars)
    long <- separate(long, col = variable, c("todump1", "data_field", "todump2", "index"))
    long[,c("todump1","todump2"):=NULL]
    long[,data_field := paste0("f.", data_field)]
    wide<- dcast(long, f.eid + index ~ data_field, value.var = "value")
    list_res[[i]]<- wide
  }
  
  wide <- merge(list_res[[1]], list_res[[2]], by = c("f.eid", "index"))
  rm(list="list_res")
  wide <- wide[grepl(tolocate, f.41270), ]
  wide <- dcast(wide, f.eid ~ f.41270, value.var = "f.41280")
  
  colnom <- colnames(wide)[!grepl(paste(c("f.eid", fieldidok),collapse = "|"), colnames(wide))]
  index <- wide[,apply(.SD, 1, function(x) any(!is.na(x))) , .SDcols = colnom]
  fidtoexclude<- wide[index, unique(f.eid)]
fidtoexclude #Those are the ID with exclusion code in the hospital inpatient record
########################
#create the data_disease data.table
list_id <- NULL
fieldid <- c("f.34.0.0", "f.52.0.0", "f.53.0.0", "f.31.0.0",  "f.21001.0.0", "f.48.0.0", 
             "f.30870.0", "f.30760.0", "f.30780.0", #Triglycerides, HDL, LDL
             "f.2443.0.0", #diabetes diagnosted by doctors at recruitment
             "f.189.0.0", #townsed deprivation index at recruitment
             "f.6177.0", "f.6153.0", #Medication for cholestero, blood pressure or diabetes Male/Female
             "f.21000.0.0",
             "f.131666.0.0", "f.131668.0", "f.131670.0", #K74, K75,K76
             "f.20116.0.0", #smoking status
             "f.40000.0.0", #date of death
             "f.1558.0.0") #alcohol intake frequency
data_disease <- fread(input = ukbpath, 
                      select = c(1, which(grepl(paste(c(unlist(list_id), fieldid ), collapse = "|"), columnname))),
                      nrows = nrows)
data_disease <- data_disease[2:.N,]

setnames(data_disease, c("f.21001.0.0", "f.48.0.0", "f.21000.0.0", "f.31.0.0", "f.30870.0.0", "f.30760.0.0", "f.30780.0.0", "f.2443.0.0", "f.189.0.0",
                         "f.20116.0.0","f.1558.0.0", "f.40000.0.0"),
         c("BMI", "WC", "eth", "sex", "tg", "hdl", "ldl", "t2d", "townsend", "smoking", "alcohol", "date_death"))
data_disease[eth %in% c(-3, -1), eth := NA]
data_disease[!is.na(eth),eth := substr(eth,1,1)]
data_disease[, eth := factor(eth, levels = 1:6, labels = paste0("Ethnicity (", c("White", "Mixed", "Southeast_Asian", "Black", "Chinese", "Other"), ")"))] #1 = white
data_disease[t2d %in% c(-3, -1), t2d := NA]
{medication<- c("f.6177.0", "f.6153.0")
colnom <- colnames(data_disease)[grepl(paste(medication, collapse = "|"), colnames(data_disease))]
data_disease[, (colnom) := lapply(.SD, function(x) x %>% ifelse(. %in% 1:3, 1, .) %>%
                                    ifelse(.==-7, 0, .) %>%
                                    ifelse(.%in%c(-1,-3,4,5), NA, .)), .SDcols = colnom]
data_disease[,(colnom) := lapply(.SD, as.logical) ,.SDcols = colnom]
myfunc<-function(x){ifelse(all(is.na(x)), NA, any(x, na.rm = TRUE))}
data_disease[, med := apply(.SD, 1, myfunc), .SDcols = colnom]
data_disease[,(colnom):=NULL ]}
data_disease[, smoking := smoking %>% ifelse(. %in% 1:2, 1, .) %>% ifelse(.%in%c(-3), NA, .)]
data_disease[, alcohol :=  ifelse(alcohol %in% -3, NA, 6-alcohol)]

dat_corr <- merge(dat, data_disease[,.(f.eid,BMI, WC)], by.x = "eid", by.y = "f.eid", all.x = TRUE, all.y = FALSE )


#correlation matrix
tmp<-cor(dat_corr, use = "pairwise.complete.obs")
saveRDS(tmp, "Data/Modified/correlationmatrix.rds")
rm(list="dat_corr")
##coxph for NAFLD
dat <- merge(dat, data_disease, by.x = "eid", by.y = "f.eid", all.x = FALSE, all.y = TRUE )
rm(list="data_disease")
dat[, date_birth := ymd(paste0(f.34.0.0, "-", f.52.0.0, "-01")) %>% as.IDate(.)]
dat[,age_enrollment := f.53.0.0 - date_birth]
dat[, nafld_date := apply(.SD, 1, function(x) min(x, na.rm = TRUE)), .SDcols = c("f.131666.0.0","f.131668.0.0", "f.131670.0.0")]
dat[, nafld_date := nafld_date %>% as.IDate(.)]
dat[eid %in% fidtoexclude,nafld_date:=NA]#Since these participants have exclusion field they are tagged as not having NAFLD.
date_end_followup <- dat[, max(nafld_date, na.rm = TRUE)] 
date_end_followup <- round_date(date_end_followup, unit = "month")+30 #Here is the date the follow up ends.
dat[,censoring_date := dplyr::if_else(is.na(date_death), ymd(date_end_followup) %>% as.IDate(.), date_death) %>%
      dplyr::if_else(!is.na(nafld_date), nafld_date, .)] 
dat[,toinclude:=1]
dat[!is.na(nafld_date) & f.53.0.0 > nafld_date, toinclude:=0 ]#Remove every participants with a disease before enrolment
dat[, nafld_censored := ifelse(is.na(nafld_date), 0, 1)]
dat[, nafld_time := censoring_date - f.53.0.0]
dat[is.na(NMR),NMR:=0]
dat[NMR==0,toinclude:=0]
col <- c("med", "smoking", "sex"); dat[, (col):=lapply(.SD, as.logical), .SDcols = col]
fwrite(dat, "Data/Modified/observationalfulldata.txt")
# summary table

k<-dat[, c("nafld_censored","WC", "sex", "tg", "hdl", "ldl", "eth", #j'enlève ethnicité pcq c,est un facteur
           "townsend", "med", "smoking", "alcohol", "age_enrollment", "toinclude")]
k[,age_enrollment := age_enrollment/365.25]

getdt_logical <- function(variable_factor) {
  data<-data.table(variable_factor = variable_factor)
  data$dummy <- TRUE
  data$id <- 1:nrow(data)
  data_wide<- dcast(data, id ~ variable_factor, value.var = "dummy")
  setnames(data_wide, "NA", "non_available")
  dat <- data_wide[, lapply(.SD, function(x) ifelse(is.na(x), FALSE, TRUE)),
                   .SDcols = data[, sapply(.SD, levels),.SDcols = "variable_factor"]]
  dat <- dat[, lapply(.SD, function(x) ifelse(data_wide[, !is.na(non_available)],NA ,x))]
  return(dat)
}

tbl_summary <- function(x) {
  if(is.numeric(x)) {
    return(paste0(round(mean(x, na.rm = TRUE), digits = 1), " (", round(sd(x, na.rm = TRUE), digits = 1), ")"))
  }
  
  if(is.logical(x)) {
    toto<-sum(x, na.rm = TRUE)
    return(paste0(toto, " (", round(100*toto/length(x), digits = 1), ")"))
  }
  
  if(is.factor(x)){
    return(NA)}
}


k[, nafld_censored := ifelse(nafld_censored == 1, "Cases", "Controls")]

toinclude<-c("age_enrollment",  "sex", "WC", "tg", "hdl", "ldl", "townsend", "med",
             "smoking", "alcohol", levels(k$eth))
k<- cbind(k, getdt_logical(k$eth))
k[,eth:=NULL]

colnom<-c("nafld_censored", toinclude)

newrowname <- data.table(ancient = colnom, 
           new = paste0(colnom, k[,.SD,.SDcols = colnom][, sapply(.SD, function(x) ifelse(is.numeric(x), " (SD)", ", n (%)"))]))

newrowname[, new := new %>% gsub("WC", "Waist circumference, cm", . ) %>% gsub("sex", "Male", .) %>% 
             gsub("tg ", "triglycerides, mmol/L ", .) %>% gsub("hdl", "HDL, mmol/L", .) %>% gsub("ldl", "LDL, mmol/L", .) %>%
             gsub("townsend", "Townsend index", . ) %>% gsub("med", "Medication", .) %>%
             gsub("smoking", "Smoking", .) %>% gsub("alcohol", "Alcohol consumption Likert scale", .) %>%
             gsub("age_enrollment", "Age at enrollment, years", .)]

bon1 <- k[,lapply(.SD, tbl_summary),by = c("nafld_censored"), .SDcols=toinclude]
bon1[,sample:="All"]
bon1 <- merge(k[,.N,by = c("nafld_censored")], bon1, by = "nafld_censored")

bon2 <- k[toinclude==1,lapply(.SD, tbl_summary),by = c("nafld_censored"), .SDcols=toinclude]
bon2[,sample:="Subsamble"]
bon2 <- merge(k[toinclude==1,.N,by = c("nafld_censored")], bon2, by = "nafld_censored")

bon<-rbind(bon1,bon2)
long <-melt(bon, id.vars = c("nafld_censored", "sample"), measure.vars = c("N", toinclude))
final <- dcast(long,variable~sample+nafld_censored, value.var = "value")
final<-rbind(as.list(colnames(final)),as.list(colnames(final)), final)
col <- colnames(final)
final[1,(col) := lapply(.SD, function(x) gsub("_.*", "", x)), .SDcols = col ]
final[2,(col) := lapply(.SD, function(x) gsub(".*_", "", x)), .SDcols = col ]
colnames(final)<-paste0("V",1:ncol(final))
final[,roworder:=1:.N]
final <- merge( final, newrowname, by.x = "V1", by.y = "ancient", all.x=TRUE)
final <- final[order(roworder)]
final[, new:=ifelse(is.na(new), V1, new)]
final[,V1:=new]
final[,c("new","roworder"):=NULL]
fwrite(final, "Data/Modified/sample_summary_table.txt")

# run the model
tidy_cox_output <- function(k) {
  k<-summary(k)
  # toselect<- rownames(k$coefficients)[1]
  k1 <- k$coefficients[,c("exp(coef)", "Pr(>|z|)")] %>% as.data.frame(.)
  k1$exposure <- rownames(k1)
  k2<- k$conf.int[, c("lower .95", "upper .95")] %>% as.data.frame(.) 
  k2$exposure <- rownames(k2)              
  k3 <- merge(k1,k2,by="exposure") %>% as.data.table(.)               
  setnames(k3, c("exp(coef)", "Pr(>|z|)", "exposure", "lower .95", "upper .95"), 
           c("HR", "pval", "exposure", "lci", "uci"))
  k3 <- k3[,c("exposure", "HR", "lci", "uci", "pval")]
  # k3<-k3[exposure == toselect]
  return(k3)
}

tidy_ph_assumption <- function(models) {
test.ph <- cox.zph(models)
dt <- test.ph$table %>% as.data.frame()
colnames(dt)<-paste0("ph_", colnames(dt))
dt$exposure <- rownames(dt)
return(dt)
}

run_coxph_wrapper<- function(dat, #The data.frame
                             IV, #The name of the variable you are interested in as a predictor
                             DV, #The variable name need to have *_time and *_censored in the data
                             cov_inc #the character of the covariate such as " + age_enrollment + sex"
) {
  message(paste0("running on ", DV, "~", IV, cov_inc))
  formulas <- as.formula(paste0("Surv(",DV,"_time,", DV, "_censored) ~ ",IV, cov_inc))
  models <-coxph(formulas,data=dat)
  res <- tidy_cox_output(models) 
  ph <- tidy_ph_assumption(models)
  res <- merge(res, ph, by = "exposure", all = TRUE)
  res<-res[grepl(paste0(IV, "|GLOBAL"), res$exposure), ]
  res$IV<-IV
  res$outcome <- DV
  res$cov_inc<-cov_inc
  return(res)}

options(future.globals.maxSize= 5e9)
plan(multisession, workers = 10, gc = TRUE) #plan(sequential)
IV_vec <- c(regname)
k<-c("WC", "tg", "hdl")
DV_vec <- "nafld"
cov_inc_vec <- paste0(" + age_enrollment + eth + sex + med + smoking + alcohol + townsend", c("", paste0(" + ", k),
                                                                          paste0(" + ",paste(k, collapse = " + "))))
arguments <- purrr::cross(list(IV_vec, DV_vec, cov_inc_vec))
run_coxph_wrapper_safely <- safely(run_coxph_wrapper)
res <- future_map(arguments, function(x) {run_coxph_wrapper_safely(dat = dat[toinclude==1,], IV = x[[1]], DV = x[[2]], cov_inc = x[[3]])},
.options = furrr_options(seed = TRUE)) 
saveRDS(res,"Data/Modified/coxHR.rds")
index<-sapply(res, function(x) is.null(x$error))
res <- lapply(res[index], function(x) x$result) %>% rbindlist(., fill = TRUE)
fwrite(res, "Data/Modified/coxHR.txt")

#Dattg
dat[!(!is.na(nafld_date) & f.53.0.0 > nafld_date), tg_quantile := dplyr::ntile(x = tg, n = 5) %>% as.factor()]
k<-getdt_logical(dat$tg_quantile)
setnames(k, colnames(k), paste0("tg_quantile",colnames(k)))
dat <- cbind(dat, k)

arguments <-purrr::cross(list(IV_vec = c("tg","tg_quantile5", "tg_quantile"), DV_vec =  DV_vec, cov_inc_vec = cov_inc_vec[2]))
rescoxtg<- map(arguments, function(x) run_coxph_wrapper(dat = dat,
                                                        IV = x$IV_vec, DV = x$DV_vec, cov_inc = x$cov_inc_vec)) %>% 
  rbindlist(., fill = TRUE)
fwrite(rescoxtg, "Data/Modified/rescoxtg.txt")

#Table 2
rescoxtg[is.na(HR) & exposure == "tg_quantile5", exposure := "tg_quantile5TRUE"]
rescoxtg <- rescoxtg[exposure %in% paste0("tg_quantile", c("",2:5)), ]
n<-5
dat[!(!is.na(nafld_date) & f.53.0.0 > nafld_date), tg_quintile:=ntile(x = tg,n = n)%>%as.factor(.)]
quant<-stats::quantile(x=dat$tg, probs = seq(0, 1, 1/n), na.rm = TRUE)
dttranslate<-data.table(var_quantile = 1:n, var_quantile_long = as.character(NA))
for(i in 1:(length(quant)-1)) {
  dttranslate[i, ]$var_quantile_long <- paste0("[",round(quant[i], digits = 2),", ", round(quant[i+1], digits = 2), "]")
}

rescoxtg<- cbind(rescoxtg, dttranslate)

k <- dat[!(!is.na(nafld_date) & f.53.0.0 > nafld_date) & !is.na(tg),
         paste0(sum(nafld_censored), "/",.N, " (", round(sum(nafld_censored)/(.N)*100, digits =1 ), ")"), by = "tg_quintile"]
setnames(k, "V1", "event_rate")
k[,tg_quintile:=as.integer(tg_quintile)]
rescoxtg<- merge(rescoxtg, k, by.x = "var_quantile", by.y ="tg_quintile")
rescoxtg[, rescox:=paste0(format(round(HR, digits = 2), nsmall = 2) , " 95% CI=", format(round(lci, digits = 2), nsmall = 2), "-",  format(round(uci, digits = 2), nsmall = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))]
rescoxtg[var_quantile==1,rescox:="1.00"]
table2 <- rescoxtg[,.(var_quantile, var_quantile_long, event_rate, rescox)]
table2<- melt(table2, id.vars = "var_quantile", measure.vars = c("var_quantile_long", "event_rate", "rescox"))
table2<- dcast(table2,  variable ~ var_quantile, value.var = "value")
table2[,variable:=variable %>% gsub("var_quantile_long", "triglyceride range, mmol/L", .) %>% 
         gsub("event_rate", "Cases/total (event rate, %)", .) %>%
         gsub("rescox", "model", .) ]
fwrite(table2, "Data/Modified/table2.txt")

message("This script finished without errors")