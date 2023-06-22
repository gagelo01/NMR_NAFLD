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
trad <- fread( "Data/Modified/trad")

dat <- fread("Data/Modified/UKB_NMR_Data_QCed_biomarkers.txt", nrows = nrows)
regname<- setdiff(colnames(dat), c("eid", "visit_index", "NMR"))
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

my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
my.lowerthan <- function(k1, k2) {
  totransformNA<- is.na(k1)&is.na(k2)
  k1[is.na(k1)]<- Inf
  k2[is.na(k2)]<- Inf
  res <- k1 < k2 
  res[totransformNA]<-NA
  return(res)
}
dtexpect <- data.table(k1 = c(1, NA, NA, 1, 2), k2 = c(NA, 1, NA, 2, 1), expect = c(TRUE, FALSE, NA, TRUE, FALSE))
map(split(dtexpect, 1:dtexpect[,.N]), function(x) my.lowerthan(k1 = x$k1, k2 = x$k2))

exclusioncode <- colnames(wide)[!grepl(paste(c("f.eid", fieldidok),collapse = "|"), colnames(wide))]
k1 <- wide[, apply(.SD, 1, function(x) my.min(x)),.SDcols = exclusioncode]
k2 <- wide[, apply(.SD, 1, function(x) my.min(x)),.SDcols = fieldidok]
wide[, toexclude := my.lowerthan(k1, k2) ]
fwrite(wide, "Data/Modified/liverdisease_icd10_code.txt")
fidtoexclude <- wide[toexclude==TRUE, f.eid] #Those are the ID with exclusion code in the hospital inpatient record
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
dat[, age_enrollment := f.53.0.0 - date_birth]
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
dat[, nafld_time := as.IDate(censoring_date) - as.IDate(f.53.0.0)]
dat[is.na(NMR),NMR:=0]
dat[NMR==0,toinclude:=0]
col <- c("med", "smoking", "sex"); dat[, (col):=lapply(.SD, as.logical), .SDcols = col]
fwrite(dat, "Data/Modified/observationalfulldata.txt")
# summary table

k<-dat[, c("nafld_censored","WC", "sex", "tg", "hdl", "ldl", "eth",
           "townsend", "med", "smoking", "alcohol", "age_enrollment", "toinclude")]
k[,age_enrollment := age_enrollment/365.25]

getdt_logical <- function(variable_factor) {
  data<-data.table(variable_factor = variable_factor)
  data$dummy <- TRUE
  data$id <- 1:nrow(data)
  data_wide<- dcast(data, id ~ variable_factor, value.var = "dummy")
  dat <- data_wide[, lapply(.SD, function(x) ifelse(is.na(x), FALSE, TRUE)),
                   .SDcols = data[, sapply(.SD, levels),.SDcols = "variable_factor"]]
  if("NA" %in% colnames(dat)) {
  setnames(data_wide, "NA", "non_available")
  dat <- dat[, lapply(.SD, function(x) ifelse(data_wide[, !is.na(non_available)],NA ,x))]
  }
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
plan(multisession, workers = 20, gc = TRUE) #plan(sequential)
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

########Tabl2e 2 and 3 sex specific
###create a cox object
dat <- fread( "Data/Modified/observationalfulldata.txt")

dat[,sex := ifelse(sex, "male", "female")]
dat<- dat[!(!is.na(nafld_date) & f.53.0.0 > nafld_date), ]

var<-c(Triglycerides = "tg", `HDL cholesterol`  = "hdl")
varsex<-list(male = "male", female = "female", all = c("male", "female"))
for(j in 1:length(varsex)) {
  for(i in 1:length(var)) {
    dat[!(!is.na(nafld_date) & f.53.0.0 > nafld_date) & sex %in% varsex[[j]], dummy := dplyr::ntile(x = get(var[i]), n = 5) %>% as.factor(.)]
    dat[!(!is.na(nafld_date) & f.53.0.0 > nafld_date) & sex %in% varsex[[j]], (paste0(var[i], "_", names(varsex)[j])) := get(var[i])]
    k<-getdt_logical(dat$dummy)
    newname<- paste0(var[i], "_", names(varsex)[j], "_quantile") #%>% gsub("_all", "", .)
    setnames(k, colnames(k), paste0(newname, colnames(k)))
    dat <- cbind(dat, k)
    setnames(dat, "dummy", newname)
  }
}


vs <- c(sapply(var, function(x) paste0(x, "_", names(varsex)))) #%>% gsub("_all", "", .)
IV_vec <- c(sapply(vs, function(x) paste0(x, c("_quantile5", "_quantile"))), var)
arguments <-data.table(IV_vec = IV_vec,
                              DV_vec =  DV_vec,
                              cov_inc_vec = ifelse(grepl("male_|female_", IV_vec), gsub("sex + ",  "", cov_inc_vec[2] , fixed = TRUE), cov_inc_vec[2])) #remove sex as a covariate since we do stratified analyses

rescoxtg<- map(split(arguments, 1:arguments[,.N]), function(x) run_coxph_wrapper(dat = dat,
                                                        IV = x$IV_vec, DV = x$DV_vec, cov_inc = x$cov_inc_vec)) %>%
  rbindlist(., fill = TRUE)
fwrite(rescoxtg, "Data/Modified/rescoxtg.txt")


#####
###Table 2-3
rescoxtg <- fread("Data/Modified/rescoxtg.txt")
rescoxtg <- rescoxtg[!(is.na(HR) & exposure %in% paste0(vs, "_quantile5")), ]
rescoxtg <- rescoxtg[exposure %in% c(sapply(vs, function(x) paste0(x, "_quantile",c("",2:5)))), ]
n<-5

rescox <- rescoxtg
list_table <- vector(mode = "list", length = length(vs))
names(list_table)<-vs
for(i in 1:length(vs)) {
  quant<-stats::quantile(x=dat[,get(vs[i])], probs = seq(0, 1, 1/n), na.rm = TRUE)
  dttranslate<-data.table(var_quantile = 1:n, var_quantile_long = as.character(NA))
  for(j in 1:(length(quant)-1)) {
    dttranslate[j, ]$var_quantile_long <- paste0("[",round(quant[j], digits = 2),", ", round(quant[j+1], digits = 2), "]")
  }
  rescoxtg<- cbind(rescox, dttranslate)
  rescoxtg <- rescoxtg[exposure %in% c(sapply(vs[i], function(x) paste0(x, "_quantile",c("",2:5)))), ]

  k <- dat[!(!is.na(nafld_date) & f.53.0.0 > nafld_date) & !is.na(get(vs[i])),
           paste0(sum(nafld_censored), "/",.N, " (", round(sum(nafld_censored)/(.N)*100, digits =1 ), ")"), by = eval(paste0(vs[i], "_quantile"))]
  setnames(k, "V1", "event_rate")
  setnames(k, paste0(vs[i], "_quantile"), "var_quantile")
  k[,var_quantile:=as.integer(var_quantile)]
  rescoxtg<- merge(rescoxtg, k, by = "var_quantile")
  rescoxtg[, rescox:=paste0(format(round(HR, digits = 2), nsmall = 2) , " 95% CI=", format(round(lci, digits = 2), nsmall = 2), "-",  format(round(uci, digits = 2), nsmall = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))]
  rescoxtg[var_quantile==1,rescox:="1.00"]
  table2 <- rescoxtg[,.(var_quantile, var_quantile_long, event_rate, rescox)]
  table2<- melt(table2, id.vars = "var_quantile", measure.vars = c("var_quantile_long", "event_rate", "rescox"))
  table2<- dcast(table2,  variable ~ var_quantile, value.var = "value")
  table2 <- cbind(Sex = c(gsub(".*_" , "", vs[i]), rep("", nrow(table2)-1)), table2)
  kk <- c(male = "Men", female = "Women", all = "All")
  table2[,Sex := kk[Sex]]
  lookup_table<-c(var_quantile_long = paste0(names(var[var==vs[i] %>% gsub("_.*","",.)]), " range, mmol/L"), event_rate = "Cases/total (event rate, %)", rescox = "model")
  table2[,variable:= lookup_table[variable]]
  list_table[[i]]<-table2
}

saveRDS(list_table, "Data/Modified/list_tablesexspecific.rds")

####
#create the var combined_dyslipidemia
dat <- fread( "Data/Modified/observationalfulldata.txt")
dt<-dat[!(!is.na(nafld_date) & f.53.0.0 > nafld_date),]
dt[, high_hdl := (sex==TRUE & hdl > 1)|(sex==FALSE & hdl > 1.3)]
dt[, low_tg := tg<1.7]
dt[, combined := 4-((high_hdl+low_tg*2))]
dt[, combined_dyslipidemia := paste0(high_hdl, "_", low_tg)]
# dt[,  combined_dyslipidemia := paste("group", combined)]
dt[is.na(combined),combined_dyslipidemia:=NA]
level <- c(sapply(c("TRUE", "FALSE"), function(x) paste0(x, "_", c("TRUE", "FALSE"))))
dt[,combined_dyslipidemia := factor(combined_dyslipidemia, levels = level)]
# dt[,combined_dyslipidemia := factor(combined_dyslipidemia, levels = paste("group",1:4))]
dt[,sex:=ifelse(sex, "Men", "Women")%>%as.factor(.)]

dt_plot <- dt[!(is.na(high_hdl)|is.na(low_tg)),]
dt_plot[,label1 := ifelse(high_hdl,"high HDL-C", "low HDL-C")]
dt_plot[,label2 := ifelse(low_tg,"low TG", "high TG")]
dt_plot[,label_var := paste0(label2, " + ", label1)]
dt_plot[,label_var := factor(label_var, levels = c("low TG + high HDL-C", "low TG + low HDL-C", "high TG + high HDL-C", "high TG + low HDL-C"))]
dt_plot[,sex:=factor(sex, levels = c("Men", "Women"))]

# fit <- survfit(Surv(nafld_time/365.25, nafld_censored) ~ label_var, data = dt_plot)
# k <- ggsurvplot_facet(
#   fit, data = dt_plot, facet.by = "sex", ylim = c(0.95,1), censor.size = 0.2,size = 0.5,
#   xlab = "Follow-up (Years)", ylab = "Diagnosis-free survival",
#   font.legend = 12,  palette = "jco", risk.table =  "absolute",
#   tables.theme = theme_cleantable())
# 
# k + theme(legend.position = "top") +
#   guides(color=guide_legend(nrow=2,byrow=TRUE))

# funcplot <- function(dt_plot,
#                      legend_code = TRUE #either FALSE or guide_legend(nrow = 2)
#                      ) {
# fit <- survfit(Surv(nafld_time/365.25, nafld_censored) ~ label_var, data = dt_plot)
# k <- ggsurvplot(
#   fit, data = dt_plot, ylim = c(0.95,1), censor.size = 0.2,size = 0.5,
#   xlab = "Follow-up (Years)", ylab = "Diagnosis-free survival",
#   font.legend = 12,  palette = "jco", risk.table =  "absolute",
#   tables.theme = theme_cleantable(), legend.labs = levels(dt_plot$label_var))
# k + guides(color = legend_code)
# }
# # 
# p1<-funcplot(dt_plot[sex=="Men",], legend_code = guide_legend(nrow = 2))
# p2<-funcplot(dt_plot[sex=="Women",], legend_code = guide_legend(nrow = 2))
# 
# 
# 
# png(filename = "Data/Modified/p1.png")
# funcplot(dt_plot[sex=="Men",], legend_code = guide_legend(nrow = 2))
# dev.off()
# png(filename = "Data/Modified/p2.png")
# p2
# dev.off()
# plot1 <- readPNG("Data/Modified/p1.png")
# plot2 <- readPNG('plot2.png')

# grid.arrange(rasterGrob(plot1),rasterGrob(plot2),ncol=1, )

# cowplot::plot_grid(p1, p2, labels = c('Men', 'Women'))
# ggpubr::ggarrange(p1, p2, 
                  # labels = c("A)", "B)"),
                  # ncol = 2, nrow = 1)


# 
source("Analysis/functionto_facetrisktable.R")
setDF(dt_plot)
fit <- survfit(Surv(nafld_time/365.25, nafld_censored) ~ label_var + sex, data = dt_plot)

res<-300
tiff(file="Results/Figure4.tiff", width = 623/72*res, height = 414/72*res, res = res)
ggsurvplot_facet_risktable(fit, dt_plot,  ylim = c(0.95,1), censor.size = 0.2,size = 0.5,
                             xlab = "Follow-up (Years)", ylab = "Diagnosis-free survival",
                             font.legend = 12,  palette = "jco")

dev.off()
rm(list="dt_plot")
# ggsave(filename = "Results/Figure4.png", plot = k, device = "png", width = 500/72, height = 360/72, units = "in")

varsex<-list(list(Men = "Men"), list(Women = "Women"), list(All = c("Men", "Women")))
rescoxtg<- map(varsex, function(x) {
  rescox <- run_coxph_wrapper(dat = dt[sex %in% x[[1]], ],
                              IV = "combined_dyslipidemia", DV = DV_vec,
                              cov_inc = ifelse(names(x)%in%"All", cov_inc_vec[2],gsub("sex + ",  "", cov_inc_vec[2], fixed = TRUE)))
  rescox[, sex := names(x)]
  return(rescox)}) %>%
  rbindlist(., fill = TRUE)


rescoxtg[, rescox:=paste0(format(round(HR, digits = 2), nsmall = 2) , " 95% CI=", format(round(lci, digits = 2), nsmall = 2), "-",  format(round(uci, digits = 2), nsmall = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))]
rescoxtg <- rescoxtg[!is.na(HR), ]
rescoxtg[,exposure := gsub("combined_dyslipidemia", "", exposure)]
rescoxtg <- rescoxtg[,.(exposure, sex, rescox)]
rescoxtg <- rbind(rescoxtg, data.table(exposure = "TRUE_TRUE", sex = rescoxtg$sex %>% unique, rescox = 1))
rescoxtg <- separate(rescoxtg,col = "exposure", into = c("high_hdl", "low_tg"))
rescoxtg[, (c("high_hdl", "low_tg")) := lapply(.SD, as.logical), .SDcols = c("high_hdl", "low_tg")]


caseevent <- dt[, paste0(sum(nafld_censored), "/",.N, " (", round(sum(nafld_censored)/(.N)*100, digits =1 ), ")"), by= c("high_hdl", "low_tg", "sex")]
caseeventall <- dt[, paste0(sum(nafld_censored), "/",.N, " (", round(sum(nafld_censored)/(.N)*100, digits =1 ), ")"), by= c("high_hdl", "low_tg")]
caseeventall[,sex:= "All"]
caseevent <- rbind(caseevent, caseeventall)
caseevent <- caseevent[!(is.na(high_hdl)|is.na(low_tg)),]
caseevent <- merge(caseevent, rescoxtg, by = c("high_hdl", "low_tg", "sex"))
caseevent[,(c("high_hdl", "low_tg")) := lapply(.SD, function(x) factor(x, levels = c("TRUE", "FALSE"))), .SDcols = c("high_hdl", "low_tg")]
caseevent <- melt(caseevent, id.vars = c("high_hdl", "low_tg", "sex"), measure.vars = c("V1", "rescox"))
caseevent <- dcast(caseevent,  sex + variable ~ high_hdl + low_tg, value.var = c("value"))
vectrad <- c(V1 = "Cases/total (event rate, %)", rescox = "Hazard ratios for incident NAFLD")
caseevent[,variable := vectrad[variable]]

fwrite(caseevent, "Data/Modified/tablecombineddyslipidemia.txt")


message("This script finished without errors")
