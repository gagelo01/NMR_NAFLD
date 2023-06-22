#!/usr/bin/env Rscript
#import instrument
library(data.table)
library(tidyverse)
library(ukbnmr)
setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Triglycerides_dis")

nmr_df <- fread("/home/couchr02/Mendel_UKB/Source/Phenotype/December_2021_Metabolomics_complete/ukb49738.tab")
nmr_data <- ukbnmr::extract_biomarkers(x = nmr_df)
nmr_processed <- ukbnmr::remove_technical_variation(nmr_df[f.eid %in% nmr_data$eid,])
dat <- nmr_processed$biomarkers
setDT(dat)
regname <- colnames(dat)[!(colnames(dat) %in% c("eid", "visit_index")) ]
dat[,M_LDL_FC_by_CE:=gsub("Inf", "NA", M_LDL_FC_by_CE)%>%as.numeric(.)]
dat[,NMR:=1]
dat<- dat[, .SD[1], by = "eid"] #I choose the first measurement for every participants. 1426 participants with two measurments. 3712 with only measurment at time 1. 116595 with only measurment at time 0
dat[,c(regname) := lapply(.SD, function(x) (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)), .SDcols = regname] #standardise

fwrite(dat, "Data/Modified/UKB_NMR_Data_QCed_biomarkers.txt")

message("This script finished without errors")
