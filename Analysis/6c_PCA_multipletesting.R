#!/usr/bin/env Rscript
#import instrument
library(data.table)
library(tidyverse)
library(furrr)
library("pcaMethods")

setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Triglycerides_dis")
nrows<-Inf#Inf
dat <- fread("Data/Modified/UKB_NMR_Data_QCed_biomarkers.txt", nrows = nrows)
trad <- fread( "Data/Modified/trad")
pca_dat <- dat[, .SD, .SDcols = trad$id]
pcares <- pcaMethods::pca(object = as.matrix(pca_dat), 
                     method = "nipals",
                     scale = "uv",
                     center = TRUE,
                     nPcs =20 )
ntest <- min(which(pcares@R2cum > 0.9))

saveRDS(pcares, "Data/Modified/pcares.rds")

message("This script finished without errors")