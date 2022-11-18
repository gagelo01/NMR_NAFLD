#!/usr/bin/env Rscript
#import instrument
library(data.table)
library(tidyverse)
library(ukbnmr)
setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Triglycerides_dis")

nmr_df <- fread("/home/couchr02/Mendel_UKB/Source/Phenotype/December_2021_Metabolomics_complete/ukb49738.tab")
nmr_data <- ukbnmr::extract_biomarkers(x = nmr_df)
nmr_processed <- ukbnmr::remove_technical_variation(nmr_df[f.eid %in% nmr_data$eid,])

fwrite(nmr_processed$biomarkers, "Data/Modified/UKB_NMR_Data_QCed_biomarkers.txt")

message("This script finished without errors")
