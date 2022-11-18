#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
library(tictoc)
library(furrr)

setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Triglycerides_dis")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao_small <- ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/"), ]
options(future.globals.maxSize= 5e9)
plan(multicore, workers = 20, gc = TRUE)#plan(sequential)

############What you need to change ###########
ID_mrbase_exp <- c("ukb-b-9405", "ieu-a-61", ao_small[pmid ==32203549,]$id, ao_small[author == "Borges CM", id])
exp_mrbase <- paste0("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/", ID_mrbase_exp, "/", ID_mrbase_exp, ".vcf.gz")
ID_server_exp <-   c("trait-1-1", "trait-14-8", paste0("trait-16-", c(1,2,4)), "dis-2-1", "dis-6-1", "dis-8-1")
exp_server <- paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ID_server_exp, "/", ID_server_exp, ".vcf.gz")
arguments <- data.table(id = c(ID_mrbase_exp,ID_server_exp), path = c(exp_mrbase,exp_server), pval = 5e-8, r2 =0.01, kb = 1000)
#################################################
inst <- future_map(split(arguments, 1:nrow(arguments)), function(x) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  GagnonMR::get_inst(x$path, r2 = x$r2, pval = x$pval, kb = x$kb)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)
setDT(inst)
out <- future_map(as.list(c(exp_mrbase,exp_server)), function(x, rsiid = unique(inst$SNP)) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  outr<-GagnonMR::extract_outcome_variant(snps = rsiid, outcomes = x)
  return(outr)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)

instmvmr <- TwoSampleMR::convert_outcome_to_exposure(out) %>% as.data.table(.)
dtn <- data.table(id = c("ieu-b-110", "ieu-b-111", "ieu-b-108", "ieu-b-109", "ieu-b-107"),
                  samplesize = c(440546, 441016, 439214, 403943, 393193))
list_res <- map(list(inst=inst,instmvmr=instmvmr), function(x) {
  x[, id.exposure := exposure]
  x[exposure %in% ao_small[author == "Borges CM", id], samplesize.exposure := 115078]
  x <- merge(x, dtn, by.x = "exposure", by.y = "id", all.x = TRUE)
  x[!is.na(samplesize), samplesize.exposure := samplesize]
  x[,samplesize:=NULL]
})

inst_all_sign_clump <- list_res$inst
split_inst <- split(inst_all_sign_clump, inst_all_sign_clump$exposure)
ressd <- map(split_inst, function(x) {
  sd <- coloc:::sdY.est(vbeta = (x$se.exposure)^2,
                        maf = x[, ifelse(eaf.exposure > 0.5,1-eaf.exposure,eaf.exposure)],
                        n = x$samplesize.exposure)
  return(data.table(exposure = x$exposure[1], sd = sd))
}) %>% rbindlist(., fill = TRUE) %>% as.data.table(.)

ressd <- ressd[!(exposure %in% df_index[id %in% ID_server_exp,trait]), ]

u<-c("logOR", "log odds", "log orr")
k <- c(df_index[unit %in% u,trait],  ao_small[unit %in% u, id]) #remove binary as exposure
ressd <- ressd[!(exposure %in% k),]
# ressd$sd %>% hist #all traits have been inverse ranked normal transform

####extract TM?SF2 and PNPLA3
inst<-GagnonMR::get_inst("/mnt/sda/gagelo01/Vcffile/Server_vcf/dis-2-1/dis-2-1.vcf.gz")
main_SNP <- inst[order(pval.exposure), ][1:2,SNP]
test<-gwasvcf::get_ld_proxies(main_SNP, bfile = ldref, tag_kb = 5000, tag_r2 = 0.8) %>% as.data.table(.)
test<-test[,.(SNP_A,SNP_B)]
test <- rbind(test, data.table(SNP_A = unique(test$SNP_A), SNP_B = unique(test$SNP_A)))
test <- merge(test, data.table(SNP_A = c("rs3747207", "rs73001065"), gene.exposure = c("PNPLA3", "TM6SF2")), by = "SNP_A")
inst_all_sign_clump <- fread( "Data/Modified/inst_all_sign_clump.txt")
inst_all_sign_clump<-merge(inst_all_sign_clump, test[,.(SNP_B, gene.exposure)], by.x = "SNP", by.y = "SNP_B", all.x = TRUE)
#####

fwrite((list_res$instmvmr), "Data/Modified/all_inst_mvmr.txt")
fwrite(GagnonMR::convert_exposure_to_outcome(list_res$instmvmr), "Data/Modified/all_outcome_mvmr.txt" )
fwrite(list_res$inst, "Data/Modified/inst_all_sign_clump.txt")

message("This script finished without errors")
