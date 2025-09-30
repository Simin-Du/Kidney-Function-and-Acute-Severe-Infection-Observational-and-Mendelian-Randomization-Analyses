#install.packages('devtools')
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library('devtools')
library("remotes")
library("TwoSampleMR")
library("MRInstruments")
library("MRPRESSO")

# devtools::install_local('./MRInstruments.zip')
# devtools::install_local('./MRPRESSO.zip')
# devtools::install_local('./TwoSampleMR.zip')
 
#工具变量
eGFR_file=read.csv("252SNP_infor.csv",header=T)
#snp32_list=read.table("31SNP_1Proxy.txt")
#eGFR_file=eGFR_file[which(eGFR_file$SNP %in% snp32_list[,]),]

eGFR_exp_dat<-format_data(
  eGFR_file,
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col ="effect_allele",
  other_allele_col = "othere_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  samplesize_col = "samplesize",
  gene_col="Closest.Gene"
)
eGFR_exp_dat$exposure="eGFR"
exp_dat_clump <-clump_data(eGFR_exp_dat,clump_r2=0.01,clump_kb=5000) #remove 35 snp


#结局变量

ao <- available_outcomes(access_token=NULL)
write.csv(ao,file="available_traits.csv",quote = F,row.names = F)

library(TwoSampleMR)
##调取数据库中BMI GWAS中显著SNP，并以“r2=0.01”/5000loci/window做clumping

dim(exp_dat_clump)
#[1] 217 14
# clumping后剩余80个SNP
#Sepsis
sepsis_out <- extract_outcome_data(
  snps=exp_dat_clump$SNP,
  outcomes='ieu-b-69',
  proxies = FALSE,
  maf_threshold = 0.01,
  access_token = NULL
)
dim(sepsis_out)
#[1] 198 16
#提取80个IV SNPs

#Harmonizing effect:
mydata <- harmonise_data(
  exposure_dat=exp_dat_clump,
  outcome_dat=sepsis_out,
  action= 2
)

#Run MR model with 4 methods:
res1 <- mr(mydata)
exp(res1$b)
res1
#mr_method_list()

##heterogeneity testing异质性检验
het <- mr_heterogeneity(mydata)
het
#IV之间存在很强的异质性(Q_pval<0.05)，需要去除outcome GWAS中P值非常小的SNP
#如果不去SNP，则可用随机效应模型估计MR效应：
mr(mydata,method_list=c('mr_ivw_mre')) 
#从随机效应模型的结果看出，BMI和T2D确实存在因果关系（pval < 0.05），
#BMI升高导致T2D发病风险增高（b > 0）

##pleiotropy testing多效性检验
pleio <- mr_pleiotropy_test(mydata)
pleio
#MR-Egger中egger_intercept和0没有显著差异（pval > 0.05）
#因此没有充分证据证明水平多效性的存在

##leave-one-out testing
single <- mr_leaveoneout(mydata)
mr_leaveoneout_plot(single)
#无论去除哪个SNP都不会对结果产生根本影响（所有CI均在0右侧）
#说明MR结果稳健


# # Urinary Tract Infection
# Urinary_tract_out <- extract_outcome_data(
#   snps=exp_dat_clump$SNP,
#   outcomes='ukb-b-6077',
#   proxies = FALSE,
#   maf_threshold = 0.01,
#   access_token = NULL
# )
# dim(Urinary_tract_out)
# #[1] 203 16
# #提取80个IV SNPs
# 
# #Harmonizing effect:
# mydata <- harmonise_data(
#   exposure_dat=exp_dat_clump,
#   outcome_dat=Urinary_tract_out,
#   action= 2
# )
# 
# #Run MR model with 4 methods:
# res2 <- mr(mydata)
# res2
# exp(res2$b-1.96*res2$se)
# #mr_method_list()
# 
# ##heterogeneity testing异质性检验
# het <- mr_heterogeneity(mydata)
# het
# #IV之间存在很强的异质性(Q_pval<0.05)，需要去除outcome GWAS中P值非常小的SNP
# #如果不去SNP，则可用随机效应模型估计MR效应：
# mr(mydata,method_list=c('mr_ivw_mre')) 
# #从随机效应模型的结果看出，BMI和T2D确实存在因果关系（pval < 0.05），
# #BMI升高导致T2D发病风险增高（b > 0）
# 
# ##pleiotropy testing多效性检验
# pleio <- mr_pleiotropy_test(mydata)
# pleio
# #MR-Egger中egger_intercept和0没有显著差异（pval > 0.05）
# #因此没有充分证据证明水平多效性的存在
# 
# ##leave-one-out testing
# single <- mr_leaveoneout(mydata)
# mr_leaveoneout_plot(single)
# #无论去除哪个SNP都不会对结果产生根本影响（所有CI均在0右侧）
# #说明MR结果稳健


# Skin and subsutaneous
Skin_out <- extract_outcome_data(
  snps=exp_dat_clump$SNP,
  outcomes='finn-b-L12_CELLULITIS',
  proxies = FALSE,
  maf_threshold = 0.01,
  access_token = NULL
)
dim(Skin_out)
#[1] 209 16
#提取80个IV SNPs

#Harmonizing effect:
mydata <- harmonise_data(
  exposure_dat=exp_dat_clump,
  outcome_dat=Skin_out,
  action= 2
)

#Run MR model with 4 methods:
res3 <- mr(mydata)
res3
exp(res3$b)
#mr_method_list()

##heterogeneity testing异质性检验
het <- mr_heterogeneity(mydata)
het
#IV之间存在很强的异质性(Q_pval<0.05)，需要去除outcome GWAS中P值非常小的SNP
#如果不去SNP，则可用随机效应模型估计MR效应：
mr(mydata,method_list=c('mr_ivw_mre')) 
#从随机效应模型的结果看出，BMI和T2D确实存在因果关系（pval < 0.05），
#BMI升高导致T2D发病风险增高（b > 0）

##pleiotropy testing多效性检验
pleio <- mr_pleiotropy_test(mydata)
pleio
#MR-Egger中egger_intercept和0没有显著差异（pval > 0.05）
#因此没有充分证据证明水平多效性的存在

##leave-one-out testing
single <- mr_leaveoneout(mydata)
mr_leaveoneout_plot(single)
#无论去除哪个SNP都不会对结果产生根本影响（所有CI均在0右侧）
#说明MR结果稳健

# Pneumonia
Pneumonia_out <- extract_outcome_data(
  snps=exp_dat_clump$SNP,
  outcomes='finn-b-J10_PNEUMONIA',
  proxies = FALSE,
  maf_threshold = 0.01,
  access_token = NULL
)
dim(Pneumonia_out)
#[1] 217 16
#提取80个IV SNPs

#Harmonizing effect:
mydata <- harmonise_data(
  exposure_dat=exp_dat_clump,
  outcome_dat=Pneumonia_out,
  action= 2
)

#Run MR model with 4 methods:
res4 <- mr(mydata)
res4
exp(res4$b)
#mr_method_list()

##heterogeneity testing异质性检验
het <- mr_heterogeneity(mydata)
het
#IV之间存在很强的异质性(Q_pval<0.05)，需要去除outcome GWAS中P值非常小的SNP
#如果不去SNP，则可用随机效应模型估计MR效应：
mr(mydata,method_list=c('mr_ivw_mre')) 
#从随机效应模型的结果看出，BMI和T2D确实存在因果关系（pval < 0.05），
#BMI升高导致T2D发病风险增高（b > 0）

##pleiotropy testing多效性检验
pleio <- mr_pleiotropy_test(mydata)
pleio
#MR-Egger中egger_intercept和0没有显著差异（pval > 0.05）
#因此没有充分证据证明水平多效性的存在

##leave-one-out testing
single <- mr_leaveoneout(mydata)
mr_leaveoneout_plot(single)
#无论去除哪个SNP都不会对结果产生根本影响（所有CI均在0右侧）
#说明MR结果稳健

two_sample_results=rbind(res1,res3,res4)
two_sample_results$beta_lower95=two_sample_results$b-1.96*two_sample_results$se
two_sample_results$beta_upper95=two_sample_results$b+1.96*two_sample_results$se
write.csv(two_sample_results,"twosampleMR_results.csv")

#######PLOTING
##散点-回归图
mr_scatter_plot(res,mydata)
#BMI升高，T2D的发病风险升高

#森林图forest plot
res_single <- mr_singlesnp(mydata)
mr_forest_plot(res_single)
#IVW方法下，BMI的升高增加糖尿病的发病风险

#漏斗图
mr_funnel_plot(res_single)
#左侧某点明显偏离了总体，建议把它去除掉。
#该结果和SNP间的较大异质性很大的结果相符。

devtools::install_github("rondolab/MR-PRESSO")
library("MRPRESSO")
data(SummaryStats) 
