setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
SNP_data=read.csv("SNP_table10.csv")

trait=c()

colnames(SNP_data)

snplist = SNP_data[,"RS.number"]
#all_SNP=read.table("ALLSNPrsID.txt")

#for (i in snplist){
#  if (!(i %in% all_SNP$V1)){
#    print(i)
#  }
#}



#data
all_data=read.table("all_data_imputation_new.txt",header = T)
PC_data=read.table("VTE_50388_PC1-20.txt",header=T)
PC_data=PC_data[,1:11]
colnames(PC_data)=c("eid","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")

geno <- geno[-1,]

for (i in 3:dim(geno)[2]) {
  geno[,i] <- as.numeric(as.character(geno[,i]))
}
all_data=merge(all_data,PC_data[,1:11],by.x="eid",by.y = "eid")
all_data$eGFR_scr_10=all_data$eGFR_scr/10
all_data$eGFR_scys_10=all_data$eGFR_scys/10
all_data$sex_cat=0
all_data$sex_cat[which(all_data$sex=="female")]=1
all_data$sex=all_data$sex_cat
#sub_all_data=all_data[,c("eid","eGFR_scr_10","eGFR_scys_10","BUN","age","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")]
#write.csv(sub_all_data,"all_data_score.txt",quote = F,row.names = F)

score_scr=read.table("eGFR_scr_10.profile_score",header=T)
score_scys=read.table("eGFR_cyst_10.profile_score",header=T)
score_BUN=read.table("BUN.profile_score",header=T)

all_data_scr=merge(all_data,score_scr,by.x="eid",by.y = "FID")
all_data_scr=all_data_scr[which(!is.na(all_data_scr$eGFR_scr_10)),]
all_data_scys=merge(all_data,score_scys,by.x="eid",by.y = "FID")
all_data_scys=all_data_scys[which(!is.na(all_data_scys$eGFR_scys_10)),]
all_data_BUN=merge(all_data,score_BUN,by.x="eid",by.y = "FID")
all_data_BUN=all_data_BUN[which(!is.na(all_data_BUN$BUN)),]


#linear MR

#R2 and F
#eGFR_scr
forml = "eGFR_scr_10 ~ SCORE"
R1=summary(lm(as.formula(forml),data=all_data_scr))$r.squared#0.009956432
F1=summary(lm(as.formula(forml),data=all_data_scr))$fstatistic[1]#86.17686 
RF1=data.frame(Exposure="eGFR_scr",R2=R1, F=F1)
#eGFR_scys
forml = "eGFR_scys_10 ~ SCORE"
R2=summary(lm(as.formula(forml),data=all_data_scys))$r.squared#0.008672422
F2=summary(lm(as.formula(forml),data=all_data_scys))$fstatistic[1]#74.99582  
RF2=data.frame(Exposure="eGFR_scys",R2=R2, F=F2)
#BUN
forml = "BUN ~ SCORE"
R3=summary(lm(as.formula(forml),data=all_data_BUN))$r.squared#0.004859592
F3=summary(lm(as.formula(forml),data=all_data_BUN))$fstatistic[1]#41.83787
RF3=data.frame(Exposure="BUN",R2=R3, F=F3)

RF=rbind(RF1,RF2,RF3)
write.csv(RF,"RF_252_new.csv")

#beta1
#eGFR_scr
fit=lm("eGFR_scr_10 ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10",data=all_data_scr)
beta1_scr=summary(fit)$coef[2,1]
#eGFR_scys
fit=lm("eGFR_scys_10 ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10",data=all_data_scys)
beta1_scys=summary(fit)$coef[2,1]
#eGFR_scr
fit=lm("BUN ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10",data=all_data_BUN)
beta1_BUN=summary(fit)$coef[2,1]

#beta2
#eGFR_scr
#All infections
all_data_scr$follow_up_time = all_data_scr$all_infection_time
all_data_scr$event.type = all_data_scr$all_infection_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scr)
beta2_scr_all_infections = summary(fit2)$coefficients[1,1]
se2_scr_all_infections = summary(fit2)$coefficients[1,3]
MR_est_scr_all_infections = beta2_scr_all_infections/beta1_scr 
MR_se_scr_all_infections = se2_scr_all_infections/beta1_scr
samples11=summary(fit2)$n
events11=summary(fit2)$nevent

#pneumonia
all_data_scr$follow_up_time = all_data_scr$pneumonia_time
all_data_scr$event.type = all_data_scr$pneumonia_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scr)
beta2_scr_pneumonia = summary(fit2)$coefficients[1,1]
se2_scr_pneumonia = summary(fit2)$coefficients[1,3]
MR_est_scr_pneumonia = beta2_scr_pneumonia/beta1_scr 
MR_se_scr_pneumonia = se2_scr_pneumonia/beta1_scr
samples12=summary(fit2)$n
events12=summary(fit2)$nevent

#sepsis
all_data_scr$follow_up_time = all_data_scr$sepsis_time
all_data_scr$event.type = all_data_scr$sepsis_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scr)
beta2_scr_sepsis = summary(fit2)$coefficients[1,1]
se2_scr_sepsis = summary(fit2)$coefficients[1,3]
MR_est_scr_sepsis = beta2_scr_sepsis/beta1_scr 
MR_se_scr_sepsis = se2_scr_sepsis/beta1_scr
samples13=summary(fit2)$n
events13=summary(fit2)$nevent


#skin_soft_tissue
all_data_scr$follow_up_time = all_data_scr$skin_soft_tissue_time
all_data_scr$event.type = all_data_scr$skin_soft_tissue_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scr)
beta2_scr_skin_soft_tissue = summary(fit2)$coefficients[1,1]
se2_scr_skin_soft_tissue = summary(fit2)$coefficients[1,3]
MR_est_scr_skin_soft_tissue = beta2_scr_skin_soft_tissue/beta1_scr 
MR_se_scr_skin_soft_tissue = se2_scr_skin_soft_tissue/beta1_scr
samples14=summary(fit2)$n
events14=summary(fit2)$nevent


#urinary_tract
all_data_scr$follow_up_time = all_data_scr$urinary_tract_time
all_data_scr$event.type = all_data_scr$urinary_tract_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scr)
beta2_scr_urinary_tract = summary(fit2)$coefficients[1,1]
se2_scr_urinary_tract = summary(fit2)$coefficients[1,3]
MR_est_scr_urinary_tract = beta2_scr_urinary_tract/beta1_scr 
MR_se_scr_urinary_tract = se2_scr_urinary_tract/beta1_scr
samples15=summary(fit2)$n
events15=summary(fit2)$nevent


#eGFR_scys
#All infections
all_data_scys$follow_up_time = all_data_scys$all_infection_time
all_data_scys$event.type = all_data_scys$all_infection_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scys)
beta2_scys_all_infections = summary(fit2)$coefficients[1,1]
se2_scys_all_infections = summary(fit2)$coefficients[1,3]
MR_est_scys_all_infections = beta2_scys_all_infections/beta1_scys 
MR_se_scys_all_infections = se2_scys_all_infections/beta1_scys
samples21=summary(fit2)$n
events21=summary(fit2)$nevent

#pneumonia
all_data_scys$follow_up_time = all_data_scys$pneumonia_time
all_data_scys$event.type = all_data_scys$pneumonia_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scys)
beta2_scys_pneumonia = summary(fit2)$coefficients[1,1]
se2_scys_pneumonia = summary(fit2)$coefficients[1,3]
MR_est_scys_pneumonia = beta2_scys_pneumonia/beta1_scys 
MR_se_scys_pneumonia = se2_scys_pneumonia/beta1_scys
samples22=summary(fit2)$n
events22=summary(fit2)$nevent

#sepsis
all_data_scys$follow_up_time = all_data_scys$sepsis_time
all_data_scys$event.type = all_data_scys$sepsis_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scys)
beta2_scys_sepsis = summary(fit2)$coefficients[1,1]
se2_scys_sepsis = summary(fit2)$coefficients[1,3]
MR_est_scys_sepsis = beta2_scys_sepsis/beta1_scys 
MR_se_scys_sepsis = se2_scys_sepsis/beta1_scys
samples23=summary(fit2)$n
events23=summary(fit2)$nevent


#skin_soft_tissue
all_data_scys$follow_up_time = all_data_scys$skin_soft_tissue_time
all_data_scys$event.type = all_data_scys$skin_soft_tissue_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scys)
beta2_scys_skin_soft_tissue = summary(fit2)$coefficients[1,1]
se2_scys_skin_soft_tissue = summary(fit2)$coefficients[1,3]
MR_est_scys_skin_soft_tissue = beta2_scys_skin_soft_tissue/beta1_scys 
MR_se_scys_skin_soft_tissue = se2_scys_skin_soft_tissue/beta1_scys
samples24=summary(fit2)$n
events24=summary(fit2)$nevent


#urinary_tract
all_data_scys$follow_up_time = all_data_scys$urinary_tract_time
all_data_scys$event.type = all_data_scys$urinary_tract_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scys)
beta2_scys_urinary_tract = summary(fit2)$coefficients[1,1]
se2_scys_urinary_tract = summary(fit2)$coefficients[1,3]
MR_est_scys_urinary_tract = beta2_scys_urinary_tract/beta1_scys 
MR_se_scys_urinary_tract = se2_scys_urinary_tract/beta1_scys
samples25=summary(fit2)$n
events25=summary(fit2)$nevent


#BUN
#All infections
all_data_BUN$follow_up_time = all_data_BUN$all_infection_time
all_data_BUN$event.type = all_data_BUN$all_infection_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_BUN)
beta2_BUN_all_infections = summary(fit2)$coefficients[1,1]
se2_BUN_all_infections = summary(fit2)$coefficients[1,3]
MR_est_BUN_all_infections = beta2_BUN_all_infections/beta1_BUN 
MR_se_BUN_all_infections = se2_BUN_all_infections/beta1_BUN
samples31=summary(fit2)$n
events31=summary(fit2)$nevent

#pneumonia
all_data_BUN$follow_up_time = all_data_BUN$pneumonia_time
all_data_BUN$event.type = all_data_BUN$pneumonia_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_BUN)
beta2_BUN_pneumonia = summary(fit2)$coefficients[1,1]
se2_BUN_pneumonia = summary(fit2)$coefficients[1,3]
MR_est_BUN_pneumonia = beta2_BUN_pneumonia/beta1_BUN 
MR_se_BUN_pneumonia = se2_BUN_pneumonia/beta1_BUN
samples32=summary(fit2)$n
events32=summary(fit2)$nevent

#sepsis
all_data_BUN$follow_up_time = all_data_BUN$sepsis_time
all_data_BUN$event.type = all_data_BUN$sepsis_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_BUN)
beta2_BUN_sepsis = summary(fit2)$coefficients[1,1]
se2_BUN_sepsis = summary(fit2)$coefficients[1,3]
MR_est_BUN_sepsis = beta2_BUN_sepsis/beta1_BUN 
MR_se_BUN_sepsis = se2_BUN_sepsis/beta1_BUN
samples33=summary(fit2)$n
events33=summary(fit2)$nevent


#skin_soft_tissue
all_data_BUN$follow_up_time = all_data_BUN$skin_soft_tissue_time
all_data_BUN$event.type = all_data_BUN$skin_soft_tissue_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_BUN)
beta2_BUN_skin_soft_tissue = summary(fit2)$coefficients[1,1]
se2_BUN_skin_soft_tissue = summary(fit2)$coefficients[1,3]
MR_est_BUN_skin_soft_tissue = beta2_BUN_skin_soft_tissue/beta1_BUN 
MR_se_BUN_skin_soft_tissue = se2_BUN_skin_soft_tissue/beta1_BUN
samples34=summary(fit2)$n
events34=summary(fit2)$nevent


#urinary_tract
all_data_BUN$follow_up_time = all_data_BUN$urinary_tract_time
all_data_BUN$event.type = all_data_BUN$urinary_tract_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_BUN)
beta2_BUN_urinary_tract = summary(fit2)$coefficients[1,1]
se2_BUN_urinary_tract = summary(fit2)$coefficients[1,3]
MR_est_BUN_urinary_tract = beta2_BUN_urinary_tract/beta1_BUN 
MR_se_BUN_urinary_tract = se2_BUN_urinary_tract/beta1_BUN
samples35=summary(fit2)$n
events35=summary(fit2)$nevent


linear_MR_result = data.frame(exposure=c("eGFR_cr","eGFR_cr","eGFR_cr","eGFR_cr","eGFR_cr",
                                         "eGFR_cys","eGFR_cys","eGFR_cys","eGFR_cys","eGFR_cys",
                                         "BUN","BUN","BUN","BUN","BUN"),
                              outcome=c("All infections","Pneumonia","Sepsis","Skin and soft tissue infections","Urinary tract infections",
                                        "All infections","Pneumonia","Sepsis","Skin and soft tissue infections","Urinary tract infections",
                                        "All infections","Pneumonia","Sepsis","Skin and soft tissue infections","Urinary tract infections"),
                              linear_MR_estimate = c(-MR_est_scr_all_infections,-MR_est_scr_pneumonia,-MR_est_scr_sepsis,-MR_est_scr_skin_soft_tissue,-MR_est_scr_urinary_tract,
                                                     -MR_est_scys_all_infections,-MR_est_scys_pneumonia,-MR_est_scys_sepsis,-MR_est_scys_skin_soft_tissue,-MR_est_scys_urinary_tract,
                                                     MR_est_BUN_all_infections,MR_est_BUN_pneumonia,MR_est_BUN_sepsis,MR_est_BUN_skin_soft_tissue,MR_est_BUN_urinary_tract),
                              linear_MR_SE = c(MR_se_scr_all_infections,MR_se_scr_pneumonia,MR_se_scr_sepsis,MR_se_scr_skin_soft_tissue,MR_se_scr_urinary_tract,
                                               MR_se_scys_all_infections,MR_se_scys_pneumonia,MR_se_scys_sepsis,MR_se_scys_skin_soft_tissue,MR_se_scys_urinary_tract,
                                               MR_se_BUN_all_infections,MR_se_BUN_pneumonia,MR_se_BUN_sepsis,MR_se_BUN_skin_soft_tissue,MR_se_BUN_urinary_tract))

linear_MR_result_df = as.data.frame(linear_MR_result)
linear_MR_result_df$est_upper = linear_MR_result_df$linear_MR_estimate + 1.96*linear_MR_result_df$linear_MR_SE
linear_MR_result_df$est_lower = linear_MR_result_df$linear_MR_estimate - 1.96*linear_MR_result_df$linear_MR_SE
linear_MR_result_df$HR = exp(linear_MR_result_df$linear_MR_estimate)
linear_MR_result_df$HR_upper = exp(linear_MR_result_df$est_upper)
linear_MR_result_df$HR_lower = exp(linear_MR_result_df$est_lower)
linear_MR_result_df$p_value = ifelse(linear_MR_result_df$linear_MR_estimate<=0,
                                     2*pnorm(linear_MR_result_df$linear_MR_estimate/linear_MR_result_df$linear_MR_SE),
                                     2*(1-pnorm(linear_MR_result_df$linear_MR_estimate/linear_MR_result_df$linear_MR_SE)))

linear_MR_result_df$sample=c(samples11,samples12,samples13,samples14,samples15,
                             samples21,samples22,samples23,samples24,samples25,
                             samples31,samples32,samples33,samples34,samples35)
linear_MR_result_df$event=c(events11,events12,events13,events14,events15,
                            events21,events22,events23,events24,events25,
                            events31,events32,events33,events34,events35)

linear_MR_result_df

write.csv(linear_MR_result_df,"linear_MR_252SNP_new.csv")


##################################################################################################
#nonlinear MR
library(metafor); library(ggplot2)
library(ggpubr)
library(metafor); library(ggplot2)
library(ggpubr)
#install.packages('ggsci')
library("RColorBrewer")
mycolors<-brewer.pal(9, "Blues")

source("nlme_summ_aes.R")  


#sTEP1&2
#eGFR_scr
all_data_scr$scr_hat=lm(eGFR_scr_10~SCORE,data=all_data_scr)$fit
all_data_scr$scr_residual=all_data_scr$eGFR_scr_10-all_data_scr$scr_hat
N=10
all_data_scr$group=cut(all_data_scr$scr_residual, labels = F,quantile(all_data_scr$scr_residual, prob = 0:N / N, names = FALSE), include = TRUE)
#eGFR_scys
all_data_scys$scys_hat=lm(eGFR_scys_10~SCORE,data=all_data_scys)$fit
all_data_scys$scys_residual=all_data_scys$eGFR_scys_10-all_data_scys$scys_hat
N=10
all_data_scys$group=cut(all_data_scys$scys_residual, labels = F,quantile(all_data_scys$scys_residual, prob = 0:N / N, names = FALSE), include = TRUE)
#BUN
all_data_BUN$BUN_hat=lm(BUN~SCORE,data=all_data_BUN)$fit
all_data_BUN$BUN_residual=all_data_BUN$BUN-all_data_BUN$BUN_hat
N=10
all_data_BUN$group=cut(all_data_BUN$BUN_residual, labels = F,quantile(all_data_BUN$BUN_residual, prob = 0:N / N, names = FALSE), include = TRUE)


#STEP345


#non-linear MR: All_infections ~ eGFR_scr

all_data_scr$follow_up_time = all_data_scr$all_infection_time
all_data_scr$event.type=all_data_scr$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = all_data_scr[which(all_data_scr$group==i),]
  temp = summary(lm(eGFR_scr_10~SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scr_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmean0= summary_stat$mean_expos # mean exposure in each quantile
xmean0
xmean1 = xmean0-6  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

a1 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=2, offset=6, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR (creatinine)", pref_y="Hazard ratio of All Infections", breaks=c(0.25,0.5,1,2,4,8))
a1

a1$coefficients
a1$p_tests
a1$p_heterogeneity

fp_a1=ifelse(a1$p_tests[1,2]<0.001," < 0.001",ifelse(a1$p_tests[1,2]>0.999," > 0.999",paste(" =",round(a1$p_tests[1,2],3))))
Q_a1=ifelse(a1$p_tests[1,4]<0.001," < 0.001",ifelse(a1$p_tests[1,4]>0.999," > 0.999",paste(" =",round(a1$p_tests[1,4],3))))

f11=ggplot() + 
  geom_histogram(aes(x = all_data_scr$eGFR_scr, y = ((..count..) / sum(..count..)*10) ), data.frame(all_data_scr$eGFR_scr), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),a1$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), a1$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),a1$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmean0)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cr]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cr"]," and All Infections")))+
  labs(x =  expression(paste(eGFR["cr"]," (","ml/min/1.73 ",m^2,")")), y="Hazard ratio")+ 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.8), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")))+
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_a1) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_a1)),  size = 3)

f11


#non-linear MR: pneumonia ~ eGFR_scr

all_data_scr$follow_up_time = all_data_scr$pneumonia_time
all_data_scr$event.type=all_data_scr$pneumonia_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = all_data_scr[which(all_data_scr$group==i),]
  temp = summary(lm(eGFR_scr_10~SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scr_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmean0 = summary_stat$mean_expos # mean exposure in each quantile
xmean0
xmean1 = xmean0-6  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

b1 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=2, offset=6, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR (creatinine)", pref_y="Hazard ratio of Pneumonia", breaks=c(0.25,0.5,1,2,4,8))
b1

b1$coefficients
b1$p_tests
b1$p_heterogeneity

fp_b1=ifelse(b1$p_tests[1,2]<0.001," < 0.001",ifelse(b1$p_tests[1,2]>0.999," > 0.999",paste(" =",round(b1$p_tests[1,2],3))))
Q_b1=ifelse(b1$p_tests[1,4]<0.001," < 0.001",ifelse(b1$p_tests[1,4]>0.999," > 0.999",paste(" =",round(b1$p_tests[1,4],3))))

f12=ggplot() + 
  geom_histogram(aes(x = all_data_scr$eGFR_scr, y = ((..count..) / sum(..count..)*10) ), data.frame(all_data_scr$eGFR_scr), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),b1$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), b1$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),b1$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmean0)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cr]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cr"]," and Pneumonia")))+
  labs(x =  expression(paste(eGFR["cr"]," (","ml/min/1.73 ",m^2,")")), y="Hazard ratio")+ 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.8), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")))+
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_b1) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_b1)),  size = 3)

f12

#non-linear MR: sepsis ~ eGFR_scr

all_data_scr$follow_up_time = all_data_scr$sepsis_time
all_data_scr$event.type=all_data_scr$sepsis_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = all_data_scr[which(all_data_scr$group==i),]
  temp = summary(lm(eGFR_scr_10~SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scr_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmean0 = summary_stat$mean_expos # mean exposure in each quantile
xmean0
xmean1 = xmean0-6  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

c1 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=2, offset=6, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR (creatinine)", pref_y="Hazard ratio of Sepsis", breaks=c(0.25,0.5,1,2,4,8))
c1

c1$coefficients
c1$p_tests
c1$p_heterogeneity
fp_c1=ifelse(c1$p_tests[1,2]<0.001," < 0.001",ifelse(c1$p_tests[1,2]>0.999," > 0.999",paste(" =",round(c1$p_tests[1,2],3))))
Q_c1=ifelse(c1$p_tests[1,4]<0.001," < 0.001",ifelse(c1$p_tests[1,4]>0.999," > 0.999",paste(" =",round(c1$p_tests[1,4],3))))

f13=ggplot() + 
  geom_histogram(aes(x = all_data_scr$eGFR_scr, y = ((..count..) / sum(..count..)*10) ), data.frame(all_data_scr$eGFR_scr), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),c1$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), c1$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),c1$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmean0)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cr]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cr"]," and Sepsis")))+
  labs(x =  expression(paste(eGFR["cr"]," (","ml/min/1.73 ",m^2,")")), y="Hazard ratio")+ 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.8), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")))+
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_c1) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_c1)),  size = 3)

f13

#non-linear MR: skin_soft_tissue ~ eGFR_scr

all_data_scr$follow_up_time = all_data_scr$skin_soft_tissue_time
all_data_scr$event.type=all_data_scr$skin_soft_tissue_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = all_data_scr[which(all_data_scr$group==i),]
  temp = summary(lm(eGFR_scr_10~SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scr_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmean0 = summary_stat$mean_expos # mean exposure in each quantile
xmean0
xmean1 = xmean0-6  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

d1 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=2, offset=6, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR (creatinine)", pref_y="Hazard ratio of Skin and Soft Tissue Infections", breaks=c(0.25,0.5,1,2,4,8))
d1

d1$coefficients
d1$p_tests
d1$p_heterogeneity
fp_d1=ifelse(d1$p_tests[1,2]<0.001," < 0.001",ifelse(d1$p_tests[1,2]>0.999," > 0.999",paste(" =",round(d1$p_tests[1,2],3))))
Q_d1=ifelse(d1$p_tests[1,4]<0.001," < 0.001",ifelse(d1$p_tests[1,4]>0.999," > 0.999",paste(" =",round(d1$p_tests[1,4],3))))

f14=ggplot() + 
  geom_histogram(aes(x = all_data_scr$eGFR_scr, y = ((..count..) / sum(..count..)*10) ), data.frame(all_data_scr$eGFR_scr), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),d1$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), d1$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),d1$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmean0)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cr]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cr"]," and Skin and Soft Tissue Infections")))+
  labs(x =  expression(paste(eGFR["cr"]," (","ml/min/1.73 ",m^2,")")), y="Hazard ratio")+ 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.8), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")))+
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_d1) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_d1)),  size = 3)

f14


#non-linear MR: urinary_tract ~ eGFR_scr

all_data_scr$follow_up_time = all_data_scr$urinary_tract_time
all_data_scr$event.type=all_data_scr$urinary_tract_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = all_data_scr[which(all_data_scr$group==i),]
  temp = summary(lm(eGFR_scr_10~SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scr_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmean0 = summary_stat$mean_expos # mean exposure in each quantile
xmean0
xmean1 = xmean0-6  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

e1 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=2, offset=6, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR (creatinine)", pref_y="Hazard ratio of Urinary Tract Infections", breaks=c(0.25,0.5,1,2,4,8))
e1

e1$coefficients
e1$p_tests
e1$p_heterogeneity

fp_e1=ifelse(e1$p_tests[1,2]<0.001," < 0.001",ifelse(e1$p_tests[1,2]>0.999," > 0.999",paste(" =",round(e1$p_tests[1,2],3))))
Q_e1=ifelse(e1$p_tests[1,4]<0.001," < 0.001",ifelse(e1$p_tests[1,4]>0.999," > 0.999",paste(" =",round(e1$p_tests[1,4],3))))

f15=ggplot() + 
  geom_histogram(aes(x = all_data_scr$eGFR_scr, y = ((..count..) / sum(..count..)*10) ), data.frame(all_data_scr$eGFR_scr), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),e1$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), e1$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),e1$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmean0)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cr]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cr"]," and Urinary Tract Infections")))+
  labs(x =  expression(paste(eGFR["cr"]," (","ml/min/1.73 ",m^2,")")), y="Hazard ratio")+ 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.8), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")))+
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_e1) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_e1)),  size = 3)

f15


#non-linear MR: All_infections ~ eGFR_scys

all_data_scys$follow_up_time = all_data_scys$all_infection_time
all_data_scys$event.type=all_data_scys$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = all_data_scys[which(all_data_scys$group==i),]
  temp = summary(lm(eGFR_scys_10~SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scys_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmean2 = summary_stat$mean_expos # mean exposure in each quantile
xmean2
xmean1 = xmean2-5  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

a2 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=2, offset=5, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_scys", pref_y="Hazard ratio of all infections", breaks=c(0.25,0.5,1,2,4,8))
a2

a2$coefficients
a2$p_tests
a2$p_heterogeneity
fp_a2=ifelse(a2$p_tests[1,2]<0.001," < 0.001",ifelse(a2$p_tests[1,2]>0.999," > 0.999",paste(" =",round(a2$p_tests[1,2],3))))
Q_a2=ifelse(a2$p_tests[1,4]<0.001," < 0.001",ifelse(a2$p_tests[1,4]>0.999," > 0.999",paste(" =",round(a2$p_tests[1,4],3))))

f21=ggplot() + 
  geom_histogram(aes(x = all_data_scys$eGFR_scys, y = ((..count..) / sum(..count..)*10) ), data.frame(all_data_scys$eGFR_scys), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),a2$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), a2$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),a2$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmean2)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cys]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cys"]," and All Infections")))+
  labs(x =  expression(paste(eGFR["cys"]," (","ml/min/1.73 ",m^2,")")), y="Hazard ratio")+ 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.8), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")))+
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_a2) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_a2)),  size = 3)

f21


#non-linear MR: pneumonia ~ eGFR_scys

all_data_scys$follow_up_time = all_data_scys$pneumonia_time
all_data_scys$event.type=all_data_scys$pneumonia_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = all_data_scys[which(all_data_scys$group==i),]
  temp = summary(lm(eGFR_scys_10~SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scys_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmean2 = summary_stat$mean_expos # mean exposure in each quantile
xmean2
xmean1 = xmean2-5  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

b2 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=2, offset=5, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_scys", pref_y="Hazard ratio of pneumonia", breaks=c(0.25,0.5,1,2,4,8))
b2

b2$coefficients
b2$p_tests
b2$p_heterogeneity
fp_b2=ifelse(b2$p_tests[1,2]<0.001," < 0.001",ifelse(b2$p_tests[1,2]>0.999," > 0.999",paste(" =",round(b2$p_tests[1,2],3))))
Q_b2=ifelse(b2$p_tests[1,4]<0.001," < 0.001",ifelse(b2$p_tests[1,4]>0.999," > 0.999",paste(" =",round(b2$p_tests[1,4],3))))

f22=ggplot() + 
  geom_histogram(aes(x = all_data_scys$eGFR_scys, y = ((..count..) / sum(..count..)*10) ), data.frame(all_data_scys$eGFR_scys), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),b2$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), b2$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),b2$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmean2)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cys]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cys"]," and Pneumonia")))+
  labs(x =  expression(paste(eGFR["cys"]," (","ml/min/1.73 ",m^2,")")), y="Hazard ratio")+ 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.8), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")))+
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_b2) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_b2)),  size = 3)

f22


#non-linear MR: sepsis ~ eGFR_scys

all_data_scys$follow_up_time = all_data_scys$sepsis_time
all_data_scys$event.type=all_data_scys$sepsis_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = all_data_scys[which(all_data_scys$group==i),]
  temp = summary(lm(eGFR_scys_10~SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scys_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmean2 = summary_stat$mean_expos # mean exposure in each quantile
xmean2
xmean1 = xmean2-5  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

c2 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=2, offset=5, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_scys", pref_y="Hazard ratio of sepsis", breaks=c(0.25,0.5,1,2,4,8))
c2

c2$coefficients
c2$p_tests
c2$p_heterogeneity

fp_c2=ifelse(c2$p_tests[1,2]<0.001," < 0.001",ifelse(c2$p_tests[1,2]>0.999," > 0.999",paste(" =",round(c2$p_tests[1,2],3))))
Q_c2=ifelse(c2$p_tests[1,4]<0.001," < 0.001",ifelse(c2$p_tests[1,4]>0.999," > 0.999",paste(" =",round(c2$p_tests[1,4],3))))

f23=ggplot() + 
  geom_histogram(aes(x = all_data_scys$eGFR_scys, y = ((..count..) / sum(..count..)*10) ), data.frame(all_data_scys$eGFR_scys), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),c2$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), c2$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),c2$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmean2)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cys]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cys"]," and Sepsis")))+
  labs(x =  expression(paste(eGFR["cys"]," (","ml/min/1.73 ",m^2,")")), y="Hazard ratio")+ 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.8), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")))+
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_c2) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_c2)),  size = 3)

f23


#non-linear MR: skin_soft_tissue ~ eGFR_scys

all_data_scys$follow_up_time = all_data_scys$skin_soft_tissue_time
all_data_scys$event.type=all_data_scys$skin_soft_tissue_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = all_data_scys[which(all_data_scys$group==i),]
  temp = summary(lm(eGFR_scys_10~SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scys_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmean2 = summary_stat$mean_expos # mean exposure in each quantile
xmean2
xmean1 = xmean2-5  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

d2 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=2, offset=5, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_scys", pref_y="Hazard ratio of Skin and Soft Tissue Infections", breaks=c(0.25,0.5,1,2,4,8))
d2

d2$coefficients
d2$p_tests
d2$p_heterogeneity

fp_d2=ifelse(d2$p_tests[1,2]<0.001," < 0.001",ifelse(d2$p_tests[1,2]>0.999," > 0.999",paste(" =",round(d2$p_tests[1,2],3))))
Q_d2=ifelse(d2$p_tests[1,4]<0.001," < 0.001",ifelse(d2$p_tests[1,4]>0.999," > 0.999",paste(" =",round(d2$p_tests[1,4],3))))

f24=ggplot() + 
  geom_histogram(aes(x = all_data_scys$eGFR_scys, y = ((..count..) / sum(..count..)*10) ), data.frame(all_data_scys$eGFR_scys), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),d2$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), d2$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),d2$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmean2)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cys]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cys"]," and Skin and Soft Tissue Infections")))+
  labs(x =  expression(paste(eGFR["cys"]," (","ml/min/1.73 ",m^2,")")), y="Hazard ratio")+ 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.8), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")))+
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_d2) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_d2)),  size = 3)

f24


#non-linear MR: urinary_tract ~ eGFR_scys

all_data_scys$follow_up_time = all_data_scys$urinary_tract_time
all_data_scys$event.type=all_data_scys$urinary_tract_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = all_data_scys[which(all_data_scys$group==i),]
  temp = summary(lm(eGFR_scys_10~SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scys_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmean2 = summary_stat$mean_expos # mean exposure in each quantile
xmean2
xmean1 = xmean2-5  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

e2 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=2, offset=5, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_scys", pref_y="Hazard ratio of urinary_tract", breaks=c(0.25,0.5,1,2,4,8))
e2

e2$coefficients
e2$p_tests
e2$p_heterogeneity

fp_e2=ifelse(e2$p_tests[1,2]<0.001," < 0.001",ifelse(e2$p_tests[1,2]>0.999," > 0.999",paste(" =",round(e2$p_tests[1,2],3))))
Q_e2=ifelse(e2$p_tests[1,4]<0.001," < 0.001",ifelse(e2$p_tests[1,4]>0.999," > 0.999",paste(" =",round(e2$p_tests[1,4],3))))

f25=ggplot() + 
  geom_histogram(aes(x = all_data_scys$eGFR_scys, y = ((..count..) / sum(..count..)*10) ), data.frame(all_data_scys$eGFR_scys), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),e2$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), e2$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),e2$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmean2)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cys]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cys"]," and Urinary Tract Infections")))+
  labs(x =  expression(paste(eGFR["cys"]," (","ml/min/1.73 ",m^2,")")), y="Hazard ratio")+ 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.8), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")))+
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_e2) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_e2)),  size = 3)

f25


#non-linear MR: All_infections ~ BUN

all_data_BUN$follow_up_time = all_data_BUN$all_infection_time
all_data_BUN$event.type=all_data_BUN$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = all_data_BUN[which(all_data_BUN$group==i),]
  temp = summary(lm(BUN~SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$BUN)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmean3 = summary_stat$mean_expos # mean exposure in each quantile
xmean3
xmean1 = xmean3-3  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

a3 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=2, offset=3, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="BUN", pref_y="Hazard ratio of all infection", breaks=c(0.25,0.5,1,2,4,8))
a3

a3$coefficients
a3$p_tests
a3$p_heterogeneity

fp_a3=ifelse(a3$p_tests[1,2]<0.001," < 0.001",ifelse(a3$p_tests[1,2]>0.999," > 0.999",paste(" =",round(a3$p_tests[1,2],3))))
Q_a3=ifelse(a3$p_tests[1,4]<0.001," < 0.001",ifelse(a3$p_tests[1,4]>0.999," > 0.999",paste(" =",round(a3$p_tests[1,4],3))))

f31=ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), all_data_BUN,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,12,1))+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x, y=log2(yest)),a3$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x, y=log2(lci)), a3$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x, y=log2(uci)),a3$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmean3), y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 4),  labels =c(0.25,0.5,1,2,4,8,16),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~BUN))) + 
  scale_x_continuous(breaks=seq(0,12,by = 1)) + 
  coord_cartesian(xlim =c(0,12)) + ggtitle("BUN and All Infections")+
  labs(x =  expression(paste(BUN," (","mmol/L)")), y="Hazard ratio")+ 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.8), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")))+
  annotate("text", x = 10, y = 3.7, label = bquote(P['Cochran\'s Q test'] ~ .(Q_a3) ), size = 3) + 
  annotate("text", x = 10, y = 3.4, label = bquote(P['nonlinearity'] ~.(fp_a3)),  size = 3)

f31


#non-linear MR: pneumonia ~ BUN

all_data_BUN$follow_up_time = all_data_BUN$pneumonia_time
all_data_BUN$event.type=all_data_BUN$pneumonia_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = all_data_BUN[which(all_data_BUN$group==i),]
  temp = summary(lm(BUN~SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$BUN)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmean3 = summary_stat$mean_expos # mean exposure in each quantile
xmean3
xmean1 = xmean3-3  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

b3 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=2, offset=3, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="BUN", pref_y="Hazard ratio of pneumonia", breaks=c(0.25,0.5,1,2,4,8))
b3

b3$coefficients
b3$p_tests
b3$p_heterogeneity

fp_b3=ifelse(b3$p_tests[1,2]<0.001," < 0.001",ifelse(b3$p_tests[1,2]>0.999," > 0.999",paste(" =",round(b3$p_tests[1,2],3))))
Q_b3=ifelse(b3$p_tests[1,4]<0.001," < 0.001",ifelse(b3$p_tests[1,4]>0.999," > 0.999",paste(" =",round(b3$p_tests[1,4],3))))

f32=ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), all_data_BUN,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,12,1))+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x, y=log2(yest)),b3$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x, y=log2(lci)), b3$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x, y=log2(uci)),b3$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmean3), y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 4),  labels =c(0.25,0.5,1,2,4,8,16),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~BUN))) + 
  scale_x_continuous(breaks=seq(0,12,by = 1)) + 
  coord_cartesian(xlim =c(0,12)) + ggtitle("BUN and Pneumonia")+
  labs(x =  expression(paste(BUN," (","mmol/L)")), y="Hazard ratio")+ 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.8), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")))+
  annotate("text", x = 10, y = 3.7, label = bquote(P['Cochran\'s Q test'] ~ .(Q_b3) ), size = 3) + 
  annotate("text", x = 10, y = 3.4, label = bquote(P['nonlinearity'] ~.(fp_b3)),  size = 3)

f32



#non-linear MR: sepsis ~ BUN

all_data_BUN$follow_up_time = all_data_BUN$sepsis_time
all_data_BUN$event.type=all_data_BUN$sepsis_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = all_data_BUN[which(all_data_BUN$group==i),]
  temp = summary(lm(BUN~SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$BUN)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmean3 = summary_stat$mean_expos # mean exposure in each quantile
xmean3
xmean1 = xmean3-3  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

c3 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=2, offset=3, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="BUN", pref_y="Hazard ratio of sepsis", breaks=c(0.25,0.5,1,2,4,8))
c3

c3$coefficients
c3$p_tests
c3$p_heterogeneity

fp_c3=ifelse(c3$p_tests[1,2]<0.001," < 0.001",ifelse(c3$p_tests[1,2]>0.999," > 0.999",paste(" =",round(c3$p_tests[1,2],3))))
Q_c3=ifelse(c3$p_tests[1,4]<0.001," < 0.001",ifelse(c3$p_tests[1,4]>0.999," > 0.999",paste(" =",round(c3$p_tests[1,4],3))))

f33=ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), all_data_BUN,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,12,1))+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x, y=log2(yest)),c3$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x, y=log2(lci)), c3$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x, y=log2(uci)),c3$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmean3), y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 4),  labels =c(0.25,0.5,1,2,4,8,16),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~BUN))) + 
  scale_x_continuous(breaks=seq(0,12,by = 1)) + 
  coord_cartesian(xlim =c(0,12)) + ggtitle("BUN and Sepsis")+
  labs(x =  expression(paste(BUN," (","mmol/L)")), y="Hazard ratio")+ 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.8), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")))+
  annotate("text", x = 10, y = 3.7, label = bquote(P['Cochran\'s Q test'] ~ .(Q_c3) ), size = 3) + 
  annotate("text", x = 10, y = 3.4, label = bquote(P['nonlinearity'] ~.(fp_c3)),  size = 3)

f33


#non-linear MR: skin_soft_tissue ~ BUN

all_data_BUN$follow_up_time = all_data_BUN$skin_soft_tissue_time
all_data_BUN$event.type=all_data_BUN$skin_soft_tissue_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = all_data_BUN[which(all_data_BUN$group==i),]
  temp = summary(lm(BUN~SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$BUN)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmean3 = summary_stat$mean_expos # mean exposure in each quantile
xmean3
xmean1 = xmean3-3  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

d3 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=2, offset=3, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="BUN", pref_y="Hazard ratio of skin_soft_tissue", breaks=c(0.25,0.5,1,2,4,8))
d3

d3$coefficients
d3$p_tests
d3$p_heterogeneity

fp_d3=ifelse(d3$p_tests[1,2]<0.001," < 0.001",ifelse(d3$p_tests[1,2]>0.999," > 0.999",paste(" =",round(d3$p_tests[1,2],3))))
Q_d3=ifelse(d3$p_tests[1,4]<0.001," < 0.001",ifelse(d3$p_tests[1,4]>0.999," > 0.999",paste(" =",round(d3$p_tests[1,4],3))))

f34=ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), all_data_BUN,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,12,1))+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x, y=log2(yest)),d3$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x, y=log2(lci)), d3$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x, y=log2(uci)),d3$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmean3), y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 4),  labels =c(0.25,0.5,1,2,4,8,16),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~BUN))) + 
  scale_x_continuous(breaks=seq(0,12,by = 1)) + 
  coord_cartesian(xlim =c(0,12)) + ggtitle("BUN and Skin and Soft Tissue Infections")+
  labs(x =  expression(paste(BUN," (","mmol/L)")), y="Hazard ratio")+ 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.8), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")))+
  annotate("text", x = 10, y = 3.7, label = bquote(P['Cochran\'s Q test'] ~ .(Q_d3) ), size = 3) + 
  annotate("text", x = 10, y = 3.4, label = bquote(P['nonlinearity'] ~.(fp_d3)),  size = 3)

f34


#non-linear MR: urinary_tract ~ BUN

all_data_BUN$follow_up_time = all_data_BUN$urinary_tract_time
all_data_BUN$event.type=all_data_BUN$urinary_tract_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = all_data_BUN[which(all_data_BUN$group==i),]
  temp = summary(lm(BUN~SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$BUN)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ SCORE+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmean3 = summary_stat$mean_expos # mean exposure in each quantile
xmean3
xmean1 = xmean3-3  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

e3 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=2, offset=3, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="BUN", pref_y="Hazard ratio of Urinary Tract Infections", breaks=c(0.25,0.5,1,2,4,8))
e3

e3$coefficients
e3$p_tests
e3$p_heterogeneity
fp_e3=ifelse(e3$p_tests[1,2]<0.001," < 0.001",ifelse(e3$p_tests[1,2]>0.999," > 0.999",paste(" =",round(e3$p_tests[1,2],3))))
Q_e3=ifelse(e3$p_tests[1,4]<0.001," < 0.001",ifelse(e3$p_tests[1,4]>0.999," > 0.999",paste(" =",round(e3$p_tests[1,4],3))))

f35=ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), all_data_BUN,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,12,1))+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x, y=log2(yest)),e3$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x, y=log2(lci)), e3$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x, y=log2(uci)),e3$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmean3), y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 4),  labels =c(0.25,0.5,1,2,4,8,16),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~BUN))) + 
  scale_x_continuous(breaks=seq(0,12,by = 1)) + 
  coord_cartesian(xlim =c(0,12)) + ggtitle("BUN and Urinary Tract Infections")+
  labs(x =  expression(paste(BUN," (","mmol/L)")), y="Hazard ratio")+ 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.8), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")))+
  annotate("text", x = 10, y = 3.7, label = bquote(P['Cochran\'s Q test'] ~ .(Q_e3) ), size = 3) + 
  annotate("text", x = 10, y = 3.4, label = bquote(P['nonlinearity'] ~.(fp_e3)),  size = 3)

f35






library(cowplot)

pdf("Nonlinear_MR_252SNP_new.pdf",family = "Times",height = 20,width = 15)

plot_grid(f11,f21,f31,
          f12,f22,f32,
          f13,f23,f33,
          f14,f24,f34,
          f15,f25,f35,
          nrow = 5,ncol =3)
dev.off()








