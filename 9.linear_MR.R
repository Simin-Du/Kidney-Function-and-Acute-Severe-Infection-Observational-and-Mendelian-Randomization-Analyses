setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
geno = read.table("31SNP_1Proxy.xmat",head=T)
all_data=read.table("all_data_imputation_new.txt",header = T)
PC_data=read.table("VTE_50388_PC1-20.txt",header=T)
PC_data=PC_data[,1:11]
colnames(PC_data)=c("eid","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")
#data_c2 = all_data
#load("SNP_afterclumping.Rdata")
#rm_SNP <- readRDS("rm_cyst.rds")
###注意：rm_SNP是需去除的SNP，注意你的list里是去除的SNP，还是剩余的SNP。

geno <- geno[-1,]

for (i in 3:dim(geno)[2]) {
  geno[,i] <- as.numeric(as.character(geno[,i]))
}

all_data=merge(all_data,geno,by.x="eid",by.y="FID")
all_data=merge(all_data,PC_data[,1:11],by.x="eid",by.y = "eid")
snp_list=read.table("31SNP_1Proxy.txt")
all_data$eGFR_scr_10=all_data$eGFR_scr/10
all_data$eGFR_scys_10=all_data$eGFR_scys/10

all_data$education=factor(all_data$education,ordered=T)
all_data$smoking=factor(all_data$smoking,ordered=T)
all_data$alcohol=factor(all_data$alcohol,ordered=T)
all_data$CVD=factor(all_data$CVD,ordered=T)
all_data$cancer=factor(all_data$cancer,ordered=T)
all_data$diabetes=factor(all_data$diabetes,ordered=T)
all_data$UACR_cat=factor(all_data$UACR_cat,ordered=T)
all_data$hypertension=factor(all_data$hypertension,ordered=T)

all_data_scr=all_data[which(!is.na(all_data$eGFR_scr_10)),]
all_data_scys=all_data[which(!is.na(all_data$eGFR_scys_10)),]
all_data_BUN=all_data[which(!is.na(all_data$BUN)),]
#####
#R2 and F
#eGFR_scr
forml = paste0("eGFR_scr_10 ~ ",paste0(snp_list$V1,collapse = "+"))
R1=summary(lm(as.formula(forml),data=all_data))$r.squared#0.009956432
F1=summary(lm(as.formula(forml),data=all_data))$fstatistic[1]#86.17686 
all_data$GIV_eGFR_scr=predict(lm(as.formula(forml),data=all_data),newdata = all_data)
RF1=data.frame(Exposure="eGFR_scr",R2=R1, F=F1)
#eGFR_scys
forml = paste0("eGFR_scys_10 ~ ",paste0(snp_list$V1,collapse = "+"))
R2=summary(lm(as.formula(forml),data=all_data))$r.squared#0.008672422
F2=summary(lm(as.formula(forml),data=all_data))$fstatistic[1]#74.99582  
all_data$GIV_eGFR_scys=predict(lm(as.formula(forml),data=all_data),newdata = all_data)
RF2=data.frame(Exposure="eGFR_scys",R2=R2, F=F2)
#BUN
forml = paste0("BUN ~ ",paste0(snp_list$V1,collapse = "+"))
R3=summary(lm(as.formula(forml),data=all_data))$r.squared#0.004859592
F3=summary(lm(as.formula(forml),data=all_data))$fstatistic[1]#41.83787
all_data$GIV_BUN=predict(lm(as.formula(forml),data=all_data),newdata = all_data)
RF3=data.frame(Exposure="BUN",R2=R3, F=F3)

RF=rbind(RF1,RF2,RF3)




#GIV
#eGFR_scr
forml = paste0("eGFR_scr_10 ~ age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+",paste0(snp_list$V1,collapse = "+"))
all_data_scr$GIV_eGFR_scr=predict(lm(as.formula(forml),data=all_data),newdata = all_data_scr)
fit1=lm(as.formula(forml),data=all_data)
res=summary(fit1)
sample1=length(all_data$eid)-length(res$n)
#eGFR_scYS
forml = paste0("eGFR_scys_10 ~ age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+",paste0(snp_list$V1,collapse = "+"))
all_data_scys$GIV_eGFR_scys=predict(lm(as.formula(forml),data=all_data),newdata = all_data_scys)
fit1=lm(as.formula(forml),data=all_data)
res=summary(fit1)
sample2=length(all_data$eid)-length(res$n)
#eGFR_BUN
forml = paste0("BUN ~ age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+",paste0(snp_list$V1,collapse = "+"))
all_data_BUN$GIV_BUN=predict(lm(as.formula(forml),data=all_data),newdata = all_data_BUN)
fit1=lm(as.formula(forml),data=all_data)
res=summary(fit1)
sample3=length(all_data$eid)-length(res$n)

samples=data.frame(exposure=c("eGFR_scr","eGFR_scys","BUN"),samples=c(sample1,sample2,sample3))


library(tableone)
forml = paste0("eGFR_scr_10 ~ ",paste0(snp_list$V1,collapse = "+"))
all_data_scr$GIV=predict(lm(as.formula(forml),data=all_data),newdata = all_data_scr)

all_data_scr[which(all_data_scr$GIV<quantile(na.omit(all_data_scr$GIV),probs = 1)),"GIV_cat"]=4
all_data_scr[which(all_data_scr$GIV<quantile(na.omit(all_data_scr$GIV),probs = 0.75)),"GIV_cat"]=3
all_data_scr[which(all_data_scr$GIV<quantile(na.omit(all_data_scr$GIV),probs = 0.5)),"GIV_cat"]=2
all_data_scr[which(all_data_scr$GIV<quantile(na.omit(all_data_scr$GIV),probs = 0.25)),"GIV_cat"]=1
#all_data_scr[which(all_data_scr$GIV<quantile(na.omit(all_data_scr$GIV),probs = 0.2)),"GIV_cat"]=0

var=c("sex","age","education","TDI","smoking","alcohol","MET","CVD",
      "cancer","diabetes","triglycerides","glucose","HbA1c","UACR_cat","BMI","WC",
      "SBP","DBP","mean_grip_strength","hypertension","unic_acid","c_reactive_protein","BUN","low_density_lipoprotein_cl",
      "high_density_lipoprotein_cl","vitamin_D" )
description=CreateTableOne(vars=var,strata="GIV_cat",all_data_scr,addOverall=T)
tab1<-print(description, nonnormal = c("sex","age","education","TDI","smoking","alcohol","MET","CVD",
                                       "cancer","diabetes","triglycerides","glucose","HbA1c","UACR_cat","BMI","WC",
                                       "SBP","DBP","mean_grip_strength","hypertension","unic_acid","c_reactive_protein","BUN","low_density_lipoprotein_cl",
                                       "high_density_lipoprotein_cl","vitamin_D"),formatOptions = list(big.mark = ","),quote = T)
write.csv(tab1,"Table1_imp_GIV_scr.csv")

#outcome~ GIV+age+sex+10PCs
#eGFR_scr

#all_infection
all_data_scr$follow_up_time = all_data_scr$all_infection_time
all_data_scr$event.type = all_data_scr$all_infection_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ GIV_eGFR_scr+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scr)
MR_CE = -summary(fit2)$coefficients[1,1]
MR_SE = summary(fit2)$coefficients[1,3]
MR_HR=1/summary(fit2)$coefficients[1,2]
MR_pvalue=summary(fit2)$coefficients[1,5]
MR_HR_lower95=1/summary(fit2)$conf.int[1,4]
MR_HR_upper95=1/summary(fit2)$conf.int[1,3]
samples=summary(fit2)$n
events=summary(fit2)$nevent
#Cox_R2=summary(fit2)$rsq[[1]]
#Cox_F=(summary(fit2)$rsq[[1]]/13)/((1-summary(fit2)$rsq[[1]])/(summary(fit2)$n-13-1))

scr_allinfection = data.frame(Exposure="eGFR_scr", Outcome="all_infections",MR_causal_effect = MR_CE,
                              MR_SE = MR_SE, MR_HR = MR_HR, MR_HR_lower95=MR_HR_lower95,MR_HR_upper95=MR_HR_upper95,
                              MR_pvalue=MR_pvalue,samples=samples,events=events)


#pneumonia
all_data_scr$follow_up_time = all_data_scr$pneumonia_time
all_data_scr$event.type = all_data_scr$pneumonia_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ GIV_eGFR_scr+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scr)
MR_CE = -summary(fit2)$coefficients[1,1]
MR_SE = summary(fit2)$coefficients[1,3]
MR_HR=1/summary(fit2)$coefficients[1,2]
MR_pvalue=summary(fit2)$coefficients[1,5]
MR_HR_lower95=1/summary(fit2)$conf.int[1,4]
MR_HR_upper95=1/summary(fit2)$conf.int[1,3]
samples=summary(fit2)$n
events=summary(fit2)$nevent
Cox_R2=summary(fit2)$rsq[[1]]
Cox_F=(summary(fit2)$rsq[[1]]/13)/((1-summary(fit2)$rsq[[1]])/(summary(fit2)$n-13-1))

scr_pneumonia = data.frame(Exposure="eGFR_scr", Outcome="pneumonia",MR_causal_effect = MR_CE,
                              MR_SE = MR_SE, MR_HR = MR_HR, MR_HR_lower95=MR_HR_lower95,MR_HR_upper95=MR_HR_upper95,
                              MR_pvalue=MR_pvalue,samples=samples,events=events)


#sepsis
all_data_scr$follow_up_time = all_data_scr$sepsis_time
all_data_scr$event.type = all_data_scr$sepsis_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ GIV_eGFR_scr+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scr)
MR_CE = -summary(fit2)$coefficients[1,1]
MR_SE = summary(fit2)$coefficients[1,3]
MR_HR=1/summary(fit2)$coefficients[1,2]
MR_pvalue=summary(fit2)$coefficients[1,5]
MR_HR_lower95=1/summary(fit2)$conf.int[1,4]
MR_HR_upper95=1/summary(fit2)$conf.int[1,3]
samples=summary(fit2)$n
events=summary(fit2)$nevent
Cox_R2=summary(fit2)$rsq[[1]]
Cox_F=(summary(fit2)$rsq[[1]]/13)/((1-summary(fit2)$rsq[[1]])/(summary(fit2)$n-13-1))

scr_sepsis = data.frame(Exposure="eGFR_scr", Outcome="sepsis",MR_causal_effect = MR_CE,
                           MR_SE = MR_SE, MR_HR = MR_HR, MR_HR_lower95=MR_HR_lower95,MR_HR_upper95=MR_HR_upper95,
                           MR_pvalue=MR_pvalue,samples=samples,events=events)


#skin_soft_tissue
all_data_scr$follow_up_time = all_data_scr$skin_soft_tissue_time
all_data_scr$event.type = all_data_scr$skin_soft_tissue_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ GIV_eGFR_scr+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scr)
MR_CE = -summary(fit2)$coefficients[1,1]
MR_SE = summary(fit2)$coefficients[1,3]
MR_HR=1/summary(fit2)$coefficients[1,2]
MR_pvalue=summary(fit2)$coefficients[1,5]
MR_HR_lower95=1/summary(fit2)$conf.int[1,4]
MR_HR_upper95=1/summary(fit2)$conf.int[1,3]
samples=summary(fit2)$n
events=summary(fit2)$nevent
Cox_R2=summary(fit2)$rsq[[1]]
Cox_F=(summary(fit2)$rsq[[1]]/13)/((1-summary(fit2)$rsq[[1]])/(summary(fit2)$n-13-1))

scr_skin_soft_tissue = data.frame(Exposure="eGFR_scr", Outcome="skin_soft_tissue",MR_causal_effect = MR_CE,
                           MR_SE = MR_SE, MR_HR = MR_HR, MR_HR_lower95=MR_HR_lower95,MR_HR_upper95=MR_HR_upper95,
                           MR_pvalue=MR_pvalue,samples=samples,events=events)


#skin_soft_tissue
all_data_scr$follow_up_time = all_data_scr$urinary_tract_time
all_data_scr$event.type = all_data_scr$urinary_tract_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ GIV_eGFR_scr+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scr)
MR_CE = -summary(fit2)$coefficients[1,1]
MR_SE = summary(fit2)$coefficients[1,3]
MR_HR=1/summary(fit2)$coefficients[1,2]
MR_pvalue=summary(fit2)$coefficients[1,5]
MR_HR_lower95=1/summary(fit2)$conf.int[1,4]
MR_HR_upper95=1/summary(fit2)$conf.int[1,3]
samples=summary(fit2)$n
events=summary(fit2)$nevent
Cox_R2=summary(fit2)$rsq[[1]]
Cox_F=(summary(fit2)$rsq[[1]]/13)/((1-summary(fit2)$rsq[[1]])/(summary(fit2)$n-13-1))

scr_urinary_tract = data.frame(Exposure="eGFR_scr", Outcome="urinary_tract",MR_causal_effect = MR_CE,
                                  MR_SE = MR_SE, MR_HR = MR_HR, MR_HR_lower95=MR_HR_lower95,MR_HR_upper95=MR_HR_upper95,
                                  MR_pvalue=MR_pvalue,samples=samples,events=events)

#eGFR_scys

#all_infection
all_data_scys$follow_up_time = all_data_scys$all_infection_time
all_data_scys$event.type = all_data_scys$all_infection_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ GIV_eGFR_scys+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scys)
MR_CE = -summary(fit2)$coefficients[1,1]
MR_SE = summary(fit2)$coefficients[1,3]
MR_HR=1/summary(fit2)$coefficients[1,2]
MR_pvalue=summary(fit2)$coefficients[1,5]
MR_HR_lower95=1/summary(fit2)$conf.int[1,4]
MR_HR_upper95=1/summary(fit2)$conf.int[1,3]
samples=summary(fit2)$n
events=summary(fit2)$nevent
Cox_R2=summary(fit2)$rsq[[1]]
Cox_F=(summary(fit2)$rsq[[1]]/13)/((1-summary(fit2)$rsq[[1]])/(summary(fit2)$n-13-1))

scys_allinfection = data.frame(Exposure="eGFR_scys", Outcome="all_infections",MR_causal_effect = MR_CE,
                              MR_SE = MR_SE, MR_HR = MR_HR, MR_HR_lower95=MR_HR_lower95,MR_HR_upper95=MR_HR_upper95,
                              MR_pvalue=MR_pvalue,samples=samples,events=events)


#pneumonia
all_data_scys$follow_up_time = all_data_scys$pneumonia_time
all_data_scys$event.type = all_data_scys$pneumonia_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ GIV_eGFR_scys+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scys)
MR_CE = -summary(fit2)$coefficients[1,1]
MR_SE = summary(fit2)$coefficients[1,3]
MR_HR=1/summary(fit2)$coefficients[1,2]
MR_pvalue=summary(fit2)$coefficients[1,5]
MR_HR_lower95=1/summary(fit2)$conf.int[1,4]
MR_HR_upper95=1/summary(fit2)$conf.int[1,3]
samples=summary(fit2)$n
events=summary(fit2)$nevent
Cox_R2=summary(fit2)$rsq[[1]]
Cox_F=(summary(fit2)$rsq[[1]]/13)/((1-summary(fit2)$rsq[[1]])/(summary(fit2)$n-13-1))

scys_pneumonia = data.frame(Exposure="eGFR_scys", Outcome="pneumonia",MR_causal_effect = MR_CE,
                           MR_SE = MR_SE, MR_HR = MR_HR, MR_HR_lower95=MR_HR_lower95,MR_HR_upper95=MR_HR_upper95,
                           MR_pvalue=MR_pvalue,samples=samples,events=events)


#sepsis
all_data_scys$follow_up_time = all_data_scys$sepsis_time
all_data_scys$event.type = all_data_scys$sepsis_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ GIV_eGFR_scys+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scys)
MR_CE = -summary(fit2)$coefficients[1,1]
MR_SE = summary(fit2)$coefficients[1,3]
MR_HR=1/summary(fit2)$coefficients[1,2]
MR_pvalue=summary(fit2)$coefficients[1,5]
MR_HR_lower95=1/summary(fit2)$conf.int[1,4]
MR_HR_upper95=1/summary(fit2)$conf.int[1,3]
samples=summary(fit2)$n
events=summary(fit2)$nevent
Cox_R2=summary(fit2)$rsq[[1]]
Cox_F=(summary(fit2)$rsq[[1]]/13)/((1-summary(fit2)$rsq[[1]])/(summary(fit2)$n-13-1))

scys_sepsis = data.frame(Exposure="eGFR_scys", Outcome="sepsis",MR_causal_effect = MR_CE,
                        MR_SE = MR_SE, MR_HR = MR_HR, MR_HR_lower95=MR_HR_lower95,MR_HR_upper95=MR_HR_upper95,
                        MR_pvalue=MR_pvalue,samples=samples,events=events)


#skin_soft_tissue
all_data_scys$follow_up_time = all_data_scys$skin_soft_tissue_time
all_data_scys$event.type = all_data_scys$skin_soft_tissue_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ GIV_eGFR_scys+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scys)
MR_CE = -summary(fit2)$coefficients[1,1]
MR_SE = summary(fit2)$coefficients[1,3]
MR_HR=1/summary(fit2)$coefficients[1,2]
MR_pvalue=summary(fit2)$coefficients[1,5]
MR_HR_lower95=1/summary(fit2)$conf.int[1,4]
MR_HR_upper95=1/summary(fit2)$conf.int[1,3]
samples=summary(fit2)$n
events=summary(fit2)$nevent
Cox_R2=summary(fit2)$rsq[[1]]
Cox_F=(summary(fit2)$rsq[[1]]/13)/((1-summary(fit2)$rsq[[1]])/(summary(fit2)$n-13-1))

scys_skin_soft_tissue = data.frame(Exposure="eGFR_scys", Outcome="skin_soft_tissue",MR_causal_effect = MR_CE,
                                  MR_SE = MR_SE, MR_HR = MR_HR, MR_HR_lower95=MR_HR_lower95,MR_HR_upper95=MR_HR_upper95,
                                  MR_pvalue=MR_pvalue,samples=samples,events=events)


#skin_soft_tissue
all_data_scys$follow_up_time = all_data_scys$urinary_tract_time
all_data_scys$event.type = all_data_scys$urinary_tract_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ GIV_eGFR_scys+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_scys)
MR_CE = -summary(fit2)$coefficients[1,1]
MR_SE = summary(fit2)$coefficients[1,3]
MR_HR=1/summary(fit2)$coefficients[1,2]
MR_pvalue=summary(fit2)$coefficients[1,5]
MR_HR_lower95=1/summary(fit2)$conf.int[1,4]
MR_HR_upper95=1/summary(fit2)$conf.int[1,3]
samples=summary(fit2)$n
events=summary(fit2)$nevent
Cox_R2=summary(fit2)$rsq[[1]]
Cox_F=(summary(fit2)$rsq[[1]]/13)/((1-summary(fit2)$rsq[[1]])/(summary(fit2)$n-13-1))

scys_urinary_tract = data.frame(Exposure="eGFR_scys", Outcome="urinary_tract",MR_causal_effect = MR_CE,
                               MR_SE = MR_SE, MR_HR = MR_HR, MR_HR_lower95=MR_HR_lower95,MR_HR_upper95=MR_HR_upper95,
                               MR_pvalue=MR_pvalue,samples=samples,events=events)



#BUN

#all_infection
all_data_BUN$follow_up_time = all_data_BUN$all_infection_time
all_data_BUN$event.type = all_data_BUN$all_infection_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ GIV_BUN+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_BUN)
MR_CE = summary(fit2)$coefficients[1,1]
MR_SE = summary(fit2)$coefficients[1,3]
MR_HR=summary(fit2)$coefficients[1,2]
MR_pvalue=summary(fit2)$coefficients[1,5]
MR_HR_lower95=summary(fit2)$conf.int[1,3]
MR_HR_upper95=summary(fit2)$conf.int[1,4]
samples=summary(fit2)$n
events=summary(fit2)$nevent
Cox_R2=summary(fit2)$rsq[[1]]
Cox_F=(summary(fit2)$rsq[[1]]/13)/((1-summary(fit2)$rsq[[1]])/(summary(fit2)$n-13-1))

BUN_allinfection = data.frame(Exposure="BUN", Outcome="all_infections",MR_causal_effect = MR_CE,
                              MR_SE = MR_SE, MR_HR = MR_HR, MR_HR_lower95=MR_HR_lower95,MR_HR_upper95=MR_HR_upper95,
                              MR_pvalue=MR_pvalue,samples=samples,events=events)


#pneumonia
all_data_BUN$follow_up_time = all_data_BUN$pneumonia_time
all_data_BUN$event.type = all_data_BUN$pneumonia_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ GIV_BUN+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_BUN)
MR_CE = summary(fit2)$coefficients[1,1]
MR_SE = summary(fit2)$coefficients[1,3]
MR_HR=summary(fit2)$coefficients[1,2]
MR_pvalue=summary(fit2)$coefficients[1,5]
MR_HR_lower95=summary(fit2)$conf.int[1,3]
MR_HR_upper95=summary(fit2)$conf.int[1,4]
samples=summary(fit2)$n
events=summary(fit2)$nevent
Cox_R2=summary(fit2)$rsq[[1]]
Cox_F=(summary(fit2)$rsq[[1]]/13)/((1-summary(fit2)$rsq[[1]])/(summary(fit2)$n-13-1))

BUN_pneumonia = data.frame(Exposure="BUN", Outcome="pneumonia",MR_causal_effect = MR_CE,
                           MR_SE = MR_SE, MR_HR = MR_HR, MR_HR_lower95=MR_HR_lower95,MR_HR_upper95=MR_HR_upper95,
                           MR_pvalue=MR_pvalue,samples=samples,events=events)


#sepsis
all_data_BUN$follow_up_time = all_data_BUN$sepsis_time
all_data_BUN$event.type = all_data_BUN$sepsis_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ GIV_BUN+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_BUN)
MR_CE = summary(fit2)$coefficients[1,1]
MR_SE = summary(fit2)$coefficients[1,3]
MR_HR=summary(fit2)$coefficients[1,2]
MR_pvalue=summary(fit2)$coefficients[1,5]
MR_HR_lower95=summary(fit2)$conf.int[1,3]
MR_HR_upper95=summary(fit2)$conf.int[1,4]
samples=summary(fit2)$n
events=summary(fit2)$nevent
Cox_R2=summary(fit2)$rsq[[1]]
Cox_F=(summary(fit2)$rsq[[1]]/13)/((1-summary(fit2)$rsq[[1]])/(summary(fit2)$n-13-1))

BUN_sepsis = data.frame(Exposure="BUN", Outcome="sepsis",MR_causal_effect = MR_CE,
                        MR_SE = MR_SE, MR_HR = MR_HR, MR_HR_lower95=MR_HR_lower95,MR_HR_upper95=MR_HR_upper95,
                        MR_pvalue=MR_pvalue,samples=samples,events=events)


#skin_soft_tissue
all_data_BUN$follow_up_time = all_data_BUN$skin_soft_tissue_time
all_data_BUN$event.type = all_data_BUN$skin_soft_tissue_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ GIV_BUN+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_BUN)
MR_CE = summary(fit2)$coefficients[1,1]
MR_SE = summary(fit2)$coefficients[1,3]
MR_HR=summary(fit2)$coefficients[1,2]
MR_pvalue=summary(fit2)$coefficients[1,5]
MR_HR_lower95=summary(fit2)$conf.int[1,3]
MR_HR_upper95=summary(fit2)$conf.int[1,4]
samples=summary(fit2)$n
events=summary(fit2)$nevent
Cox_R2=summary(fit2)$rsq[[1]]
Cox_F=(summary(fit2)$rsq[[1]]/13)/((1-summary(fit2)$rsq[[1]])/(summary(fit2)$n-13-1))

BUN_skin_soft_tissue = data.frame(Exposure="BUN", Outcome="skin_soft_tissue",MR_causal_effect = MR_CE,
                                  MR_SE = MR_SE, MR_HR = MR_HR, MR_HR_lower95=MR_HR_lower95,MR_HR_upper95=MR_HR_upper95,
                                  MR_pvalue=MR_pvalue,samples=samples,events=events)


#skin_soft_tissue
all_data_BUN$follow_up_time = all_data_BUN$urinary_tract_time
all_data_BUN$event.type = all_data_BUN$urinary_tract_event
fit2 = coxph(Surv(follow_up_time, event.type) ~ GIV_BUN+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = all_data_BUN)
MR_CE = summary(fit2)$coefficients[1,1]
MR_SE = summary(fit2)$coefficients[1,3]
MR_HR=summary(fit2)$coefficients[1,2]
MR_pvalue=summary(fit2)$coefficients[1,5]
MR_HR_lower95=summary(fit2)$conf.int[1,3]
MR_HR_upper95=summary(fit2)$conf.int[1,4]
samples=summary(fit2)$n
events=summary(fit2)$nevent
Cox_R2=summary(fit2)$rsq[[1]]
Cox_F=(summary(fit2)$rsq[[1]]/13)/((1-summary(fit2)$rsq[[1]])/(summary(fit2)$n-13-1))

BUN_urinary_tract = data.frame(Exposure="BUN", Outcome="urinary_tract",MR_causal_effect = MR_CE,
                               MR_SE = MR_SE, MR_HR = MR_HR, MR_HR_lower95=MR_HR_lower95,MR_HR_upper95=MR_HR_upper95,
                               MR_pvalue=MR_pvalue,samples=samples,events=events)


MR_result=rbind(scys_allinfection,scys_pneumonia,scys_sepsis,scys_skin_soft_tissue,scys_urinary_tract,
                scr_allinfection,scr_pneumonia,scr_sepsis,scr_skin_soft_tissue,scr_urinary_tract,
                BUN_allinfection,BUN_pneumonia,BUN_sepsis,BUN_skin_soft_tissue,BUN_urinary_tract)

write.csv(RF,"RF_new.csv")
write.csv(MR_result,"linear_MR_result_new.csv")

