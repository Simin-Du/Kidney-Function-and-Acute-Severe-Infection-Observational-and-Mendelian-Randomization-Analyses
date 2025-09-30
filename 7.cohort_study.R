library(dplyr)
library(survival)
library(cmprsk)
library(riskRegression)
library(car)
library(Publish)
library(Hmisc)
library(rms)
library(survminer)

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
all_data=read.table("all_data_imputation_cohort.txt", header=T)

all_data$education=factor(all_data$education,order=T)
all_data$smoking=factor(all_data$smoking,order=F)
all_data$alcohol=factor(all_data$alcohol,order=F)
all_data$UACR_cat=factor(all_data$UACR_cat,order=T)
all_data$eGFR_scr_cat=factor(all_data$eGFR_scr_cat,order=F)
all_data$eGFR_scys_cat=factor(all_data$eGFR_scys_cat,order=F)
all_data$BUN_cat=factor(all_data$BUN_cat,order=F)
all_data$eGFR_scr_10=all_data$eGFR_scr/10
all_data$eGFR_scys_10=all_data$eGFR_scys/10

all_data_eGFR_scr=all_data[which(!is.na(all_data$eGFR_scr)),]
all_data_eGFR_scys=all_data[which(!is.na(all_data$eGFR_scys)),]
all_data_BUN=all_data[which(!is.na(all_data$BUN)),]


#log_eGFR_scr
model1.1 <- as.formula(paste("Surv(all_infection_time, all_infection_event) ~  eGFR_scr_10+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))

sm1.1 <- coxph(as.formula(model1.1), data = all_data)
vif <- vif(sm1.1)
vif[which(vif>5)]
print("model 1.1")
publish(sm1.1)
prop_test1.1 <- cox.zph(sm1.1,transform="identity",terms=T)
prop_test1.1
row.names(prop_test1.1$table)[which(prop_test1.1$table[,3] < 0.05)]
PH_fail <- setdiff(rownames(prop_test1.1$table)[which(prop_test1.1$table[,3]<0.05)],"GLOBAL")
#pdf("model1.1.pdf",width = 25,height = 20)
#ggcoxzph(prop_test1.1,var = PH_fail,nsmo = dim(all_data)[1],point.alpha = 1,point.size = 0.25)
#dev.off()

res=summary(sm1.1)
coef=-round(res$coefficients[1,1],3)
HR=1/round(res$coefficients[1,2],3)
se=round(res$coefficients[1,3],3)
PValue=round(res$coefficients[1,5],3)
CI5=1/round(res$conf.int[1,4],3)
CI95=1/round(res$conf.int[1,3],3)
n=res$n
nevent=res$nevent
model1.1_result=data.frame("model"="model1.1","exposure"="eGFR_scr","outcome"="all_infection","coef"=coef,"HR"=HR,"CI5"=CI5,"CI95"=CI95,"se"=se,"p-value"=PValue,"number of participant"=n,"number of outcomes"=nevent,
                           "time_at_risk"=sum(all_data_eGFR_scr$all_infection_time/365),"incidence"=sum(all_data_eGFR_scr$all_infection_event)/sum(all_data_eGFR_scr$all_infection_time/365)*100000)


model1.2 <- as.formula(paste("Surv(pneumonia_time, pneumonia_event) ~  eGFR_scr_10+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm1.2 <- coxph(as.formula(model1.2), data = all_data)
vif <- vif(sm1.2)
vif[which(vif>5)]
print("model 1.2")
publish(sm1.2)
prop_test1.2 <- cox.zph(sm1.2,transform="identity",terms=T)
prop_test1.2
row.names(prop_test1.2$table)[which(prop_test1.2$table[,3] < 0.05)]
PH_fail <- setdiff(rownames(prop_test1.2$table)[which(prop_test1.2$table[,3]<0.05)],"GLOBAL")
#pdf("model1.2.pdf",width = 25,height = 20)
#ggcoxzph(prop_test1.2,var = PH_fail,nsmo = dim(all_data)[1],point.alpha = 1,point.size = 0.25)
#dev.off()

res=summary(sm1.2)
coef=-round(res$coefficients[1,1],3)
HR=1/round(res$coefficients[1,2],3)
se=round(res$coefficients[1,3],3)
PValue=round(res$coefficients[1,5],3)
CI5=1/round(res$conf.int[1,4],3)
CI95=1/round(res$conf.int[1,3],3)
n=res$n
nevent=res$nevent
model1.2_result=data.frame("model"="model1.2","exposure"="eGFR_scr","outcome"="Pneumonia","coef"=coef,"HR"=HR,"CI5"=CI5,"CI95"=CI95,"se"=se,"p-value"=PValue,"number of participant"=n,"number of outcomes"=nevent,
                           "time_at_risk"=sum(all_data_eGFR_scr$pneumonia_time/365),"incidence"=sum(all_data_eGFR_scr$pneumonia_event)/sum(all_data_eGFR_scr$pneumonia_time/365)*100000)


model1.3 <- as.formula(paste("Surv(sepsis_time, sepsis_event) ~  eGFR_scr_10+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm1.3 <- coxph(as.formula(model1.3), data = all_data)
vif <- vif(sm1.3)
vif[which(vif>5)]
print("model 1.3")
publish(sm1.3)
prop_test1.3 <- cox.zph(sm1.3,transform="identity",terms=T)
prop_test1.3
row.names(prop_test1.3$table)[which(prop_test1.3$table[,3] < 0.05)]
PH_fail <- setdiff(rownames(prop_test1.3$table)[which(prop_test1.3$table[,3]<0.05)],"GLOBAL")
#pdf("model1.3.pdf",width = 25,height = 20)
#ggcoxzph(prop_test1.3,var = PH_fail,nsmo = dim(all_data)[1],point.alpha = 1,point.size = 0.25)
#dev.off()

res=summary(sm1.3)
coef=-round(res$coefficients[1,1],3)
HR=1/round(res$coefficients[1,2],3)
se=round(res$coefficients[1,3],3)
PValue=round(res$coefficients[1,5],3)
CI5=1/round(res$conf.int[1,4],3)
CI95=1/round(res$conf.int[1,3],3)
n=res$n
nevent=res$nevent
model1.3_result=data.frame("model"="model1.3","exposure"="eGFR_scr","outcome"="Sepsis","coef"=coef,"HR"=HR,"CI5"=CI5,"CI95"=CI95,"se"=se,"p-value"=PValue,"number of participant"=n,"number of outcomes"=nevent,
                           "time_at_risk"=sum(all_data_eGFR_scr$sepsis_time/365),"incidence"=sum(all_data_eGFR_scr$sepsis_event)/sum(all_data_eGFR_scr$sepsis_time/365)*100000)


model1.4 <- as.formula(paste("Surv(skin_soft_tissue_time, skin_soft_tissue_event) ~  eGFR_scr_10+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm1.4 <- coxph(as.formula(model1.4), data = all_data)
vif <- vif(sm1.4)
vif[which(vif>5)]
print("model 1.4")
publish(sm1.4)
prop_test1.4 <- cox.zph(sm1.4,transform="identity",terms=T)
prop_test1.4
row.names(prop_test1.4$table)[which(prop_test1.4$table[,3] < 0.05)]
PH_fail <- setdiff(rownames(prop_test1.4$table)[which(prop_test1.4$table[,3]<0.05)],"GLOBAL")
#pdf("model1.4.pdf",width = 25,height = 20)
#ggcoxzph(prop_test1.4,var = PH_fail,nsmo = dim(all_data)[1],point.alpha = 1,point.size = 0.25)
#dev.off()

res=summary(sm1.4)
coef=-round(res$coefficients[1,1],3)
HR=1/round(res$coefficients[1,2],3)
se=round(res$coefficients[1,3],3)
PValue=round(res$coefficients[1,5],3)
CI5=1/round(res$conf.int[1,4],3)
CI95=1/round(res$conf.int[1,3],3)
n=res$n
nevent=res$nevent
model1.4_result=data.frame("model"="model1.4","exposure"="eGFR_scr","outcome"="skin_soft_tissue","coef"=coef,"HR"=HR,"CI5"=CI5,"CI95"=CI95,"se"=se,"p-value"=PValue,"number of participant"=n,"number of outcomes"=nevent,
                           "time_at_risk"=sum(all_data_eGFR_scr$skin_soft_tissue_time/365),"incidence"=sum(all_data_eGFR_scr$skin_soft_tissue_event)/sum(all_data_eGFR_scr$skin_soft_tissue_time/365)*100000)

all_data$urinary_tract_time
model1.5 <- as.formula(paste("Surv(urinary_tract_time, urinary_tract_event) ~  eGFR_scr_10+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm1.5 <- coxph(as.formula(model1.5), data = all_data)
vif <- vif(sm1.5)
vif[which(vif>5)]
print("model 1.5")
publish(sm1.5)
prop_test1.5 <- cox.zph(sm1.5,transform="identity",terms=T)
prop_test1.5
row.names(prop_test1.5$table)[which(prop_test1.5$table[,3] < 0.05)]
PH_fail <- setdiff(rownames(prop_test1.5$table)[which(prop_test1.5$table[,3]<0.05)],"GLOBAL")
#pdf("model1.5.pdf",width = 25,height = 20)
#ggcoxzph(prop_test1.5,var = PH_fail,nsmo = dim(all_data)[1],point.alpha = 1,point.size = 0.25)
#dev.off()

res=summary(sm1.5)
coef=-round(res$coefficients[1,1],3)
HR=1/round(res$coefficients[1,2],3)
se=round(res$coefficients[1,3],3)
PValue=round(res$coefficients[1,5],3)
CI5=1/round(res$conf.int[1,4],3)
CI95=1/round(res$conf.int[1,3],3)
n=res$n
nevent=res$nevent
model1.5_result=data.frame("model"="model1.5","exposure"="eGFR_scr","outcome"="urinary_tract","coef"=coef,"HR"=HR,"CI5"=CI5,"CI95"=CI95,"se"=se,"p-value"=PValue,"number of participant"=n,"number of outcomes"=nevent,
                           "time_at_risk"=sum(all_data_eGFR_scr$urinary_tract_time/365),"incidence"=sum(all_data_eGFR_scr$urinary_tract_event)/sum(all_data_eGFR_scr$urinary_tract_time/365)*100000)

#log_eGFR_scys
model2.1 <- as.formula(paste("Surv(all_infection_time, all_infection_event) ~  eGFR_scys_10+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm2.1 <- coxph(as.formula(model2.1), data = all_data)
vif <- vif(sm2.1)
vif[which(vif>5)]
print("model 2.1")
publish(sm2.1)
prop_test2.1 <- cox.zph(sm2.1,transform="identity",terms=T)
prop_test2.1
row.names(prop_test2.1$table)[which(prop_test2.1$table[,3] < 0.05)]
PH_fail <- setdiff(rownames(prop_test2.1$table)[which(prop_test2.1$table[,3]<0.05)],"GLOBAL")
#pdf("model2.1.pdf",width = 25,height = 20)
#ggcoxzph(prop_test2.1,var = PH_fail,nsmo = dim(all_data)[1],point.alpha = 1,point.size = 0.25)
#dev.off()

res=summary(sm2.1)
coef=-round(res$coefficients[1,1],3)
HR=1/round(res$coefficients[1,2],3)
se=round(res$coefficients[1,3],3)
PValue=round(res$coefficients[1,5],3)
CI5=1/round(res$conf.int[1,4],3)
CI95=1/round(res$conf.int[1,3],3)
n=res$n
nevent=res$nevent
model2.1_result=data.frame("model"="model2.1","exposure"="eGFR_scys","outcome"="all_infection","coef"=coef,"HR"=HR,"CI5"=CI5,"CI95"=CI95,"se"=se,"p-value"=PValue,"number of participant"=n,"number of outcomes"=nevent,
                           "time_at_risk"=sum(all_data_eGFR_scys$all_infection_time/365),"incidence"=sum(all_data_eGFR_scys$all_infection_event)/sum(all_data_eGFR_scys$all_infection_time/365)*100000)



model2.2 <- as.formula(paste("Surv(pneumonia_time, pneumonia_event) ~  eGFR_scys_10+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm2.2 <- coxph(as.formula(model2.2), data = all_data)
vif <- vif(sm2.2)
vif[which(vif>5)]
print("model 2.2")
publish(sm2.2)
prop_test2.2 <- cox.zph(sm2.2,transform="identity",terms=T)
prop_test2.2
row.names(prop_test2.2$table)[which(prop_test2.2$table[,3] < 0.05)]
PH_fail <- setdiff(rownames(prop_test2.2$table)[which(prop_test2.2$table[,3]<0.05)],"GLOBAL")
#pdf("model2.2.pdf",width = 25,height = 20)
#ggcoxzph(prop_test2.2,var = PH_fail,nsmo = dim(all_data)[1],point.alpha = 1,point.size = 0.25)
#dev.off()

res=summary(sm2.2)
coef=-round(res$coefficients[1,1],3)
HR=1/round(res$coefficients[1,2],3)
se=round(res$coefficients[1,3],3)
PValue=round(res$coefficients[1,5],3)
CI5=1/round(res$conf.int[1,4],3)
CI95=1/round(res$conf.int[1,3],3)
n=res$n
nevent=res$nevent
model2.2_result=data.frame("model"="model2.2","exposure"="eGFR_scys","outcome"="Pneumonia","coef"=coef,"HR"=HR,"CI5"=CI5,"CI95"=CI95,"se"=se,"p-value"=PValue,"number of participant"=n,"number of outcomes"=nevent,
                           "time_at_risk"=sum(all_data_eGFR_scys$pneumonia_time/365),"incidence"=sum(all_data_eGFR_scys$pneumonia_event)/sum(all_data_eGFR_scys$pneumonia_time/365)*100000)


model2.3 <- as.formula(paste("Surv(sepsis_time, sepsis_event) ~  eGFR_scys_10+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm2.3 <- coxph(as.formula(model2.3), data = all_data)
vif <- vif(sm2.3)
vif[which(vif>5)]
print("model 2.3")
publish(sm2.3)
prop_test2.3 <- cox.zph(sm2.3,transform="identity",terms=T)
prop_test2.3
row.names(prop_test2.3$table)[which(prop_test2.3$table[,3] < 0.05)]
PH_fail <- setdiff(rownames(prop_test2.3$table)[which(prop_test2.3$table[,3]<0.05)],"GLOBAL")
#pdf("model2.3.pdf",width = 25,height = 20)
#ggcoxzph(prop_test2.3,var = PH_fail,nsmo = dim(all_data)[1],point.alpha = 1,point.size = 0.25)
#dev.off()

res=summary(sm2.3)
coef=-round(res$coefficients[1,1],3)
HR=1/round(res$coefficients[1,2],3)
se=round(res$coefficients[1,3],3)
PValue=round(res$coefficients[1,5],3)
CI5=1/round(res$conf.int[1,4],3)
CI95=1/round(res$conf.int[1,3],3)
n=res$n
nevent=res$nevent
model2.3_result=data.frame("model"="model2.3","exposure"="eGFR_scys","outcome"="Sepsis","coef"=coef,"HR"=HR,"CI5"=CI5,"CI95"=CI95,"se"=se,"p-value"=PValue,"number of participant"=n,"number of outcomes"=nevent,
                           "time_at_risk"=sum(all_data_eGFR_scys$sepsis_time/365),"incidence"=sum(all_data_eGFR_scys$sepsis_event)/sum(all_data_eGFR_scys$sepsis_time/365)*100000)


model2.4 <- as.formula(paste("Surv(skin_soft_tissue_time, skin_soft_tissue_event) ~  eGFR_scys_10+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm2.4 <- coxph(as.formula(model2.4), data = all_data)
vif <- vif(sm2.4)
vif[which(vif>5)]
print("model 2.4")
publish(sm2.4)
prop_test2.4 <- cox.zph(sm2.4,transform="identity",terms=T)
prop_test2.4
row.names(prop_test2.4$table)[which(prop_test2.4$table[,3] < 0.05)]
PH_fail <- setdiff(rownames(prop_test2.4$table)[which(prop_test2.4$table[,3]<0.05)],"GLOBAL")
#pdf("model2.4.pdf",width = 25,height = 20)
#ggcoxzph(prop_test2.4,var = PH_fail,nsmo = dim(all_data)[1],point.alpha = 1,point.size = 0.25)
#dev.off()

res=summary(sm2.4)
coef=-round(res$coefficients[1,1],3)
HR=1/round(res$coefficients[1,2],3)
se=round(res$coefficients[1,3],3)
PValue=round(res$coefficients[1,5],3)
CI5=1/round(res$conf.int[1,4],3)
CI95=1/round(res$conf.int[1,3],3)
n=res$n
nevent=res$nevent
model2.4_result=data.frame("model"="model2.4","exposure"="eGFR_scys","outcome"="skin_soft_tissue","coef"=coef,"HR"=HR,"CI5"=CI5,"CI95"=CI95,"se"=se,"p-value"=PValue,"number of participant"=n,"number of outcomes"=nevent,
                           "time_at_risk"=sum(all_data_eGFR_scys$skin_soft_tissue_time/365),"incidence"=sum(all_data_eGFR_scys$skin_soft_tissue_event)/sum(all_data_eGFR_scys$skin_soft_tissue_time/365)*100000)


model2.5 <- as.formula(paste("Surv(urinary_tract_time, urinary_tract_event) ~  eGFR_scys_10+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm2.5 <- coxph(as.formula(model2.5), data = all_data)
vif <- vif(sm2.5)
vif[which(vif>5)]
print("model 2.5")
publish(sm2.5)
prop_test2.5 <- cox.zph(sm2.5,transform="identity",terms=T)
prop_test2.5
row.names(prop_test2.5$table)[which(prop_test2.5$table[,3] < 0.05)]
PH_fail <- setdiff(rownames(prop_test2.5$table)[which(prop_test2.5$table[,3]<0.05)],"GLOBAL")
#pdf("model2.5.pdf",width = 25,height = 20)
#ggcoxzph(prop_test2.5,var = PH_fail,nsmo = dim(all_data)[1],point.alpha = 1,point.size = 0.25)
#dev.off()

res=summary(sm2.5)
coef=-round(res$coefficients[1,1],3)
HR=1/round(res$coefficients[1,2],3)
se=round(res$coefficients[1,3],3)
PValue=round(res$coefficients[1,5],3)
CI5=1/round(res$conf.int[1,4],3)
CI95=1/round(res$conf.int[1,3],3)
n=res$n
nevent=res$nevent
model2.5_result=data.frame("model"="model2.5","exposure"="eGFR_scys","outcome"="urinary_tract","coef"=coef,"HR"=HR,"CI5"=CI5,"CI95"=CI95,"se"=se,"p-value"=PValue,"number of participant"=n,"number of outcomes"=nevent,
                           "time_at_risk"=sum(all_data_eGFR_scys$urinary_tract_time/365),"incidence"=sum(all_data_eGFR_scys$urinary_tract_event)/sum(all_data_eGFR_scys$urinary_tract_time/365)*100000)


#BUN
model3.1 <- as.formula(paste("Surv(all_infection_time, all_infection_event) ~  BUN+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm3.1 <- coxph(as.formula(model3.1), data = all_data)
vif <- vif(sm3.1)
vif[which(vif>5)]
print("model 3.1")
publish(sm3.1)
prop_test3.1 <- cox.zph(sm3.1,transform="identity",terms=T)
prop_test3.1
row.names(prop_test3.1$table)[which(prop_test3.1$table[,3] < 0.05)]
PH_fail <- setdiff(rownames(prop_test3.1$table)[which(prop_test3.1$table[,3]<0.05)],"GLOBAL")
#pdf("model3.1.pdf",width = 25,height = 20)
#ggcoxzph(prop_test3.1,var = PH_fail,nsmo = dim(all_data)[1],point.alpha = 1,point.size = 0.25)
#dev.off()

res=summary(sm3.1)
coef=round(res$coefficients[1,1],3)
HR=round(res$coefficients[1,2],3)
se=round(res$coefficients[1,3],3)
PValue=round(res$coefficients[1,5],3)
CI5=round(res$conf.int[1,3],3)
CI95=round(res$conf.int[1,4],3)
n=res$n
nevent=res$nevent
model3.1_result=data.frame("model"="model3.1","exposure"="BUN","outcome"="all_infection","coef"=coef,"HR"=HR,"CI5"=CI5,"CI95"=CI95,"se"=se,"p-value"=PValue,"number of participant"=n,"number of outcomes"=nevent,
                           "time_at_risk"=sum(all_data_BUN$all_infection_time/365),"incidence"=sum(all_data_BUN$all_infection_event)/sum(all_data_BUN$all_infection_time/365)*100000)



model3.2 <- as.formula(paste("Surv(pneumonia_time, pneumonia_event) ~  BUN+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm3.2 <- coxph(as.formula(model3.2), data = all_data)
vif <- vif(sm3.2)
vif[which(vif>5)]
print("model 3.2")
publish(sm3.2)
prop_test3.2 <- cox.zph(sm3.2,transform="identity",terms=T)
prop_test3.2
row.names(prop_test3.2$table)[which(prop_test3.2$table[,3] < 0.05)]
PH_fail <- setdiff(rownames(prop_test3.2$table)[which(prop_test3.2$table[,3]<0.05)],"GLOBAL")
#pdf("model3.2.pdf",width = 25,height = 20)
#ggcoxzph(prop_test3.2,var = PH_fail,nsmo = dim(all_data)[1],point.alpha = 1,point.size = 0.25)
#dev.off()

res=summary(sm3.2)
coef=round(res$coefficients[1,1],3)
HR=round(res$coefficients[1,2],3)
se=round(res$coefficients[1,3],3)
PValue=round(res$coefficients[1,5],3)
CI5=round(res$conf.int[1,3],3)
CI95=round(res$conf.int[1,4],3)
n=res$n
nevent=res$nevent
model3.2_result=data.frame("model"="model3.2","exposure"="BUN","outcome"="Pneumonia","coef"=coef,"HR"=HR,"CI5"=CI5,"CI95"=CI95,"se"=se,"p-value"=PValue,"number of participant"=n,"number of outcomes"=nevent,
                           "time_at_risk"=sum(all_data_BUN$pneumonia_time/365),"incidence"=sum(all_data_BUN$pneumonia_event)/sum(all_data_BUN$pneumonia_time/365)*100000)


model3.3 <- as.formula(paste("Surv(sepsis_time, sepsis_event) ~  BUN+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm3.3 <- coxph(as.formula(model3.3), data = all_data)
vif <- vif(sm3.3)
vif[which(vif>5)]
print("model 3.3")
publish(sm3.3)
prop_test3.3 <- cox.zph(sm3.3,transform="identity",terms=T)
prop_test3.3
row.names(prop_test3.3$table)[which(prop_test3.3$table[,3] < 0.05)]
PH_fail <- setdiff(rownames(prop_test3.3$table)[which(prop_test3.3$table[,3]<0.05)],"GLOBAL")
#pdf("model2.3.pdf",width = 25,height = 20)
#ggcoxzph(prop_test3.3,var = PH_fail,nsmo = dim(all_data)[1],point.alpha = 1,point.size = 0.25)
#dev.off()

res=summary(sm3.3)
coef=round(res$coefficients[1,1],3)
HR=round(res$coefficients[1,2],3)
se=round(res$coefficients[1,3],3)
PValue=round(res$coefficients[1,5],3)
CI5=round(res$conf.int[1,3],3)
CI95=round(res$conf.int[1,4],3)
n=res$n
nevent=res$nevent
model3.3_result=data.frame("model"="model3.3","exposure"="BUN","outcome"="Sepsis","coef"=coef,"HR"=HR,"CI5"=CI5,"CI95"=CI95,"se"=se,"p-value"=PValue,"number of participant"=n,"number of outcomes"=nevent,
                           "time_at_risk"=sum(all_data_BUN$sepsis_time/365),"incidence"=sum(all_data_BUN$sepsis_event)/sum(all_data_BUN$sepsis_time/365)*100000)


model3.4 <- as.formula(paste("Surv(skin_soft_tissue_time, skin_soft_tissue_event) ~  BUN+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm3.4 <- coxph(as.formula(model3.4), data = all_data)
vif <- vif(sm3.4)
vif[which(vif>5)]
print("model 3.4")
publish(sm3.4)
prop_test3.4 <- cox.zph(sm3.4,transform="identity",terms=T)
prop_test3.4
row.names(prop_test3.4$table)[which(prop_test3.4$table[,3] < 0.05)]
PH_fail <- setdiff(rownames(prop_test3.4$table)[which(prop_test3.4$table[,3]<0.05)],"GLOBAL")
#pdf("model3.4.pdf",width = 25,height = 20)
#ggcoxzph(prop_test3.4,var = PH_fail,nsmo = dim(all_data)[1],point.alpha = 1,point.size = 0.25)
#dev.off()

res=summary(sm3.4)
coef=round(res$coefficients[1,1],3)
HR=round(res$coefficients[1,2],3)
se=round(res$coefficients[1,3],3)
PValue=round(res$coefficients[1,5],3)
CI5=round(res$conf.int[1,3],3)
CI95=round(res$conf.int[1,4],3)
n=res$n
nevent=res$nevent
model3.4_result=data.frame("model"="model3.4","exposure"="BUN","outcome"="skin_soft_tissue","coef"=coef,"HR"=HR,"CI5"=CI5,"CI95"=CI95,"se"=se,"p-value"=PValue,"number of participant"=n,"number of outcomes"=nevent,
                           "time_at_risk"=sum(all_data_BUN$skin_soft_tissue_time/365),"incidence"=sum(all_data_BUN$skin_soft_tissue_event)/sum(all_data_BUN$skin_soft_tissue_time/365)*100000)


model3.5 <- as.formula(paste("Surv(urinary_tract_time, urinary_tract_event) ~  BUN+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm3.5 <- coxph(as.formula(model3.5), data = all_data)
vif <- vif(sm3.5)
vif[which(vif>5)]
print("model 3.5")
publish(sm3.5)
prop_test3.5 <- cox.zph(sm3.5,transform="identity",terms=T)
prop_test3.5
row.names(prop_test3.5$table)[which(prop_test3.5$table[,3] < 0.05)]
PH_fail <- setdiff(rownames(prop_test3.5$table)[which(prop_test3.5$table[,3]<0.05)],"GLOBAL")
#pdf("model3.5.pdf",width = 25,height = 20)
#ggcoxzph(prop_test3.5,var = PH_fail,nsmo = dim(all_data)[1],point.alpha = 1,point.size = 0.25)
#dev.off()

res=summary(sm3.5)
coef=round(res$coefficients[1,1],3)
HR=round(res$coefficients[1,2],3)
se=round(res$coefficients[1,3],3)
PValue=round(res$coefficients[1,5],3)
CI5=round(res$conf.int[1,3],3)
CI95=round(res$conf.int[1,4],3)
n=res$n
nevent=res$nevent
model3.5_result=data.frame("model"="model3.5","exposure"="BUN","outcome"="urinary_tract","coef"=coef,"HR"=HR,"CI5"=CI5,"CI95"=CI95,"se"=se,"p-value"=PValue,"number of participant"=n,"number of outcomes"=nevent,
                           "time_at_risk"=sum(all_data_BUN$urinary_tract_time/365),"incidence"=sum(all_data_BUN$urinary_tract_event)/sum(all_data_BUN$urinary_tract_time/365)*100000)

model_result=rbind(model1.1_result,model1.2_result,model1.3_result,model1.4_result,model1.5_result,model2.1_result,
                   model2.2_result,model2.3_result,model2.4_result,model2.5_result,model3.1_result,
                   model3.2_result,model3.3_result,model3.4_result,model3.5_result)




write.csv(model_result,"cohort_study.csv")


#category
get_time_risk<-function(exposure,outcome_time){
  new_data=all_data[which(!is.na(all_data[,exposure])),]
  num_level=length(table(new_data[,exposure]))
  time_at_risk=c()
  for(i in 1:num_level){
    time_at_risk=c(time_at_risk,sum(new_data[which(new_data[,exposure]==i),outcome_time]/365))
  }
  return(time_at_risk)
}

get_incidence<-function(exposure,outcome_event,outcome_time){
  new_data=all_data[which(!is.na(all_data[,exposure])),]
  num_level=length(table(new_data[,exposure]))
  incidence=c()
  for(i in 1:num_level){
    incidence=c(incidence,sum(new_data[which(new_data[,exposure]==i),outcome_event])/sum(new_data[which(new_data[,exposure]==i),outcome_time]/365)*100000)
  }
  return(incidence)
}


model1.1 <- as.formula(paste("Surv(all_infection_time, all_infection_event) ~  eGFR_scr_cat+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm1.1 <- coxph(as.formula(model1.1), data = all_data)
res=summary(sm1.1)

coef=round(res$coefficients[1:4,1],3)
HR=round(res$coefficients[1:4,2],3)
se=round(res$coefficients[1:4,3],3)
PValue=round(res$coefficients[1:4,5],3)
CI5=round(res$conf.int[1:4,3],3)
CI95=round(res$conf.int[1:4,4],3)
model1.1_result=data.frame("model"=c("model1.1","model1.1","model1.1","model1.1","model1.1"),
                           "exposure"=c("eGFR_scr(90-105)","eGFR_scr(<60)","eGFR_scr(60-75)","eGFR_scr(75-90)","eGFR_scr(>=105)"),
                           "outcome"=c("all_infection","all_infection","all_infection","all_infection","all_infection"),
                           "coef"=c(NA,coef),"HR"=c(NA,HR),"CI5"=c(NA,CI5),"CI95"=c(NA,CI95),"se"=c(NA,se),"p-value"=c(NA,PValue),
                           "number of participant"=table(all_data$eGFR_scr_cat)[1:5],
                           "number of outcomes"=table(all_data$all_infection_event,all_data$eGFR_scr_cat)[2,1:5],
                           "time at risk"=get_time_risk("eGFR_scr_cat","all_infection_time"),
                           "incidence"=get_incidence("eGFR_scr_cat","all_infection_event","all_infection_time"))



model1.2 <- as.formula(paste("Surv(pneumonia_time, pneumonia_event) ~  eGFR_scr_cat+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm1.2 <- coxph(as.formula(model1.2), data = all_data)

res=summary(sm1.2)
coef=round(res$coefficients[1:4,1],3)
HR=round(res$coefficients[1:4,2],3)
se=round(res$coefficients[1:4,3],3)
PValue=round(res$coefficients[1:4,5],3)
CI5=round(res$conf.int[1:4,3],3)
CI95=round(res$conf.int[1:4,4],3)
model1.2_result=data.frame("model"=c("model1.2","model1.2","model1.2","model1.2","model1.2"),
                           "exposure"=c("eGFR_scr(90-105)","eGFR_scr(<60)","eGFR_scr(60-75)","eGFR_scr(75-90)","eGFR_scr(>=105)"),
                           "outcome"=c("Pneumonia","Pneumonia","Pneumonia","Pneumonia","Pneumonia"),
                           "coef"=c(NA,coef),"HR"=c(NA,HR),"CI5"=c(NA,CI5),"CI95"=c(NA,CI95),"se"=c(NA,se),"p-value"=c(NA,PValue),
                           "number of participant"=table(all_data$eGFR_scr_cat)[1:5],
                           "number of outcomes"=table(all_data$pneumonia_event,all_data$eGFR_scr_cat)[2,1:5],
                           "time at risk"=get_time_risk("eGFR_scr_cat","pneumonia_time"),
                           "incidence"=get_incidence("eGFR_scr_cat","pneumonia_event","pneumonia_time"))


model1.3 <- as.formula(paste("Surv(sepsis_time, sepsis_event) ~  eGFR_scr_cat+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm1.3 <- coxph(as.formula(model1.3), data = all_data)

res=summary(sm1.3)
coef=round(res$coefficients[1:4,1],3)
HR=round(res$coefficients[1:4,2],3)
se=round(res$coefficients[1:4,3],3)
PValue=round(res$coefficients[1:4,5],3)
CI5=round(res$conf.int[1:4,3],3)
CI95=round(res$conf.int[1:4,4],3)
model1.3_result=data.frame("model"=c("model1.3","model1.3","model1.3","model1.3","model1.3"),
                           "exposure"=c("eGFR_scr(90-105)","eGFR_scr(<60)","eGFR_scr(60-75)","eGFR_scr(75-90)","eGFR_scr(>=105)"),
                           "outcome"=c("Sepsis","Sepsis","Sepsis","Sepsis","Sepsis"),
                           "coef"=c(NA,coef),"HR"=c(NA,HR),"CI5"=c(NA,CI5),"CI95"=c(NA,CI95),"se"=c(NA,se),"p-value"=c(NA,PValue),
                           "number of participant"=table(all_data$eGFR_scr_cat)[1:5],
                           "number of outcomes"=table(all_data$sepsis_event,all_data$eGFR_scr_cat)[2,1:5],
                           "time at risk"=get_time_risk("eGFR_scr_cat","sepsis_time"),
                           "incidence"=get_incidence("eGFR_scr_cat","sepsis_event","sepsis_time"))



model1.4 <- as.formula(paste("Surv(skin_soft_tissue_time, skin_soft_tissue_event) ~  eGFR_scr_cat+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm1.4 <- coxph(as.formula(model1.4), data = all_data)

res=summary(sm1.4)
coef=round(res$coefficients[1:4,1],3)
HR=round(res$coefficients[1:4,2],3)
se=round(res$coefficients[1:4,3],3)
PValue=round(res$coefficients[1:4,5],3)
CI5=round(res$conf.int[1:4,3],3)
CI95=round(res$conf.int[1:4,4],3)
model1.4_result=data.frame("model"=c("model1.4","model1.4","model1.4","model1.4","model1.4"),
                           "exposure"=c("eGFR_scr(90-105)","eGFR_scr(<60)","eGFR_scr(60-75)","eGFR_scr(75-90)","eGFR_scr(>=105)"),
                           "outcome"=c("skin_soft_tissue","skin_soft_tissue","skin_soft_tissue","skin_soft_tissue","skin_soft_tissue"),
                           "coef"=c(NA,coef),"HR"=c(NA,HR),"CI5"=c(NA,CI5),"CI95"=c(NA,CI95),"se"=c(NA,se),"p-value"=c(NA,PValue),
                           "number of participant"=table(all_data$eGFR_scr_cat)[1:5],
                           "number of outcomes"=table(all_data$skin_soft_tissue_event,all_data$eGFR_scr_cat)[2,1:5],
                           "time at risk"=get_time_risk("eGFR_scr_cat","skin_soft_tissue_time"),
                           "incidence"=get_incidence("eGFR_scr_cat","skin_soft_tissue_event","skin_soft_tissue_time"))


model1.5 <- as.formula(paste("Surv(urinary_tract_time, urinary_tract_event) ~  eGFR_scr_cat+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm1.5 <- coxph(as.formula(model1.5), data = all_data)

res=summary(sm1.5)
coef=round(res$coefficients[1:4,1],3)
HR=round(res$coefficients[1:4,2],3)
se=round(res$coefficients[1:4,3],3)
PValue=round(res$coefficients[1:4,5],3)
CI5=round(res$conf.int[1:4,3],3)
CI95=round(res$conf.int[1:4,4],3)
model1.5_result=data.frame("model"=c("model1.5","model1.5","model1.5","model1.5","model1.5"),
                           "exposure"=c("eGFR_scr(90-105)","eGFR_scr(<60)","eGFR_scr(60-75)","eGFR_scr(75-90)","eGFR_scr(>=105)"),
                           "outcome"=c("urinary_tract","urinary_tract","urinary_tract","urinary_tract","urinary_tract"),
                           "coef"=c(NA,coef),"HR"=c(NA,HR),"CI5"=c(NA,CI5),"CI95"=c(NA,CI95),"se"=c(NA,se),"p-value"=c(NA,PValue),
                           "number of participant"=table(all_data$eGFR_scr_cat)[1:5],
                           "number of outcomes"=table(all_data$urinary_tract_event,all_data$eGFR_scr_cat)[2,1:5],
                           "time at risk"=get_time_risk("eGFR_scr_cat","urinary_tract_time"),
                           "incidence"=get_incidence("eGFR_scr_cat","urinary_tract_event","urinary_tract_time"))



#log_eGFR_scys
model2.1 <- as.formula(paste("Surv(all_infection_time, all_infection_event) ~  eGFR_scys_cat+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm2.1 <- coxph(as.formula(model2.1), data = all_data)

res=summary(sm2.1)
coef=round(res$coefficients[1:4,1],3)
HR=round(res$coefficients[1:4,2],3)
se=round(res$coefficients[1:4,3],3)
PValue=round(res$coefficients[1:4,5],3)
CI5=round(res$conf.int[1:4,3],3)
CI95=round(res$conf.int[1:4,4],3)
model2.1_result=data.frame("model"=c("model2.1","model2.1","model2.1","model2.1","model2.1"),
                           "exposure"=c("eGFR_scys(90-105)","eGFR_scys(<60)","eGFR_scys(60-75)","eGFR_scys(75-90)","eGFR_scys(>=105)"),
                           "outcome"=c("all_infection","all_infection","all_infection","all_infection","all_infection"),
                           "coef"=c(NA,coef),"HR"=c(NA,HR),"CI5"=c(NA,CI5),"CI95"=c(NA,CI95),"se"=c(NA,se),"p-value"=c(NA,PValue),
                           "number of participant"=table(all_data$eGFR_scys_cat)[1:5],
                           "number of outcomes"=table(all_data$all_infection_event,all_data$eGFR_scys_cat)[2,1:5],
                           "time at risk"=get_time_risk("eGFR_scys_cat","all_infection_time"),
                           "incidence"=get_incidence("eGFR_scys_cat","all_infection_event","all_infection_time"))


model2.2 <- as.formula(paste("Surv(pneumonia_time, pneumonia_event) ~  eGFR_scys_cat+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm2.2 <- coxph(as.formula(model2.2), data = all_data)

res=summary(sm2.2)
coef=round(res$coefficients[1:4,1],3)
HR=round(res$coefficients[1:4,2],3)
se=round(res$coefficients[1:4,3],3)
PValue=round(res$coefficients[1:4,5],3)
CI5=round(res$conf.int[1:4,3],3)
CI95=round(res$conf.int[1:4,4],3)
model2.2_result=data.frame("model"=c("model2.2","model2.2","model2.2","model2.2","model2.2"),
                           "exposure"=c("eGFR_scys(90-105)","eGFR_scys(<60)","eGFR_scys(60-75)","eGFR_scys(75-90)","eGFR_scys(>=105)"),
                           "outcome"=c("Pneumonia","Pneumonia","Pneumonia","Pneumonia","Pneumonia"),
                           "coef"=c(NA,coef),"HR"=c(NA,HR),"CI5"=c(NA,CI5),"CI95"=c(NA,CI95),"se"=c(NA,se),"p-value"=c(NA,PValue),
                           "number of participant"=table(all_data$eGFR_scys_cat)[1:5],
                           "number of outcomes"=table(all_data$pneumonia_event,all_data$eGFR_scys_cat)[2,1:5],
                           "time at risk"=get_time_risk("eGFR_scys_cat","pneumonia_time"),
                           "incidence"=get_incidence("eGFR_scys_cat","pneumonia_event","pneumonia_time"))


model2.3 <- as.formula(paste("Surv(sepsis_time, sepsis_event) ~  eGFR_scys_cat+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm2.3 <- coxph(as.formula(model2.3), data = all_data)

res=summary(sm2.3)
coef=round(res$coefficients[1:4,1],3)
HR=round(res$coefficients[1:4,2],3)
se=round(res$coefficients[1:4,3],3)
PValue=round(res$coefficients[1:4,5],3)
CI5=round(res$conf.int[1:4,3],3)
CI95=round(res$conf.int[1:4,4],3)
model2.3_result=data.frame("model"=c("model2.3","model2.3","model2.3","model2.3","model2.3"),
                           "exposure"=c("eGFR_scys(90-105)","eGFR_scys(<60)","eGFR_scys(60-75)","eGFR_scys(75-90)","eGFR_scys(>=105)"),
                           "outcome"=c("Sepsis","Sepsis","Sepsis","Sepsis","Sepsis"),
                           "coef"=c(NA,coef),"HR"=c(NA,HR),"CI5"=c(NA,CI5),"CI95"=c(NA,CI95),"se"=c(NA,se),"p-value"=c(NA,PValue),
                           "number of participant"=table(all_data$eGFR_scys_cat)[1:5],
                           "number of outcomes"=table(all_data$sepsis_event,all_data$eGFR_scys_cat)[2,1:5],
                           "time at risk"=get_time_risk("eGFR_scys_cat","sepsis_time"),
                           "incidence"=get_incidence("eGFR_scys_cat","sepsis_event","sepsis_time"))


model2.4 <- as.formula(paste("Surv(skin_soft_tissue_time, skin_soft_tissue_event) ~  eGFR_scys_cat+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm2.4 <- coxph(as.formula(model2.4), data = all_data)

res=summary(sm2.4)
coef=round(res$coefficients[1:4,1],3)
HR=round(res$coefficients[1:4,2],3)
se=round(res$coefficients[1:4,3],3)
PValue=round(res$coefficients[1:4,5],3)
CI5=round(res$conf.int[1:4,3],3)
CI95=round(res$conf.int[1:4,4],3)
model2.4_result=data.frame("model"=c("model2.4","model2.4","model2.4","model2.4","model2.4"),
                           "exposure"=c("eGFR_scys(90-105)","eGFR_scys(<60)","eGFR_scys(60-75)","eGFR_scys(75-90)","eGFR_scys(>=105)"),
                           "outcome"=c("skin_soft_tissue","skin_soft_tissue","skin_soft_tissue","skin_soft_tissue","skin_soft_tissue"),
                           "coef"=c(NA,coef),"HR"=c(NA,HR),"CI5"=c(NA,CI5),"CI95"=c(NA,CI95),"se"=c(NA,se),"p-value"=c(NA,PValue),
                           "number of participant"=table(all_data$eGFR_scys_cat)[1:5],
                           "number of outcomes"=table(all_data$skin_soft_tissue_event,all_data$eGFR_scys_cat)[2,1:5],
                           "time at risk"=get_time_risk("eGFR_scys_cat","skin_soft_tissue_time"),
                           "incidence"=get_incidence("eGFR_scys_cat","skin_soft_tissue_event","skin_soft_tissue_time"))



model2.5 <- as.formula(paste("Surv(urinary_tract_time, urinary_tract_event) ~  eGFR_scys_cat+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm2.5 <- coxph(as.formula(model2.5), data = all_data)

res=summary(sm2.5)
coef=round(res$coefficients[1:4,1],3)
HR=round(res$coefficients[1:4,2],3)
se=round(res$coefficients[1:4,3],3)
PValue=round(res$coefficients[1:4,5],3)
CI5=round(res$conf.int[1:4,3],3)
CI95=round(res$conf.int[1:4,4],3)
model2.5_result=data.frame("model"=c("model2.5","model2.5","model2.5","model2.5","model2.5"),
                           "exposure"=c("eGFR_scys(90-105)","eGFR_scys(<60)","eGFR_scys(60-75)","eGFR_scys(75-90)","eGFR_scys(>=105)"),
                           "outcome"=c("urinary_tract","urinary_tract","urinary_tract","urinary_tract","urinary_tract"),
                           "coef"=c(NA,coef),"HR"=c(NA,HR),"CI5"=c(NA,CI5),"CI95"=c(NA,CI95),"se"=c(NA,se),"p-value"=c(NA,PValue),
                           "number of participant"=table(all_data$eGFR_scys_cat)[1:5],
                           "number of outcomes"=table(all_data$urinary_tract_event,all_data$eGFR_scys_cat)[2,1:5],
                           "time at risk"=get_time_risk("eGFR_scys_cat","urinary_tract_time"),
                           "incidence"=get_incidence("eGFR_scys_cat","urinary_tract_event","urinary_tract_time"))



#BUN
model3.1 <- as.formula(paste("Surv(all_infection_time, all_infection_event) ~  BUN_cat+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm3.1 <- coxph(as.formula(model3.1), data = all_data)


res=summary(sm3.1)
coef=round(res$coefficients[1:3,1],3)
HR=round(res$coefficients[1:3,2],3)
se=round(res$coefficients[1:3,3],3)
PValue=round(res$coefficients[1:3,5],3)
CI5=round(res$conf.int[1:3,3],3)
CI95=round(res$conf.int[1:3,4],3)
model3.1_result=data.frame("model"=c("model3.1","model3.1","model3.1","model3.1"),
                           "exposure"=c("BUN_IQR1(0.81-4.52)","BUN_IQR2(4.52-5.29)","BUN_IQR3(5.29-6.15)","BUN_IQR4(6.15-41.83)"),
                           "outcome"=c("all_infection","all_infection","all_infection","all_infection"),
                           "coef"=c(NA,coef),"HR"=c(NA,HR),"CI5"=c(NA,CI5),"CI95"=c(NA,CI95),"se"=c(NA,se),"p-value"=c(NA,PValue),
                           "number of participant"=table(all_data$BUN_cat)[1:4],
                           "number of outcomes"=table(all_data$all_infection_event,all_data$BUN_cat)[2,1:4],
                           "time at risk"=get_time_risk("BUN_cat","all_infection_time"),
                           "incidence"=get_incidence("BUN_cat","all_infection_event","all_infection_time"))



model3.2 <- as.formula(paste("Surv(pneumonia_time, pneumonia_event) ~  BUN_cat+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm3.2 <- coxph(as.formula(model3.2), data = all_data)

res=summary(sm3.2)
coef=round(res$coefficients[1:3,1],3)
HR=round(res$coefficients[1:3,2],3)
se=round(res$coefficients[1:3,3],3)
PValue=round(res$coefficients[1:3,5],3)
CI5=round(res$conf.int[1:3,3],3)
CI95=round(res$conf.int[1:3,4],3)
model3.2_result=data.frame("model"=c("model3.2","model3.2","model3.2","model3.2"),
                           "exposure"=c("BUN_IQR1(0.81-4.52)","BUN_IQR2(4.52-5.29)","BUN_IQR3(5.29-6.15)","BUN_IQR4(6.15-41.83)"),
                           "outcome"=c("Pneumonia","Pneumonia","Pneumonia","Pneumonia"),
                           "coef"=c(NA,coef),"HR"=c(NA,HR),"CI5"=c(NA,CI5),"CI95"=c(NA,CI95),"se"=c(NA,se),"p-value"=c(NA,PValue),
                           "number of participant"=table(all_data$BUN_cat)[1:4],
                           "number of outcomes"=table(all_data$pneumonia_event,all_data$BUN_cat)[2,1:4],
                           "time at risk"=get_time_risk("BUN_cat","pneumonia_time"),
                           "incidence"=get_incidence("BUN_cat","pneumonia_event","pneumonia_time"))


model3.3 <- as.formula(paste("Surv(sepsis_time, sepsis_event) ~  BUN_cat+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm3.3 <- coxph(as.formula(model3.3), data = all_data)

res=summary(sm3.3)
coef=round(res$coefficients[1:3,1],3)
HR=round(res$coefficients[1:3,2],3)
se=round(res$coefficients[1:3,3],3)
PValue=round(res$coefficients[1:3,5],3)
CI5=round(res$conf.int[1:3,3],3)
CI95=round(res$conf.int[1:3,4],3)
n=res$n
nevent=res$nevent
model3.3_result=data.frame("model"=c("model3.3","model3.3","model3.3","model3.3"),
                           "exposure"=c("BUN_IQR1(0.81-4.52)","BUN_IQR2(4.52-5.29)","BUN_IQR3(5.29-6.15)","BUN_IQR4(6.15-41.83)"),
                           "outcome"=c("Sepsis","Sepsis","Sepsis","Sepsis"),
                           "coef"=c(NA,coef),"HR"=c(NA,HR),"CI5"=c(NA,CI5),"CI95"=c(NA,CI95),"se"=c(NA,se),"p-value"=c(NA,PValue),
                           "number of participant"=table(all_data$BUN_cat)[1:4],
                           "number of outcomes"=table(all_data$sepsis_event,all_data$BUN_cat)[2,1:4],
                           "time at risk"=get_time_risk("BUN_cat","sepsis_time"),
                           "incidence"=get_incidence("BUN_cat","sepsis_event","sepsis_time"))


model3.4 <- as.formula(paste("Surv(skin_soft_tissue_time, skin_soft_tissue_event) ~  BUN_cat+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm3.4 <- coxph(as.formula(model3.4), data = all_data)

res=summary(sm3.4)
coef=round(res$coefficients[1:3,1],3)
HR=round(res$coefficients[1:3,2],3)
se=round(res$coefficients[1:3,3],3)
PValue=round(res$coefficients[1:3,5],3)
CI5=round(res$conf.int[1:3,3],3)
CI95=round(res$conf.int[1:3,4],3)
model3.4_result=data.frame("model"=c("model3.4","model3.4","model3.4","model3.4"),
                           "exposure"=c("BUN_IQR1(0.81-4.52)","BUN_IQR2(4.52-5.29)","BUN_IQR3(5.29-6.15)","BUN_IQR4(6.15-41.83)"),
                           "outcome"=c("skin_soft_tissue","skin_soft_tissue","skin_soft_tissue","skin_soft_tissue"),
                           "coef"=c(NA,coef),"HR"=c(NA,HR),"CI5"=c(NA,CI5),"CI95"=c(NA,CI95),"se"=c(NA,se),"p-value"=c(NA,PValue),
                           "number of participant"=table(all_data$BUN_cat)[1:4],
                           "number of outcomes"=table(all_data$skin_soft_tissue_event,all_data$BUN_cat)[2,1:4],
                           "time at risk"=get_time_risk("BUN_cat","skin_soft_tissue_time"),
                           "incidence"=get_incidence("BUN_cat","skin_soft_tissue_event","skin_soft_tissue_time"))


model3.5 <- as.formula(paste("Surv(urinary_tract_time, urinary_tract_event) ~  BUN_cat+
                             age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
                             cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i"))
sm3.5 <- coxph(as.formula(model3.5), data = all_data)

res=summary(sm3.5)
coef=round(res$coefficients[1:3,1],3)
HR=round(res$coefficients[1:3,2],3)
se=round(res$coefficients[1:3,3],3)
PValue=round(res$coefficients[1:3,5],3)
CI5=round(res$conf.int[1:3,3],3)
CI95=round(res$conf.int[1:3,4],3)
model3.5_result=data.frame("model"=c("model3.5","model3.5","model3.5","model3.5"),
                           "exposure"=c("BUN_IQR1(0.81-4.52)","BUN_IQR2(4.52-5.29)","BUN_IQR3(5.29-6.15)","BUN_IQR4(6.15-41.83)"),
                           "outcome"=c("urinary_tract","urinary_tract","urinary_tract","urinary_tract"),
                           "coef"=c(NA,coef),"HR"=c(NA,HR),"CI5"=c(NA,CI5),"CI95"=c(NA,CI95),"se"=c(NA,se),"p-value"=c(NA,PValue),
                           "number of participant"=table(all_data$BUN_cat)[1:4],
                           "number of outcomes"=table(all_data$urinary_tract_event,all_data$BUN_cat)[2,1:4],
                           "time at risk"=get_time_risk("BUN_cat","urinary_tract_time"),
                           "incidence"=get_incidence("BUN_cat","urinary_tract_event","urinary_tract_time"))

model_result_cat=rbind(model1.1_result,model1.2_result,model1.3_result,model1.4_result,model1.5_result,model2.1_result,
                   model2.2_result,model2.3_result,model2.4_result,model2.5_result,model3.1_result,
                   model3.2_result,model3.3_result,model3.4_result,model3.5_result)



write.csv(model_result_cat,"cohort_study_category_new.csv")






