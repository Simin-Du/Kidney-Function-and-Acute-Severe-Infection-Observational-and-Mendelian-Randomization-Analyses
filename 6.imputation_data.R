setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(mice)
#all_data=read.table("all_data_exclusion_sensitive.txt",header=T)
all_data=read.table("all_data_exclusion.txt",header=T)

#covariates: age, sex, education, TDI, smoking, alcohol, CVD, chronic_pulmonary_disease,
#cancer, diabetes, oral_hypoglycemic_drugs,insulin, glucose, HbA1c, BMI, WC, vitamin_D, 
#c_reactive_protein,mean_grip_strength, MET

#do not use characters data type for the categorical data:
all_data$education[which(all_data$education=="level1")] = 1
all_data$education[which(all_data$education=="level2")] = 2
all_data$education[which(all_data$education=="level3")] = 3
all_data$education[which(all_data$education=="level4")] = 4
#ordinal categorical type:
all_data$education = factor(all_data$education,ordered=T)

#smoking
all_data$smoking[which(all_data$smoking=="never")]=0
all_data$smoking[which(all_data$smoking=="previous")]=1
all_data$smoking[which(all_data$smoking=="current")]=2
all_data$smoking=factor(all_data$smoking,order=F)

#alcohol
all_data$alcohol[which(all_data$alcohol=="never")]=0
all_data$alcohol[which(all_data$alcohol=="previous")]=1
all_data$alcohol[which(all_data$alcohol=="current")]=2
all_data$alcohol=factor(all_data$alcohol,order=F)

#the columns need imputation
col_imp = c("education", "TDI", "smoking", "alcohol",
            "glucose", "HbA1c", "BMI", "WC", "vitamin_D", 
            "c_reactive_protein", "mean_grip_strength", "MET", "UACR_cat", "hypertension")

kknn_c=all_data
kknn_c_ori = kknn_c
#if only impute once
for(i in 1:1){
  imp = kknn_c[,col_imp]
  imp_res <- mice(imp,m=1,maxit = 10,seed = i+100) #do not define the method
  imp_res_c = complete(imp_res)
  kknn_c[,col_imp] = imp_res_c
  kknn_imp = kknn_c
  kknn_c = kknn_c_ori
  save(kknn_imp,file = paste("kknn_imp",i,".Rdata",sep=""))
}

kknn_imp$eGFR_scr=exp(kknn_imp$log_eGFR_scr)
kknn_imp$eGFR_scys=exp(kknn_imp$log_eGFR_scys)

all_data=kknn_imp

all_data$eGFR_scr=exp(all_data$log_eGFR_scr)
all_data$eGFR_scys=exp(all_data$log_eGFR_scys)

all_data[which(all_data$eGFR_scr<60),"eGFR_scr_cat"]=2
all_data[which(all_data$eGFR_scr>=60&all_data$eGFR_scr<75),"eGFR_scr_cat"]=3
all_data[which(all_data$eGFR_scr>=75&all_data$eGFR_scr<90),"eGFR_scr_cat"]=4
all_data[which(all_data$eGFR_scr>=90&all_data$eGFR_scr<105),"eGFR_scr_cat"]=1
all_data[which(all_data$eGFR_scr>=105),"eGFR_scr_cat"]=5

all_data[which(all_data$eGFR_scys<60),"eGFR_scys_cat"]=2
all_data[which(all_data$eGFR_scys>=60&all_data$eGFR_scys<75),"eGFR_scys_cat"]=3
all_data[which(all_data$eGFR_scys>=75&all_data$eGFR_scys<90),"eGFR_scys_cat"]=4
all_data[which(all_data$eGFR_scys>=90&all_data$eGFR_scys<105),"eGFR_scys_cat"]=1
all_data[which(all_data$eGFR_scys>=105),"eGFR_scys_cat"]=5

all_data[which(all_data$BUN<quantile(na.omit(all_data$BUN),probs = 1)),"BUN_cat"]=4
all_data[which(all_data$BUN<quantile(na.omit(all_data$BUN),probs = 0.75)),"BUN_cat"]=3
all_data[which(all_data$BUN<quantile(na.omit(all_data$BUN),probs = 0.5)),"BUN_cat"]=2
all_data[which(all_data$BUN<quantile(na.omit(all_data$BUN),probs = 0.25)),"BUN_cat"]=1

all_data$lost_to_follow_up=all_data$lost_to_follow_up.x
all_data$death_date=all_data$death_date.x


all_data=all_data[,-c(which(colnames(all_data)%in% c("lost_to_follow_up.x","lost_to_follow_up.y","death_date.x","death_date.y")))]
all_data_add1=read.csv("medication_data_clean.csv")
all_data1 <- merge(all_data, all_data_add1, by = "eid", all.x = TRUE)
all_data_add2=read.csv("data_need.csv")
all_data2 <- merge(all_data1, all_data_add2[,c("eid","eGFR_scr_scys")], by = "eid", all.x = TRUE)
write.table(all_data2,"all_data_imputation_cohort.txt")
#write.table(all_data,"all_data_imputation_new.txt")
#write.table(all_data,"all_data_imputation_sensitive.txt")

#############################################################################
library(tableone)
#all_data=read.table("all_data_imputation_cohort.txt",header=T)

#all_data1=read.table("all_data_exclusion_new.txt",header=T)
all_data1=read.table("all_data_imputation_new.txt",header = T)
all_data_add1=read.csv("medication_data_clean.csv")
all_data2 <- merge(all_data1, all_data_add1, by = "eid", all.x = TRUE)
all_data_add2=read.csv("data_need.csv")
all_data <- merge(all_data2, all_data_add2[,c("eid","eGFR_scr_scys")], by = "eid", all.x = TRUE)

all_data=all_data[which(!is.na(all_data$eGFR_scr)),]
all_data$eGFR_scr_cat=factor(all_data$eGFR_scr_cat,order=F)
all_data$eGFR_scys_cat=factor(all_data$eGFR_scys_cat,order=F)
all_data$BUN_cat=factor(all_data$BUN_cat,order=F)
all_data$education=factor(all_data$education)
all_data$sex=factor(all_data$sex)
all_data$UACR_cat=factor(all_data$UACR_cat)
all_data$smoking=factor(all_data$smoking)
all_data$alcohol=factor(all_data$alcohol)
all_data$CVD=factor(all_data$CVD)
all_data$cancer=factor(all_data$cancer)
all_data$diabetes=factor(all_data$diabetes)
all_data$all_infection_event=factor(all_data$all_infection_event)
all_data$sepsis_event=factor(all_data$sepsis_event)
all_data$pneumonia_event=factor(all_data$pneumonia_event)
all_data$skin_soft_tissue_event=factor(all_data$skin_soft_tissue_event)
all_data$urinary_tract_event=factor(all_data$urinary_tract_event)
all_data$hypertension=factor(all_data$hypertension)
all_data$medicine_g_i=factor(all_data$medicine_g_i,order=F)

var=c("sex","age","education","TDI","smoking","alcohol","MET","CVD",
      "cancer","diabetes","triglycerides","glucose","HbA1c","UACR_cat","BMI","WC","eGFR_scr","eGFR_scys",
      "SBP","DBP","mean_grip_strength","hypertension","unic_acid","c_reactive_protein","BUN","low_density_lipoprotein_cl",
      "high_density_lipoprotein_cl","vitamin_D","all_infection_event","sepsis_event",
      "pneumonia_event","skin_soft_tissue_event","urinary_tract_event","medicine_g_i","eGFR_scr_scys")
description=CreateTableOne(vars=var,strata="eGFR_scr_cat",all_data,addOverall=T,includeNA = TRUE)
tab1<-print(description, nonnormal = c("sex","age","education","TDI","smoking","alcohol","MET","CVD",
                                       "cancer","diabetes","triglycerides","glucose","HbA1c","UACR_cat","BMI","WC","eGFR_scr","eGFR_scys",
                                       "SBP","DBP","mean_grip_strength","hypertension","unic_acid","c_reactive_protein","BUN","low_density_lipoprotein_cl",
                                       "high_density_lipoprotein_cl","vitamin_D","all_infection_event","sepsis_event",
                                       "pneumonia_event","skin_soft_tissue_event","urinary_tract_event","medicine_g_i","eGFR_scr_scys"),formatOptions = list(big.mark = ","),quote = T)
write.csv(tab1,"Table1_MR_imputation.csv")
#write.csv(tab1,"Table1_before_imputation_new1.csv")
#write.csv(tab1,"Table1_cohort_imputation.csv")







####################################################################
all_data=read.table("all_data_imputation_new.txt",header = T)
all_data3=read.table("infection_data.txt",header=T)
time_data=merge(all_data[,c("eid","ESRD_date")],all_data3,by.x="eid",by.y="eid")
all_data$all_infection_time=time_data$all_infection_time
all_data$all_infection_date=time_data$all_infection_date
all_data$pneumonia_date=time_data$pneumonia_date
all_data$pneumonia_time=time_data$pneumonia_time
all_data$sepsis_date=time_data$sepsis_date
all_data$sepsis_time=time_data$sepsis_time
all_data$skin_soft_tissue_date=time_data$skin_soft_tissue_date
all_data$skin_soft_tissue_time=time_data$skin_soft_tissue_time
all_data$urinary_tract_date=time_data$urinary_tract_date
all_data$urinary_tract_time=time_data$urinary_tract_time

exclude_index8=unique(c(which(all_data$all_infection_lost_time==1),which(all_data$urinary_tract_lost_time==1),which(all_data$pneumonia_lost_time==1),
                      which(all_data$sepsis_lost_time==1),which(all_data$skin_soft_tissue_lost_time==1)))
#all_data=all_data[-exclude_index8,]

dim(all_data)[1]#423040 #436074

write.table(all_data,"all_data_imputation_new.txt")

all_data=read.table("all_data_imputation.txt",header=T)
colnames(all_data)


all_data=read.table("all_data_imputation.txt",header=T)
hist(all_data$eGFR_scr, breaks=100,xlim=c(0,150),col="lightsteelblue",main="Histogram of eGFRcr",xlab="eGFRcr")
max(all_data$eGFR_scr[which(!is.na(all_data$eGFR_scr))])
hist(all_data$eGFR_scys, breaks=100,xlim=c(0,150),col="peachpuff",main="Histogram of eGFRcys",xlab="eGFRcys")
hist(all_data$BUN, breaks=100,xlim=c(0,20),col="palegreen",main="Histogram of BUN",xlab="BUN")

