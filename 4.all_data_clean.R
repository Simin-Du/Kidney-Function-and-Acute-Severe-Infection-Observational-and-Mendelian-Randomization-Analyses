setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
all_data_source=read.table("eGFR_ASCVD_cohort_data_clean.txt",header = T)
all_data1=read.table("infection_FID-0.0_50388.txt",header = T)
all_data2=read.table("infection_FID-0.0_51001.txt",header = T)
all_data3=read.table("infection_data.txt",header=T)
all_data4=read.table("merge_result_f.txt",header=T)
all_data5=read.table("all_data_ESRD.txt",header=T)
death_ESRD=read.table("death_ESRD.txt",header=T)
all_data_30880=read.table("MAFLD_FID-0.0_50388.txt",header = T)
medication_data=read.table("medication_data_clean_V3.txt",header=T)
colnames(all_data1)=c("FID","BUN","low_density_lipoprotein_cl","high_density_lipoprotein_cl","vitamin_D","attending_assessment_centre_date","treatment_code",
                      "lost_to_follow_up","kidney_disease","renal_disease","heterozygosity_missing_rate","death_date","cause_of_death")
#colnames(all_data2)=c("outliers_heterozygosity","recommended_genomic_analysis_exclusions","genetic_relatedness_exclusions")

#####Chronic pulmonary disease
all_data_chronic=merge(all_data2[,c("eid","X131486.0.0","X131487.0.0","X131488.0.0","X131489.0.0","X131492.0.0","X131493.0.0")],death_ESRD,by.x = "eid",by.y = "eid")
chronic_pulmonary_disease=c()
for(i in 1:length(all_data_chronic$eid)){
  has_disease=0
  if(!is.na(all_data_chronic$X131486.0.0[i]) & all_data_chronic$X131486.0.0[i]<=all_data_chronic$recruitment[i]){
    has_disease=1
  }else if(!is.na(all_data_chronic$X131488.0.0[i])& all_data_chronic$X131488.0.0[i]<=all_data_chronic$recruitment[i]){
    has_disease=1
  }else if(!is.na(all_data_chronic$X131492.0.0[i])& all_data_chronic$X131492.0.0[i]<=all_data_chronic$recruitment[i]){
    has_disease=1
  }
  chronic_pulmonary_disease=c(chronic_pulmonary_disease,has_disease)
}
sum(chronic_pulmonary_disease)
all_data_chronic$chronic_pulmonary_disease=chronic_pulmonary_disease
write.table(all_data_chronic,"all_data_chronic.txt")
all_data_chronic=read.table("all_data_chronic.txt",header=T)
#############################3
all_data_30880=all_data_30880[,c("eid","X30880.0.0","X30710.0.0")]
colnames(all_data_30880)=c("eid","unic_acid","c_reactive_protein")
all_data=merge(all_data_source[,c("eid","log_eGFR_scr","log_eGFR_scys","sex","age","ethnicity",
                                  "education","TDI","smoking","alcohol","MET","CVD","cancer",
                                  "diabetes","triglycerides","glucose","HbA1c","UACR_cat","BMI","WC",
                                  "SBP","DBP","mean_grip_strength","hypertension")],all_data_30880,by.x = "eid",by.y = "eid")
all_data=merge(all_data,all_data1,by.x = "eid",by.y = "FID")
all_data2_merge=all_data2[,c("eid","X22027.0.0","X22010.0.0","X22018.0.0")]
colnames(all_data2_merge)=c("eid","outliers_heterozygosity","recommended_genomic_analysis_exclusions","genetic_relatedness_exclusions")
all_data=merge(all_data,all_data2_merge,by.x = "eid",by.y = "eid")
#all_data_chronic
all_data=merge(all_data,all_data_chronic[,c("eid","chronic_pulmonary_disease")],by.x="eid",by.y = "eid")
all_data=merge(all_data,all_data3,by.x = "eid",by.y="eid")
all_data=merge(all_data,medication_data[,c("eid","oral_hypoglycemic_drugs","insulin")],by.x = "eid",by.y = "eid")
all_data=merge(all_data,all_data4[,c("eid","check_pregnancy")],by.x = "eid",by.y = "eid")
all_data5$WHR=all_data5$WC/all_data5$HC
all_data=merge(all_data,all_data5[,c("eid","WHR")],by.x="eid",by.y = "eid")

all_data=read.table("all_Data.txt",header=T)
dim(all_data)[1]
all_data3=read.table("all_data3.txt",header=T)
all_data3$renal_90_event=all_data3$renal_90_even
all_data=merge(all_data,all_data3[,c("eid","renal_90_event")],by.x="eid",by.y="eid")
write.table(all_data,"all_Data.txt")
colnames(all_data)


########Infection
CDK_FID_extract=read.table("CDK_FID_ICD9_10_extract.txt",head=T)
newdata=merge(CDK_FID_extract, death_ESRD, by.x = "eid", by.y = "eid")
all_data3=newdata[,c("eid","ESRD_date","death_date","lost_to_follow_up","recruitment")]
#all_data3=read.table("infection_data.txt",header=T)
uni_all_infections=c("A","B")
uni_pneumonia=c("J100", "J110", "J111", "J120", "J121", "J122", "J123", "J128", "J129", 
                "J13", "J14", "J150", "J151", "J152", "J153", "J154", "J155", "J156", 
                "J157", "J158", "J159", "J160", "J168", "J170", "J172", "J173", "J178", 
                "J180", "J181", "J188", "J189")
uni_sepsis=c("A021", "A227", "A327", "A427", "A400", "A401", "A402", "A403", "A408", 
             "A409", "A410", "A411", "A412", "A413", "A414", "A415", "A418", "A419", 
             "B377", "O85", "R572", "R650", "R65.1")
uni_skin_soft_tissue=c("L00", "L010", "L011", "L020", "L021", "L022", "L023", "L024", 
                       "L028", "L029", "L030", "L031", "L032", "L033", "L038", "L039", 
                       "L040", "L041", "L042", "L043", "L048", "L049", "L080", "L081", 
                       "L088", "L089")
uni_urinary_tract=c("N10", "N390")


#colname_ICD9_disease=paste("X41271.0.",0:46,sep="")
#colname_ICD9_date=paste("X41281.0.",0:46,sep="")
colname_disease=paste("X41270.0.",0:242,sep="")
colname_date=paste("X41280.0.",0:242,sep="")

all_infections_event=c()
all_infections_time=c()
all_infections_beforebaseline=c()
all_infections_date=c()
all_infections_lost_time=c()

for(j in c(1:dim(newdata)[1])){
  time="9999-12-31"
  for(i in 1:length(uni_all_infections)){
    sign=paste("^",uni_all_infections[i],sep="")#regular expression -- to make sure that it is the prefix of following strings
    is_in=grep(sign,newdata[j,colname_disease])
    if(any(is_in)){
      time1=min(c(t(newdata[j,colname_date[is_in]])))
      time=min(c(time,time1))
    }
  }
  if(is.na(time)){
    all_infections_event=c(all_infections_event,1)
    #all_infections_time=c(all_infections_time,NA)
    all_infections_beforebaseline=c(all_infections_beforebaseline,0)
    #all_infections_date=c(all_infections_date,NA)
    all_infections_lost_time=c(all_infections_lost_time,1)
    if(is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      all_infections_date=c(all_infections_date,"2021-11-12")
    }else if(is.na(newdata[j,"death_date"])&!is.na(newdata[j,"lost_to_follow_up"])){
      all_infections_date=c(all_infections_date,newdata[j,"lost_to_follow_up"])
    }else if(!is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      all_infections_date=c(all_infections_date,newdata[j,"death_date"])
    }else{
      all_infections_date=c(all_infections_date,newdata[j,"death_date"])
    }
    all_infections_time=c(all_infections_time,difftime(all_infections_date[j],newdata$recruitment[j]))
  }else if(time=="9999-12-31"){
    all_infections_event=c(all_infections_event,0)
    all_infections_beforebaseline=c(all_infections_beforebaseline,0)
    all_infections_lost_time=c(all_infections_lost_time,0)
    if(is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      all_infections_date=c(all_infections_date,"2021-11-12")
    }else if(is.na(newdata[j,"death_date"])&!is.na(newdata[j,"lost_to_follow_up"])){
      all_infections_date=c(all_infections_date,newdata[j,"lost_to_follow_up"])
    }else if(!is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      all_infections_date=c(all_infections_date,newdata[j,"death_date"])
    }else{
      all_infections_date=c(all_infections_date,newdata[j,"death_date"])
    }
    all_infections_time=c(all_infections_time,difftime(all_infections_date[j],newdata$recruitment[j]))
  }else{
    all_infections_event=c(all_infections_event,1)
    all_infections_time=c(all_infections_time,difftime(time,newdata$recruitment[j]))
    all_infections_date=c(all_infections_date,time)
    all_infections_lost_time=c(all_infections_lost_time,0)
    if(time<=newdata[j,"recruitment"]){
      all_infections_beforebaseline=c(all_infections_beforebaseline,1)
    }else{
      all_infections_beforebaseline=c(all_infections_beforebaseline,0)
    }
  }
}

all_data3$all_infection_event=all_infections_event
all_data3$all_infection_time=all_infections_time
all_data3$all_infection_beforebaseline=all_infections_beforebaseline
all_data3$all_infection_date=all_infections_date
all_data3$all_infection_lost_time=all_infections_lost_time

###########################

pneumonia_event=c()
pneumonia_time=c()
pneumonia_beforebaseline=c()
pneumonia_date=c()
pneumonia_lost_time=c()

for(j in c(1:dim(newdata)[1])){
  time="9999-12-31"
  for(i in 1:length(uni_pneumonia)){
    sign=paste("^",uni_pneumonia[i],sep="")#regular expression -- to make sure that it is the prefix of following strings
    is_in=grep(sign,newdata[j,colname_disease])
    if(any(is_in)){
      time1=min(c(t(newdata[j,colname_date[is_in]])))
      time=min(c(time,time1))
    }
  }
  if(is.na(time)){
    pneumonia_event=c(pneumonia_event,1)
    #pneumonia_time=c(pneumonia_time,NA)
    pneumonia_beforebaseline=c(pneumonia_beforebaseline,0)
    #pneumonia_date=c(pneumonia_date,NA)
    pneumonia_lost_time=c(pneumonia_lost_time,1)
    if(is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      pneumonia_date=c(pneumonia_date,"2021-11-12")
    }else if(is.na(newdata[j,"death_date"])&!is.na(newdata[j,"lost_to_follow_up"])){
      pneumonia_date=c(pneumonia_date,newdata[j,"lost_to_follow_up"])
    }else if(!is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      pneumonia_date=c(pneumonia_date,newdata[j,"death_date"])
    }else{
      pneumonia_date=c(pneumonia_date,newdata[j,"death_date"])
    }
    pneumonia_time=c(pneumonia_time,difftime(pneumonia_date[j],newdata$recruitment[j]))
  }else if(time=="9999-12-31"){
    pneumonia_event=c(pneumonia_event,0)
    #pneumonia_time=c(pneumonia_time,NA)
    pneumonia_beforebaseline=c(pneumonia_beforebaseline,0)
    #pneumonia_date=c(pneumonia_date,NA)
    pneumonia_lost_time=c(pneumonia_lost_time,0)
    if(is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      pneumonia_date=c(pneumonia_date,"2021-11-12")
    }else if(is.na(newdata[j,"death_date"])&!is.na(newdata[j,"lost_to_follow_up"])){
      pneumonia_date=c(pneumonia_date,newdata[j,"lost_to_follow_up"])
    }else if(!is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      pneumonia_date=c(pneumonia_date,newdata[j,"death_date"])
    }else{
      pneumonia_date=c(pneumonia_date,newdata[j,"death_date"])
    }
    pneumonia_time=c(pneumonia_time,difftime(pneumonia_date[j],newdata$recruitment[j]))
  }else{
    pneumonia_event=c(pneumonia_event,1)
    pneumonia_time=c(pneumonia_time,difftime(time,newdata$recruitment[j]))
    pneumonia_date=c(pneumonia_date,time)
    pneumonia_lost_time=c(pneumonia_lost_time,0)
    if(time<=newdata[j,"recruitment"]){
      pneumonia_beforebaseline=c(pneumonia_beforebaseline,1)
    }else{
      pneumonia_beforebaseline=c(pneumonia_beforebaseline,0)
    }
  }
}

all_data3$pneumonia_event=pneumonia_event
all_data3$pneumonia_time=pneumonia_time
all_data3$pneumonia_beforebaseline=pneumonia_beforebaseline
all_data3$pneumonia_date=pneumonia_date
all_data3$pneumonia_lost_time=pneumonia_lost_time


#############################################


sepsis_event=c()
sepsis_time=c()
sepsis_beforebaseline=c()
sepsis_date=c()
sepsis_lost_time=c()

for(j in c(1:dim(newdata)[1])){
  time="9999-12-31"
  for(i in 1:length(uni_sepsis)){
    sign=paste("^",uni_sepsis[i],sep="")#regular expression -- to make sure that it is the prefix of following strings
    is_in=grep(sign,newdata[j,colname_disease])
    if(any(is_in)){
      time1=min(c(t(newdata[j,colname_date[is_in]])))
      time=min(c(time,time1))
    }
  }
  if(is.na(time)){
    sepsis_event=c(sepsis_event,1)
    #sepsis_time=c(sepsis_time,NA)
    sepsis_beforebaseline=c(sepsis_beforebaseline,0)
    #sepsis_date=c(sepsis_date,NA)
    sepsis_lost_time=c(sepsis_lost_time,1)
    if(is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      sepsis_date=c(sepsis_date,"2021-11-12")
    }else if(is.na(newdata[j,"death_date"])&!is.na(newdata[j,"lost_to_follow_up"])){
      sepsis_date=c(sepsis_date,newdata[j,"lost_to_follow_up"])
    }else if(!is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      sepsis_date=c(sepsis_date,newdata[j,"death_date"])
    }else{
      sepsis_date=c(sepsis_date,newdata[j,"death_date"])
    }
    sepsis_time=c(sepsis_time,difftime(sepsis_date[j],newdata$recruitment[j]))
  }else if(time=="9999-12-31"){
    sepsis_event=c(sepsis_event,0)
    #sepsis_time=c(sepsis_time,NA)
    sepsis_beforebaseline=c(sepsis_beforebaseline,0)
    #sepsis_date=c(sepsis_date,NA)
    sepsis_lost_time=c(sepsis_lost_time,0)
    if(is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      sepsis_date=c(sepsis_date,"2021-11-12")
    }else if(is.na(newdata[j,"death_date"])&!is.na(newdata[j,"lost_to_follow_up"])){
      sepsis_date=c(sepsis_date,newdata[j,"lost_to_follow_up"])
    }else if(!is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      sepsis_date=c(sepsis_date,newdata[j,"death_date"])
    }else{
      sepsis_date=c(sepsis_date,newdata[j,"death_date"])
    }
    sepsis_time=c(sepsis_time,difftime(sepsis_date[j],newdata$recruitment[j]))
  }else{
    sepsis_event=c(sepsis_event,1)
    sepsis_time=c(sepsis_time,difftime(time,newdata$recruitment[j]))
    sepsis_date=c(sepsis_date,time)
    sepsis_lost_time=c(sepsis_lost_time,0)
    if(time<=newdata[j,"recruitment"]){
      sepsis_beforebaseline=c(sepsis_beforebaseline,1)
    }else{
      sepsis_beforebaseline=c(sepsis_beforebaseline,0)
    }
  }
}

all_data3$sepsis_event=sepsis_event
all_data3$sepsis_time=sepsis_time
all_data3$sepsis_beforebaseline=sepsis_beforebaseline
all_data3$sepsis_date=sepsis_date
all_data3$sepsis_lost_time=sepsis_lost_time
#####################################################

skin_soft_tissue_event=c()
skin_soft_tissue_time=c()
skin_soft_tissue_beforebaseline=c()
skin_soft_tissue_date=c()
skin_soft_tissue_lost_time=c()

for(j in c(1:dim(newdata)[1])){
  time="9999-12-31"
  for(i in 1:length(uni_skin_soft_tissue)){
    sign=paste("^",uni_skin_soft_tissue[i],sep="")#regular expression -- to make sure that it is the prefix of following strings
    is_in=grep(sign,newdata[j,colname_disease])
    if(any(is_in)){
      time1=min(c(t(newdata[j,colname_date[is_in]])))
      time=min(c(time,time1))
    }
  }
  if(is.na(time)){
    skin_soft_tissue_event=c(skin_soft_tissue_event,1)
    #skin_soft_tissue_time=c(skin_soft_tissue_time,NA)
    skin_soft_tissue_beforebaseline=c(skin_soft_tissue_beforebaseline,0)
    #skin_soft_tissue_date=c(skin_soft_tissue_date,NA)
    skin_soft_tissue_lost_time=c(skin_soft_tissue_lost_time,1)
    if(is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      skin_soft_tissue_date=c(skin_soft_tissue_date,"2021-11-12")
    }else if(is.na(newdata[j,"death_date"])&!is.na(newdata[j,"lost_to_follow_up"])){
      skin_soft_tissue_date=c(skin_soft_tissue_date,newdata[j,"lost_to_follow_up"])
    }else if(!is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      skin_soft_tissue_date=c(skin_soft_tissue_date,newdata[j,"death_date"])
    }else{
      skin_soft_tissue_date=c(skin_soft_tissue_date,newdata[j,"death_date"])
    }
    skin_soft_tissue_time=c(skin_soft_tissue_time,difftime(skin_soft_tissue_date[j],newdata$recruitment[j]))
  }else if(time=="9999-12-31"){
    skin_soft_tissue_event=c(skin_soft_tissue_event,0)
    #skin_soft_tissue_time=c(skin_soft_tissue_time,NA)
    skin_soft_tissue_beforebaseline=c(skin_soft_tissue_beforebaseline,0)
    #skin_soft_tissue_date=c(skin_soft_tissue_date,NA)
    skin_soft_tissue_lost_time=c(skin_soft_tissue_lost_time,0)
    if(is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      skin_soft_tissue_date=c(skin_soft_tissue_date,"2021-11-12")
    }else if(is.na(newdata[j,"death_date"])&!is.na(newdata[j,"lost_to_follow_up"])){
      skin_soft_tissue_date=c(skin_soft_tissue_date,newdata[j,"lost_to_follow_up"])
    }else if(!is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      skin_soft_tissue_date=c(skin_soft_tissue_date,newdata[j,"death_date"])
    }else{
      skin_soft_tissue_date=c(skin_soft_tissue_date,newdata[j,"death_date"])
    }
    skin_soft_tissue_time=c(skin_soft_tissue_time,difftime(skin_soft_tissue_date[j],newdata$recruitment[j]))
  }else{
    skin_soft_tissue_event=c(skin_soft_tissue_event,1)
    skin_soft_tissue_time=c(skin_soft_tissue_time,difftime(time,newdata$recruitment[j]))
    skin_soft_tissue_date=c(skin_soft_tissue_date,time)
    skin_soft_tissue_lost_time=c(skin_soft_tissue_lost_time,0)
    if(time<=newdata[j,"recruitment"]){
      skin_soft_tissue_beforebaseline=c(skin_soft_tissue_beforebaseline,1)
    }else{
      skin_soft_tissue_beforebaseline=c(skin_soft_tissue_beforebaseline,0)
    }
  }
}

all_data3$skin_soft_tissue_event=skin_soft_tissue_event
all_data3$skin_soft_tissue_time=skin_soft_tissue_time
all_data3$skin_soft_tissue_beforebaseline=skin_soft_tissue_beforebaseline
all_data3$skin_soft_tissue_date=skin_soft_tissue_date
all_data3$skin_soft_tissue_lost_time=skin_soft_tissue_lost_time

############################################################33



urinary_tract_event=c()
urinary_tract_time=c()
urinary_tract_beforebaseline=c()
urinary_tract_date=c()
urinary_tract_lost_time=c()

for(j in c(1:dim(newdata)[1])){
  time="9999-12-31"
  for(i in 1:length(uni_urinary_tract)){
    sign=paste("^",uni_urinary_tract[i],sep="")#regular expression -- to make sure that it is the prefix of following strings
    is_in=grep(sign,newdata[j,colname_disease])
    if(any(is_in)){
      time1=min(c(t(newdata[j,colname_date[is_in]])))
      time=min(c(time,time1))
    }
  }
  if(is.na(time)){
    urinary_tract_event=c(urinary_tract_event,1)
    #urinary_tract_time=c(urinary_tract_time,NA)
    urinary_tract_beforebaseline=c(urinary_tract_beforebaseline,0)
    #urinary_tract_date=c(urinary_tract_date,NA)
    urinary_tract_lost_time=c(urinary_tract_lost_time,1)
    if(is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      urinary_tract_date=c(urinary_tract_date,"2021-11-12")
    }else if(is.na(newdata[j,"death_date"])&!is.na(newdata[j,"lost_to_follow_up"])){
      urinary_tract_date=c(urinary_tract_date,newdata[j,"lost_to_follow_up"])
    }else if(!is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      urinary_tract_date=c(urinary_tract_date,newdata[j,"death_date"])
    }else{
      urinary_tract_date=c(urinary_tract_date,newdata[j,"death_date"])
    }
    urinary_tract_time=c(urinary_tract_time,difftime(urinary_tract_date[j],newdata$recruitment[j]))
  }else if(time=="9999-12-31"){
    urinary_tract_event=c(urinary_tract_event,0)
    #urinary_tract_time=c(urinary_tract_time,NA)
    urinary_tract_beforebaseline=c(urinary_tract_beforebaseline,0)
    #urinary_tract_date=c(urinary_tract_date,NA)
    urinary_tract_lost_time=c(urinary_tract_lost_time,0)
    if(is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      urinary_tract_date=c(urinary_tract_date,"2021-11-12")
    }else if(is.na(newdata[j,"death_date"])&!is.na(newdata[j,"lost_to_follow_up"])){
      urinary_tract_date=c(urinary_tract_date,newdata[j,"lost_to_follow_up"])
    }else if(!is.na(newdata[j,"death_date"])&is.na(newdata[j,"lost_to_follow_up"])){
      urinary_tract_date=c(urinary_tract_date,newdata[j,"death_date"])
    }else{
      urinary_tract_date=c(urinary_tract_date,newdata[j,"death_date"])
    }
    urinary_tract_time=c(urinary_tract_time,difftime(urinary_tract_date[j],newdata$recruitment[j]))
  }else{
    urinary_tract_event=c(urinary_tract_event,1)
    urinary_tract_time=c(urinary_tract_time,difftime(time,newdata$recruitment[j]))
    urinary_tract_date=c(urinary_tract_date,time)
    urinary_tract_lost_time=c(urinary_tract_lost_time,0)
    if(time<=newdata[j,"recruitment"]){
      urinary_tract_beforebaseline=c(urinary_tract_beforebaseline,1)
    }else{
      urinary_tract_beforebaseline=c(urinary_tract_beforebaseline,0)
    }
  }
}

all_data3$urinary_tract_event=urinary_tract_event
all_data3$urinary_tract_time=urinary_tract_time
all_data3$urinary_tract_beforebaseline=urinary_tract_beforebaseline
all_data3$urinary_tract_date=urinary_tract_date
all_data3$urinary_tract_lost_time=urinary_tract_lost_time

write.table(all_data3,"infection_data.txt")

all_data=read.table("all_Data.txt",header = T)
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
write.table(all_data,"all_Data.txt")


####################################

CDK_FID_extract=read.table("CDK_FID_ICD9_10_extract.txt",head=T)
death_ESRD=read.table("death_ESRD.txt",header=T)
newdata=merge(CDK_FID_extract, death_ESRD, by.x = "eid", by.y = "eid")
#all_data3=newdata[,c("eid","ESRD_date","death_date","lost_to_follow_up","recruitment")]
all_data3=read.table("infection_data.txt",header=T)

uni_renal=c("N17","N19")
colname_disease=paste("X41270.0.",0:242,sep="")
colname_date=paste("X41280.0.",0:242,sep="")

renal_event=c()
renal_90_event=c()
renal_date=c()

for(j in c(1:dim(newdata)[1])){
  time="9999-12-31"
  for(i in 1:length(uni_renal)){
    sign=paste("^",uni_renal[i],sep="")#regular expression -- to make sure that it is the prefix of following strings
    is_in=grep(sign,newdata[j,colname_disease])
    if(any(is_in)){
      time1=min(c(t(newdata[j,colname_date[is_in]])))
      time=min(c(time,time1))
    }
  }
  if(is.na(time)){
    renal_event=c(renal_event,1)
    renal_90_event=c(renal_90_event,0)
    renal_date=c(renal_date,NA)
  }else if(time=="9999-12-31"){
    renal_event=c(renal_event,0)
    renal_90_event=c(renal_90_event,0)
    renal_date=c(renal_date,NA)
  }else{
    renal_event=c(renal_event,1)
    renal_date=c(renal_date,time)
    if(abs(difftime(time,newdata$recruitment[j]))<=90){
      renal_90_event=c(renal_90_event,1)
    }else{
      renal_90_event=c(renal_90_event,0)
    }
  }
}

all_data3$renal_90_event=renal_90_event
all_data3$renal_event=renal_event
all_data3$renal_date=renal_date

write.table(all_data3,"all_data3.txt")














































