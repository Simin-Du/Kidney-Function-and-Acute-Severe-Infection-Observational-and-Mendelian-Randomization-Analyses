library(Matching)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
all_data=read.table("all_Data.txt",header=T)

####Inclusion
#critical 1
exclude_index1=which(is.na(all_data$log_eGFR_scys)&is.na(all_data$log_eGFR_scr)&is.na(all_data$BUN))
length(exclude_index1)#32586
all_data=all_data[-exclude_index1,]#469824

#critical 3
length(which(all_data$ethnicity!="white"))
length(which(is.na(all_data$ethnicity)))
all_data=all_data[which(all_data$ethnicity=="white"),]
dim(all_data)[1]#442817

####Exclusion
#critical 1 2
snp_exist=read.table("eid_wGenotype.txt")
all_data=merge(all_data,snp_exist,by.x = "eid",by.y = "V1")
length(all_data$eid)#436668

#critical 3
exclude_index3=unique(which(is.na(all_data$sex)),which(is.na(all_data$age)))#0

#critical 4
#kidney_beforebaseline=c()
#for (i in 1:length(all_data$eid)){
#  if(is.na(all_data$kidney_disease[i])){
#    kidney_beforebaseline=c(kidney_beforebaseline,0)
#  }else if(all_data$kidney_disease[i]<=all_data$recruitment[i]){
#    kidney_beforebaseline=c(kidney_beforebaseline,1)
#  }else{
#    kidney_beforebaseline=c(kidney_beforebaseline,0)
#  }
#}
all_data$kidney_beforebaseline=kidney_beforebaseline
exclude_index4=which(all_data$kidney_beforebaseline==1)
sum(kidney_beforebaseline)#476
all_data=all_data[-exclude_index4,]

#critical 5
#exclude_index5=which(all_data$renal_90_event==1)
#length(exclude_index5)#85
#all_data=all_data[-exclude_index5,]

#critical 6
exclude_index6=which(all_data$check_pregnancy==1)#116
length(exclude_index6)
all_data=all_data[-exclude_index6,]

#critical 7
#exclude_index7=unique(which(all_data$all_infection_beforebaseline==1))#12950
#length(exclude_index7)
#all_data=all_data[-exclude_index7,]

#critical 8
exclude_index8=unique(c(which(all_data$all_infection_lost_time==1),which(all_data$urinary_tract_lost_time==1),which(all_data$pneumonia_lost_time==1),
                      which(all_data$sepsis_lost_time==1),which(all_data$skin_soft_tissue_lost_time==1)))
all_data=all_data[-exclude_index8,]

dim(all_data)[1]#423039 #436074

sum(!is.na(all_data$log_eGFR_scr))
sum(!is.na(all_data$log_eGFR_scys))
sum(!is.na(all_data$BUN))

#matching
all_data$sex=factor(all_data$sex,ordered = F)
fit1<-glm(all_infection_event~sex+age,data=all_data,family=binomial(link="logit"))
X=fit1$fitted
Tr=all_data$all_infection_event
Match(X=X, Tr=Tr,M=4)

write.table(all_data,"all_data_exclusion_new.txt")
write.table(all_data,"all_data_exclusion.txt")



#############################################
#Sensitive Analysis for MR
all_data=read.table("all_Data.txt",header=T)

####Inclusion
#critical 1
exclude_index1=which(is.na(all_data$log_eGFR_scys)&is.na(all_data$log_eGFR_scr)&is.na(all_data$BUN))
length(exclude_index1)#32586
all_data=all_data[-exclude_index1,]#469824

#critical 3
length(which(all_data$ethnicity!="white"))
length(which(is.na(all_data$ethnicity)))
all_data=all_data[which(all_data$ethnicity=="white"),]
dim(all_data)[1]#442817

####Exclusion
#critical 1 2
snp_exist=read.table("eid_wGenotype.txt")
all_data=merge(all_data,snp_exist,by.x = "eid",by.y = "V1")
length(all_data$eid)#436668

#critical 3
exclude_index3=unique(which(is.na(all_data$sex)),which(is.na(all_data$age)))#0

#critical 4
#kidney_beforebaseline=c()
#for (i in 1:length(all_data$eid)){
#  if(is.na(all_data$kidney_disease[i])){
#    kidney_beforebaseline=c(kidney_beforebaseline,0)
#  }else if(all_data$kidney_disease[i]<=all_data$recruitment[i]){
#    kidney_beforebaseline=c(kidney_beforebaseline,1)
#  }else{
#    kidney_beforebaseline=c(kidney_beforebaseline,0)
#  }
#}
all_data$kidney_beforebaseline=kidney_beforebaseline
all_data$log_eGFR_scr[which(all_data$kidney_beforebaseline==1)]=log(3)
all_data$log_eGFR_scys[which(all_data$kidney_beforebaseline==1)]=log(3)

#critical 5
#exclude_index5=which(all_data$renal_90_event==1)
#length(exclude_index5)#85
#all_data=all_data[-exclude_index5,]

#critical 6
exclude_index6=which(all_data$check_pregnancy==1)#116
length(exclude_index6)
all_data=all_data[-exclude_index6,]

#critical 7
#exclude_index7=unique(which(all_data$all_infection_beforebaseline==1))#12950
#length(exclude_index7)
#all_data=all_data[-exclude_index7,]

#critical 8
exclude_index8=unique(c(which(all_data$all_infection_lost_time==1),which(all_data$urinary_tract_lost_time==1),which(all_data$pneumonia_lost_time==1),
                        which(all_data$sepsis_lost_time==1),which(all_data$skin_soft_tissue_lost_time==1)))
all_data=all_data[-exclude_index8,]

dim(all_data)[1]#436550

sum(!is.na(all_data$log_eGFR_scr))#436038
sum(!is.na(all_data$log_eGFR_scys))#436221
sum(!is.na(all_data$BUN))#435957


write.table(all_data,"all_data_exclusion_sensitive.txt")

sensitive_data=read.table("all_data_exclusion_sensitive.txt")



