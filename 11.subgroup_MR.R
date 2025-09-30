library(dplyr)
library(survival)
library(cmprsk)
library(riskRegression)
library(car)
library(Publish)
library(Hmisc)
library(rms)
library(survival)
library(metafor); library(ggplot2)
library(ggpubr)
library("RColorBrewer")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
mycolors<-brewer.pal(9, "Blues")


source("nlme_summ_aes.R")


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#all_data
geno = read.table("31SNP_1Proxy.xmat",head=T)
all_data=read.table("all_data_imputation.txt",header = T)
PC_data=read.table("VTE_50388_PC1-20.txt",header=T)
PC_data=PC_data[,1:11]
colnames(PC_data)=c("eid","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")

geno <- geno[-1,]

for (i in 3:dim(geno)[2]) {
  geno[,i] <- as.numeric(as.character(geno[,i]))
}

all_data=merge(all_data,geno,by.x="eid",by.y="FID")
all_data=merge(all_data,PC_data[,1:11],by.x="eid",by.y = "eid")
snp_list=read.table("31SNP_1Proxy.txt")
all_data$eGFR_scr_10=all_data$eGFR_scr/10
all_data$eGFR_scys_10=all_data$eGFR_scys/10



#PRS
#eGFR_scr
forml = paste0("eGFR_scr_10 ~ ",paste0(snp_list$V1,collapse = "+"))
all_data$prs_eGFR_scr=predict(lm(as.formula(forml),data=all_data),newdata = all_data)
all_data_scr=all_data[which(!is.na(all_data$eGFR_scr_10 & all_data$prs_eGFR_scr)),]
#eGFR_scys
all_data_scys=all_data[which(!is.na(all_data$eGFR_scys_10 & all_data$prs_eGFR_scys)),]
forml = paste0("eGFR_scys_10 ~ ",paste0(snp_list$V1,collapse = "+"))
all_data$prs_eGFR_scys=predict(lm(as.formula(forml),data=all_data),newdata = all_data)
#BUN
all_data_BUN=all_data[which(!is.na(all_data$BUN & all_data$prs_BUN)),]
forml = paste0("BUN ~ ",paste0(snp_list$V1,collapse = "+"))
all_data$prs_BUN=predict(lm(as.formula(forml),data=all_data),newdata = all_data)




# #252snp prs
# score_scr=read.table("eGFR_scr_10.profile_score",header=T)
# score_scys=read.table("eGFR_cyst_10.profile_score",header=T)
# score_BUN=read.table("BUN.profile_score",header=T)
# 
# colnames(score_scr)=c("FID","prs_eGFR_scr")
# all_data_scr=merge(all_data,score_scr,by.x="eid",by.y = "FID")
# all_data_scr=all_data_scr[which(!is.na(all_data_scr$eGFR_scr_10)),]
# colnames(score_scys)=c("FID","prs_eGFR_scys")
# all_data_scys=merge(all_data,score_scys,by.x="eid",by.y = "FID")
# all_data_scys=all_data_scys[which(!is.na(all_data_scys$eGFR_scys_10)),]
# colnames(score_BUN)=c("FID","prs_BUN")
# all_data_BUN=merge(all_data,score_BUN,by.x="eid",by.y = "FID")
# all_data_BUN=all_data_BUN[which(!is.na(all_data_BUN$BUN)),]



#sTEP1&2
#eGFR_scr
all_data_scr$scr_hat=lm(eGFR_scr_10~prs_eGFR_scr,data=all_data_scr)$fit
all_data_scr$scr_residual=all_data_scr$eGFR_scr_10-all_data_scr$scr_hat
N=10
all_data_scr$group=cut(all_data_scr$scr_residual, labels = F,quantile(all_data_scr$scr_residual, prob = 0:N / N, names = FALSE), include = TRUE)
#eGFR_scys
all_data_scys$scys_hat=lm(eGFR_scys_10~prs_eGFR_scys,data=all_data_scys)$fit
all_data_scys$scys_residual=all_data_scys$eGFR_scys_10-all_data_scys$scys_hat
N=10
all_data_scys$group=cut(all_data_scys$scys_residual, labels = F,quantile(all_data_scys$scys_residual, prob = 0:N / N, names = FALSE), include = TRUE)
#BUN
all_data_BUN$BUN_hat=lm(BUN~prs_BUN,data=all_data_BUN)$fit
all_data_BUN$BUN_residual=all_data_BUN$BUN-all_data_BUN$BUN_hat
N=10
all_data_BUN$group=cut(all_data_BUN$BUN_residual, labels = F,quantile(all_data_BUN$BUN_residual, prob = 0:N / N, names = FALSE), include = TRUE)


##############################################subgroup3###############################################################
#subgroup3 male

#eGFR_scr
data=all_data_scr

data_eGFR_scr_male=data[which(data$sex=="male"),]
data_eGFR_scr_female=data[which(data$sex=="female"),]

#All_infection
#data_eGFR_scr_male
data_eGFR_scr_male$scr_hat=lm(eGFR_scr_10~prs_eGFR_scr,data=data_eGFR_scr_male)$fit
data_eGFR_scr_male$scr_residual=data_eGFR_scr_male$eGFR_scr_10-data_eGFR_scr_male$scr_hat
N=10
data_eGFR_scr_male$group=cut(data_eGFR_scr_male$scr_residual, labels = F,quantile(data_eGFR_scr_male$scr_residual, prob = 0:N / N, names = FALSE), include = TRUE)

data_eGFR_scr_male$follow_up_time = data_eGFR_scr_male$all_infection_time
data_eGFR_scr_male$event.type=data_eGFR_scr_male$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scr_male[which(data_eGFR_scr_male$group==i),]
  temp = summary(lm(eGFR_scr_10~prs_eGFR_scr+age+pc1+pc2,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scr_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scr+age+pc1+pc2,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
summary_stat_scr_allinfection=summary_stat
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeana02 = summary_stat$mean_expos # mean exposure in each quantile
xmeana02
xmean1 = xmeana02-6  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

aa12 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=6, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cr", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
aa12

aa12$coefficients
aa12$p_tests
aa12$p_heterogeneity

fp_aa12=ifelse(aa12$p_tests[1,2]<0.001," < 0.001",ifelse(aa12$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa12$p_tests[1,2],3))))
Q_aa12=ifelse(aa12$p_tests[1,4]<0.001," < 0.001",ifelse(aa12$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa12$p_tests[1,4],3))))

fa112=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scr_male$eGFR_scr, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scr_male$eGFR_scr), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),aa12$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), aa12$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),aa12$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeana02)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cr]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(130,50)) + ggtitle(expression(paste(eGFR["cr"]," and All Infections (Male)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_aa12) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_aa12)),  size = 3)

fa112
infa112=data.frame(subgroup="male",exposure="eGFR_cr",n=dim(data_eGFR_scr_male)[1],event_n=sum(data_eGFR_scr_male$event.type),
                  P_Cochrans_Q=aa12$p_tests[1,4],p_non_linear=aa12$p_tests[1,2])

#data_eGFR_scr_female
data_eGFR_scr_female$scr_hat=lm(eGFR_scr_10~prs_eGFR_scr,data=data_eGFR_scr_female)$fit
data_eGFR_scr_female$scr_residual=data_eGFR_scr_female$eGFR_scr_10-data_eGFR_scr_female$scr_hat
N=10
data_eGFR_scr_female$group=cut(data_eGFR_scr_female$scr_residual, labels = F,quantile(data_eGFR_scr_female$scr_residual, prob = 0:N / N, names = FALSE), include = TRUE)

data_eGFR_scr_female$follow_up_time = data_eGFR_scr_female$all_infection_time
data_eGFR_scr_female$event.type=data_eGFR_scr_female$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scr_female[which(data_eGFR_scr_female$group==i),]
  temp = summary(lm(eGFR_scr_10~prs_eGFR_scr+age+pc1+pc2,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scr_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scr+age+pc1+pc2,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeanb02 = summary_stat$mean_expos # mean exposure in each quantile
xmeanb02
xmean1 = xmeanb02-6  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

ba12 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=6, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cr", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
ba12

ba12$coefficients
ba12$p_tests
ba12$p_heterogeneity

fp_ba12=ifelse(ba12$p_tests[1,2]<0.001," < 0.001",ifelse(ba12$p_tests[1,2]>0.999," > 0.999",paste(" =",round(ba12$p_tests[1,2],3))))
Q_ba12=ifelse(ba12$p_tests[1,4]<0.001," < 0.001",ifelse(ba12$p_tests[1,4]>0.999," > 0.999",paste(" =",round(ba12$p_tests[1,4],3))))

fb112=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scr_female$eGFR_scr, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scr_female$eGFR_scr), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),ba12$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), ba12$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),ba12$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeanb02)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cr]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(130,50)) + ggtitle(expression(paste(eGFR["cr"]," and All Infections (Female)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_ba12) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_ba12)),  size = 3)

fb112
infb112=data.frame(subgroup="female",exposure="eGFR_cr",n=dim(data_eGFR_scr_female)[1],event_n=sum(data_eGFR_scr_female$event.type),
                  P_Cochrans_Q=ba12$p_tests[1,4],p_non_linear=ba12$p_tests[1,2])

#eGFR_scys

data=all_data_scys

data_eGFR_scys_male=data[which(data$sex=="male"),]
data_eGFR_scys_female=data[which(data$sex=="female"),]

#All_infection
#data_eGFR_scys_male
data_eGFR_scys_male$scys_hat=lm(eGFR_scys_10~prs_eGFR_scys,data=data_eGFR_scys_male)$fit
data_eGFR_scys_male$scys_residual=data_eGFR_scys_male$eGFR_scys_10-data_eGFR_scys_male$scys_hat
N=10
data_eGFR_scys_male$group=cut(data_eGFR_scys_male$scys_residual, labels = F,quantile(data_eGFR_scys_male$scys_residual, prob = 0:N / N, names = FALSE), include = TRUE)

data_eGFR_scys_male$follow_up_time = data_eGFR_scys_male$all_infection_time
data_eGFR_scys_male$event.type=data_eGFR_scys_male$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scys_male[which(data_eGFR_scys_male$group==i),]
  temp = summary(lm(eGFR_scys_10~prs_eGFR_scys+age+pc1+pc2,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scys_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scys+age+pc1+pc2,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeana22 = summary_stat$mean_expos # mean exposure in each quantile
xmeana22
xmean1 = xmeana22-5  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

aa22 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=5, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cys", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
aa22

aa22$coefficients
aa22$p_tests
aa22$p_heterogeneity

fp_aa22=ifelse(aa22$p_tests[1,2]<0.001," < 0.001",ifelse(aa22$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa22$p_tests[1,2],3))))
Q_aa22=ifelse(aa22$p_tests[1,4]<0.001," < 0.001",ifelse(aa22$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa22$p_tests[1,4],3))))

fa212=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scys_male$eGFR_scys, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scys_male$eGFR_scys), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),aa22$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), aa22$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),aa22$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeana22)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cys]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cys"]," and All Infections (Male)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_aa22) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_aa22)),  size = 3)

fa212
infa212=data.frame(subgroup="male",exposure="eGFR_cys",n=dim(data_eGFR_scys_male)[1],event_n=sum(data_eGFR_scys_male$event.type),
                  P_Cochrans_Q=aa22$p_tests[1,4],p_non_linear=aa22$p_tests[1,2])

#data_eGFR_scys_female
data_eGFR_scys_female$scys_hat=lm(eGFR_scys_10~prs_eGFR_scys,data=data_eGFR_scys_female)$fit
data_eGFR_scys_female$scys_residual=data_eGFR_scys_female$eGFR_scys_10-data_eGFR_scys_female$scys_hat
N=10
data_eGFR_scys_female$group=cut(data_eGFR_scys_female$scys_residual, labels = F,quantile(data_eGFR_scys_female$scys_residual, prob = 0:N / N, names = FALSE), include = TRUE)

data_eGFR_scys_female$follow_up_time = data_eGFR_scys_female$all_infection_time
data_eGFR_scys_female$event.type=data_eGFR_scys_female$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scys_female[which(data_eGFR_scys_female$group==i),]
  temp = summary(lm(eGFR_scys_10~prs_eGFR_scys+age+pc1+pc2,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scys_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scys+age+pc1+pc2,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeanb22 = summary_stat$mean_expos # mean exposure in each quantile
xmeanb22
xmean1 = xmeanb22-5  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

ba22 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=5, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cys", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
ba22

ba22$coefficients
ba22$p_tests
ba22$p_heterogeneity

fp_ba22=ifelse(ba22$p_tests[1,2]<0.001," < 0.001",ifelse(ba22$p_tests[1,2]>0.999," > 0.999",paste(" =",round(ba22$p_tests[1,2],3))))
Q_ba22=ifelse(ba22$p_tests[1,4]<0.001," < 0.001",ifelse(ba22$p_tests[1,4]>0.999," > 0.999",paste(" =",round(ba22$p_tests[1,4],3))))

fb212=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scys_female$eGFR_scys, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scys_female$eGFR_scys), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),ba22$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), ba22$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),ba22$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeanb22)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cys]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cys"]," and All Infections (Female)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_ba22) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_ba22)),  size = 3)

fb212
infb212=data.frame(subgroup="female",exposure="eGFR_cys",n=dim(data_eGFR_scys_female)[1],event_n=sum(data_eGFR_scys_female$event.type),
                  P_Cochrans_Q=ba22$p_tests[1,4],p_non_linear=ba22$p_tests[1,2])

#BUN
data=all_data_BUN
data_BUN_male=data[which(data$sex=="male"),]
data_BUN_female=data[which(data$sex=="female"),]

#All_infection
#data_BUN_male
data_BUN_male$BUN_hat=lm(BUN~prs_BUN,data=data_BUN_male)$fit
data_BUN_male$BUN_residual=data_BUN_male$BUN-data_BUN_male$BUN_hat
N=10
data_BUN_male$group=cut(data_BUN_male$BUN_residual, labels = F,quantile(data_BUN_male$BUN_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_BUN_male$follow_up_time = data_BUN_male$all_infection_time
data_BUN_male$event.type=data_BUN_male$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_BUN_male[which(data_BUN_male$group==i),]
  temp = summary(lm(BUN~prs_BUN+age+pc1+pc2,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$BUN)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_BUN+age+pc1+pc2,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeana32 = summary_stat$mean_expos # mean exposure in each quantile
xmeana32
xmean1 = xmeana32-3  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

aa32 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=3, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="BUN", pref_y="Hazard ratio of all infections", breaks=c(0.25,0.5,1,2,4,8))
aa32

aa32$coefficients
aa32$p_tests
aa32$p_heterogeneity

fp_aa32=ifelse(aa32$p_tests[1,2]<0.001," < 0.001",ifelse(aa32$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa32$p_tests[1,2],3))))
Q_aa32=ifelse(aa32$p_tests[1,4]<0.001," < 0.001",ifelse(aa32$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa32$p_tests[1,4],3))))

fa312=ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), data_BUN_male,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,12,1))+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x, y=log2(yest)),aa32$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x, y=log2(lci)), aa32$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x, y=log2(uci)),aa32$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeana32), y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 4),  labels =c(0.25,0.5,1,2,4,8,16),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~BUN))) + 
  scale_x_continuous(breaks=seq(0,12,by = 1)) + 
  coord_cartesian(xlim =c(0,12)) + ggtitle("BUN and All Infections (Male)")+
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
  annotate("text", x = 10, y = 3.7, label = bquote(P['Cochran\'s Q test'] ~ .(Q_aa32) ), size = 3) + 
  annotate("text", x = 10, y = 3.4, label = bquote(P['nonlinearity'] ~.(fp_aa32)),  size = 3)

fa312
infa312=data.frame(subgroup="male",exposure="BUN",n=dim(data_BUN_male)[1],event_n=sum(data_BUN_male$event.type),
                  P_Cochrans_Q=aa32$p_tests[1,4],p_non_linear=aa32$p_tests[1,2])

#data_BUN_female
data_BUN_female$BUN_hat=lm(BUN~prs_BUN,data=data_BUN_female)$fit
data_BUN_female$BUN_residual=data_BUN_female$BUN-data_BUN_female$BUN_hat
N=10
data_BUN_female$group=cut(data_BUN_female$BUN_residual, labels = F,quantile(data_BUN_female$BUN_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_BUN_female$follow_up_time = data_BUN_female$all_infection_time
data_BUN_female$event.type=data_BUN_female$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_BUN_female[which(data_BUN_female$group==i),]
  temp = summary(lm(BUN~prs_BUN+age+pc1+pc2,data=pheno))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$BUN)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_BUN+age+pc1+pc2,data = pheno))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeanb32 = summary_stat$mean_expos # mean exposure in each quantile
xmeanb32
xmean1 = xmeanb32-3  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

ba32 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=3, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="BUN", pref_y="Hazard ratio of all infections", breaks=c(0.25,0.5,1,2,4,8))
ba32

ba32$coefficients
ba32$p_tests
ba32$p_heterogeneity

fp_ba32=ifelse(ba32$p_tests[1,2]<0.001," < 0.001",ifelse(ba32$p_tests[1,2]>0.999," > 0.999",paste(" =",round(ba32$p_tests[1,2],3))))
Q_ba32=ifelse(ba32$p_tests[1,4]<0.001," < 0.001",ifelse(ba32$p_tests[1,4]>0.999," > 0.999",paste(" =",round(ba32$p_tests[1,4],3))))

fb312=ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), data_BUN_female,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,12,1))+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x, y=log2(yest)),ba32$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x, y=log2(lci)), ba32$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x, y=log2(uci)),ba32$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeanb32), y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 4),  labels =c(0.25,0.5,1,2,4,8,16),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~BUN))) + 
  scale_x_continuous(breaks=seq(0,12,by = 1)) + 
  coord_cartesian(xlim =c(0,12)) + ggtitle("BUN and All Infections (Female)")+
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
  annotate("text", x = 10, y = 3.7, label = bquote(P['Cochran\'s Q test'] ~ .(Q_ba32) ), size = 3) + 
  annotate("text", x = 10, y = 3.4, label = bquote(P['nonlinearity'] ~.(fp_ba32)),  size = 3)

fb312
infb312=data.frame(subgroup="female",exposure="BUN",n=dim(data_BUN_female)[1],event_n=sum(data_BUN_female$event.type),
                   P_Cochrans_Q=ba32$p_tests[1,4],p_non_linear=ba32$p_tests[1,2])

inf_sub1=rbind(infa112,infb112,infa212,infb212,infa312,infb312)
##############################################subgroup4#########################################################################
#subgroup4 diabete

#eGFR_scr
data=all_data_scr

#logistic_regression
fit = glm(diabetes~prs_eGFR_scr+BMI+age+sex+smoking+alcohol+MET+education+TDI+cancer,data = data,family = 'binomial')
data$prob = predict(fit,type="response")

#delete 1%
quantile = quantile(data$prob,probs=seq(0,1,0.01))
index = which(data$prob<quantile[2])
data_reduced = data[-index,]
data=data_reduced
data_eGFR_scr_diabetes=data[which(data$diabetes==1),]
data_eGFR_scr_nodiabetes=data[which(data$diabetes==0),]


#All_infection
#data_eGFR_scr_diabetes
data_eGFR_scr_diabetes$scr_hat=lm(eGFR_scr_10~prs_eGFR_scr,data=data_eGFR_scr_diabetes)$fit
data_eGFR_scr_diabetes$scr_residual=data_eGFR_scr_diabetes$eGFR_scr_10-data_eGFR_scr_diabetes$scr_hat
N=10
data_eGFR_scr_diabetes$group=cut(data_eGFR_scr_diabetes$scr_residual, labels = F,quantile(data_eGFR_scr_diabetes$scr_residual, prob = 0:N / N, names = FALSE), include = TRUE)

data_eGFR_scr_diabetes$follow_up_time = data_eGFR_scr_diabetes$all_infection_time
data_eGFR_scr_diabetes$event.type=data_eGFR_scr_diabetes$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scr_diabetes[which(data_eGFR_scr_diabetes$group==i),]
  temp = summary(lm(eGFR_scr_10~prs_eGFR_scr+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scr_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scr+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeana03 = summary_stat$mean_expos # mean exposure in each quantile
xmeana03
xmean1 = xmeana03-5  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

aa13 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d=1, offset=5, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cr", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
aa13

aa13$coefficients
aa13$p_tests
aa13$p_heterogeneity

fp_aa13=ifelse(aa13$p_tests[1,2]<0.001," < 0.001",ifelse(aa13$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa13$p_tests[1,2],3))))
Q_aa13=ifelse(aa13$p_tests[1,4]<0.001," < 0.001",ifelse(aa13$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa13$p_tests[1,4],3))))


fp_aa13=ifelse(aa13$p_tests[1,2]<0.001," < 0.001",ifelse(aa13$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa13$p_tests[1,2],3))))
Q_aa13=ifelse(aa13$p_tests[1,4]<0.001," < 0.001",ifelse(aa13$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa13$p_tests[1,4],3))))

fa113=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scr_diabetes$eGFR_scr, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scr_diabetes$eGFR_scr), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),aa13$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), aa13$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),aa13$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeana03)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cr]))) + 
  scale_x_continuous(breaks=seq(130,40,by = -10)) + 
  coord_cartesian(xlim =c(  130,40)) + ggtitle(expression(paste(eGFR["cr"]," and All Infections (With Diabetes)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_aa13) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_aa13)),  size = 3)

fa113
infa113=data.frame(subgroup="Diabetes",exposure="eGFR_scr",n=dim(data_eGFR_scr_diabetes)[1],event_n=sum(data_eGFR_scr_diabetes$event.type),
                   P_Cochrans_Q=aa13$p_tests[1,4],p_non_linear=aa13$p_tests[1,2])

#data_eGFR_scr_nodiabetes
data_eGFR_scr_nodiabetes$scr_hat=lm(eGFR_scr_10~prs_eGFR_scr,data=data_eGFR_scr_nodiabetes)$fit
data_eGFR_scr_nodiabetes$scr_residual=data_eGFR_scr_nodiabetes$eGFR_scr_10-data_eGFR_scr_nodiabetes$scr_hat
N=10
data_eGFR_scr_nodiabetes$group=cut(data_eGFR_scr_nodiabetes$scr_residual, labels = F,quantile(data_eGFR_scr_nodiabetes$scr_residual, prob = 0:N / N, names = FALSE), include = TRUE)

data_eGFR_scr_nodiabetes$follow_up_time = data_eGFR_scr_nodiabetes$all_infection_time
data_eGFR_scr_nodiabetes$event.type=data_eGFR_scr_nodiabetes$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scr_nodiabetes[which(data_eGFR_scr_nodiabetes$group==i),]
  temp = summary(lm(eGFR_scr_10~prs_eGFR_scr+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scr_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scr+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeanb03 = summary_stat$mean_expos # mean exposure in each quantile
xmeanb03
xmean1 = xmeanb03-6  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

ba13 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=6, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cr", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
ba13

ba13$coefficients
ba13$p_tests
ba13$p_heterogeneity

fp_ba13=ifelse(ba13$p_tests[1,2]<0.001," < 0.001",ifelse(ba13$p_tests[1,2]>0.999," > 0.999",paste(" =",round(ba13$p_tests[1,2],3))))
Q_ba13=ifelse(ba13$p_tests[1,4]<0.001," < 0.001",ifelse(ba13$p_tests[1,4]>0.999," > 0.999",paste(" =",round(ba13$p_tests[1,4],3))))

fb113=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scr_nodiabetes$eGFR_scr, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scr_nodiabetes$eGFR_scr), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),ba13$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), ba13$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),ba13$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeanb03)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cr]))) + 
  scale_x_continuous(breaks=seq(130,40,by = -10)) + 
  coord_cartesian(xlim =c(  130,40)) + ggtitle(expression(paste(eGFR["cr"]," and All Infections (Without Diabetes)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_ba13) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_ba13)),  size = 3)

fb113
infb113=data.frame(subgroup="noDiabetes",exposure="eGFR_scr",n=dim(data_eGFR_scr_nodiabetes)[1],event_n=sum(data_eGFR_scr_nodiabetes$event.type),
                   P_Cochrans_Q=ba13$p_tests[1,4],p_non_linear=ba13$p_tests[1,2])
#eGFR_scys

data=all_data_scys
#logistic_regression
fit = glm(diabetes~prs_eGFR_scys+BMI+age+smoking+alcohol+MET+education+TDI+cancer,data = data,family = 'binomial')
data$prob = predict(fit,type="response")

#delete 1%
quantile = quantile(data$prob,probs=seq(0,1,0.01))
index = which(data$prob<quantile[2])
data_reduced = data[-index,]
data=data_reduced
data_eGFR_scys_diabetes=data[which(data$diabetes==1),]
data_eGFR_scys_nodiabetes=data[which(data$diabetes==0),]

#All_infection
#data_eGFR_scys_diabetes
data_eGFR_scys_diabetes$scys_hat=lm(eGFR_scys_10~prs_eGFR_scys,data=data_eGFR_scys_diabetes)$fit
data_eGFR_scys_diabetes$scys_residual=data_eGFR_scys_diabetes$eGFR_scys_10-data_eGFR_scys_diabetes$scys_hat
N=10
data_eGFR_scys_diabetes$group=cut(data_eGFR_scys_diabetes$scys_residual, labels = F,quantile(data_eGFR_scys_diabetes$scys_residual, prob = 0:N / N, names = FALSE), include = TRUE)

data_eGFR_scys_diabetes$follow_up_time = data_eGFR_scys_diabetes$all_infection_time
data_eGFR_scys_diabetes$event.type=data_eGFR_scys_diabetes$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scys_diabetes[which(data_eGFR_scys_diabetes$group==i),]
  temp = summary(lm(eGFR_scys_10~prs_eGFR_scys+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scys_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scys+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeana23 = summary_stat$mean_expos # mean exposure in each quantile
xmeana23
xmean1 = xmeana23-4  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

aa23 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=4, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cys", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
aa23

aa23$coefficients
aa23$p_tests
aa23$p_heterogeneity  

fp_aa23=ifelse(aa23$p_tests[1,2]<0.001," < 0.001",ifelse(aa23$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa23$p_tests[1,2],3)))) 
Q_aa23=ifelse(aa23$p_tests[1,4]<0.001," < 0.001",ifelse(aa23$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa23$p_tests[1,4],3))))

fa213=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scys_diabetes$eGFR_scys, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scys_diabetes$eGFR_scys), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),aa23$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), aa23$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),aa23$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeana23)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cys]))) + 
  scale_x_continuous(breaks=seq(130,40,by = -10)) + 
  coord_cartesian(xlim =c(  130,40)) + ggtitle(expression(paste(eGFR["cys"]," and All Infections (With Diabetes)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_aa23) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_aa23)),  size = 3)

fa213
infa213=data.frame(subgroup="Diabetes",exposure="eGFR_scys",n=dim(data_eGFR_scys_diabetes)[1],event_n=sum(data_eGFR_scys_diabetes$event.type),
                   P_Cochrans_Q=aa23$p_tests[1,4],p_non_linear=aa23$p_tests[1,2])

#data_eGFR_scys_nodiabetes
data_eGFR_scys_nodiabetes$scys_hat=lm(eGFR_scys_10~prs_eGFR_scys,data=data_eGFR_scys_nodiabetes)$fit
data_eGFR_scys_nodiabetes$scys_residual=data_eGFR_scys_nodiabetes$eGFR_scys_10-data_eGFR_scys_nodiabetes$scys_hat
N=10
data_eGFR_scys_nodiabetes$group=cut(data_eGFR_scys_nodiabetes$scys_residual, labels = F,quantile(data_eGFR_scys_nodiabetes$scys_residual, prob = 0:N / N, names = FALSE), include = TRUE)

data_eGFR_scys_nodiabetes$follow_up_time = data_eGFR_scys_nodiabetes$all_infection_time
data_eGFR_scys_nodiabetes$event.type=data_eGFR_scys_nodiabetes$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scys_nodiabetes[which(data_eGFR_scys_nodiabetes$group==i),]
  temp = summary(lm(eGFR_scys_10~prs_eGFR_scys+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scys_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scys+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeanb23 = summary_stat$mean_expos # mean exposure in each quantile
xmeanb23
xmean1 = xmeanb23-5  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

ba23 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=5, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cys", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
ba23

ba23$coefficients
ba23$p_tests
ba23$p_heterogeneity  
fp_ba23=ifelse(ba23$p_tests[1,2]<0.001," < 0.001",ifelse(ba23$p_tests[1,2]>0.999," > 0.999",paste(" =",round(ba23$p_tests[1,2],3)))) 
Q_ba23=ifelse(ba23$p_tests[1,4]<0.001," < 0.001",ifelse(ba23$p_tests[1,4]>0.999," > 0.999",paste(" =",round(ba23$p_tests[1,4],3))))

fb213=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scys_nodiabetes$eGFR_scys, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scys_nodiabetes$eGFR_scys), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),ba23$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), ba23$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),ba23$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeanb23)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cys]))) + 
  scale_x_continuous(breaks=seq(130,40,by = -10)) + 
  coord_cartesian(xlim =c(  130,40)) + ggtitle(expression(paste(eGFR["cys"]," and All Infections (Without Diabetes)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_ba23) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_ba23)),  size = 3)

fb213
infb213=data.frame(subgroup="noDiabetes",exposure="eGFR_scys",n=dim(data_eGFR_scys_nodiabetes)[1],event_n=sum(data_eGFR_scys_nodiabetes$event.type),
                   P_Cochrans_Q=ba23$p_tests[1,4],p_non_linear=ba23$p_tests[1,2])


#BUN
data=all_data_BUN
#logistic_regression
fit = glm(diabetes~prs_BUN+BMI+age+smoking+alcohol+MET+education+TDI+cancer,data = data,family = 'binomial')
data$prob = predict(fit,type="response")

#delete 1%
quantile = quantile(data$prob,probs=seq(0,1,0.01))
index = which(data$prob<quantile[2])
data_reduced = data[-index,]
data=data_reduced
data_BUN_diabetes=data[which(data$diabetes==1),]
data_BUN_nodiabetes=data[which(data$diabetes==0),]

#All_infection
#data_BUN_diabetes
data_BUN_diabetes$BUN_hat=lm(BUN~prs_BUN,data=data_BUN_diabetes)$fit
data_BUN_diabetes$BUN_residual=data_BUN_diabetes$BUN-data_BUN_diabetes$BUN_hat
N=10
data_BUN_diabetes$group=cut(data_BUN_diabetes$BUN_residual, labels = F,quantile(data_BUN_diabetes$BUN_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_BUN_diabetes$follow_up_time = data_BUN_diabetes$all_infection_time
data_BUN_diabetes$event.type=data_BUN_diabetes$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_BUN_diabetes[which(data_BUN_diabetes$group==i),]
  temp = summary(lm(BUN~prs_BUN+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$BUN)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_BUN+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeana33 = summary_stat$mean_expos # mean exposure in each quantile
xmeana33
xmean1 = xmeana33-3  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

aa33 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=3, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="BUN", pref_y="Hazard ratio of all infections", breaks=c(0.25,0.5,1,2,4,8))
aa33

aa33$coefficients
aa33$p_tests
aa33$p_heterogeneity  
fp_aa33=ifelse(aa33$p_tests[1,2]<0.001," < 0.001",ifelse(aa33$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa33$p_tests[1,2],3)))) 
Q_aa33=ifelse(aa33$p_tests[1,4]<0.001," < 0.001",ifelse(aa33$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa33$p_tests[1,4],3))))

fa313=ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), data_BUN_diabetes,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,12,1))+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x, y=log2(yest)),aa33$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x, y=log2(lci)), aa33$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x, y=log2(uci)),aa33$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeana33), y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 4),  labels =c(0.25,0.5,1,2,4,8,16),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~BUN))) + 
  scale_x_continuous(breaks=seq(0,12,by = 1)) + 
  coord_cartesian(xlim =c(0,12)) + ggtitle("BUN and All Infections (With Diabetes)")+
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
  annotate("text", x = 10, y = 3.7, label = bquote(P['Cochran\'s Q test'] ~ .(Q_aa33) ), size = 3) + 
  annotate("text", x = 10, y = 3.4, label = bquote(P['nonlinearity'] ~.(fp_aa33)),  size = 3)

fa313
infa313=data.frame(subgroup="Diabetes",exposure="BUN",n=dim(data_BUN_diabetes)[1],event_n=sum(data_BUN_diabetes$event.type),
                   P_Cochrans_Q=aa33$p_tests[1,4],p_non_linear=aa33$p_tests[1,2])

#data_BUN_nodiabetes
data_BUN_nodiabetes$BUN_hat=lm(BUN~prs_BUN,data=data_BUN_nodiabetes)$fit
data_BUN_nodiabetes$BUN_residual=data_BUN_nodiabetes$BUN-data_BUN_nodiabetes$BUN_hat
N=10
data_BUN_nodiabetes$group=cut(data_BUN_nodiabetes$BUN_residual, labels = F,quantile(data_BUN_nodiabetes$BUN_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_BUN_nodiabetes$follow_up_time = data_BUN_nodiabetes$all_infection_time
data_BUN_nodiabetes$event.type=data_BUN_nodiabetes$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_BUN_nodiabetes[which(data_BUN_nodiabetes$group==i),]
  temp = summary(lm(BUN~prs_BUN+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$BUN)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_BUN+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeanb33 = summary_stat$mean_expos # mean exposure in each quantile
xmeanb33
xmean1 = xmeanb33-3  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

ba33 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=3, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="BUN", pref_y="Hazard ratio of all infections", breaks=c(0.25,0.5,1,2,4,8))
ba33

ba33$coefficients
ba33$p_tests
ba33$p_heterogeneity  
fp_ba33=ifelse(ba33$p_tests[1,2]<0.001," < 0.001",ifelse(ba33$p_tests[1,2]>0.999," > 0.999",paste(" =",round(ba33$p_tests[1,2],3)))) 
Q_ba33=ifelse(ba33$p_tests[1,4]<0.001," < 0.001",ifelse(ba33$p_tests[1,4]>0.999," > 0.999",paste(" =",round(ba33$p_tests[1,4],3))))

fb313=ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), data_BUN_nodiabetes,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,12,1))+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x, y=log2(yest)),ba33$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x, y=log2(lci)), ba33$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x, y=log2(uci)),ba33$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeanb33), y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 4),  labels =c(0.25,0.5,1,2,4,8,16),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~BUN))) + 
  scale_x_continuous(breaks=seq(0,12,by = 1)) + 
  coord_cartesian(xlim =c(0,12)) + ggtitle("BUN and All Infections (Without Diabetes)")+
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
  annotate("text", x = 10, y = 3.7, label = bquote(P['Cochran\'s Q test'] ~ .(Q_ba33) ), size = 3) + 
  annotate("text", x = 10, y = 3.4, label = bquote(P['nonlinearity'] ~.(fp_ba33)),  size = 3)

fb313
infb313=data.frame(subgroup="noDiabetes",exposure="BUN",n=dim(data_BUN_nodiabetes)[1],event_n=sum(data_BUN_nodiabetes$event.type),
                   P_Cochrans_Q=ba33$p_tests[1,4],p_non_linear=ba33$p_tests[1,2])
inf_sub2=rbind(infa113,infb113,infa213,infb213,infa313,infb313)
##############################################subgroup5#################################################################
#subgroup5 cancer

#eGFR_scr
data=all_data_scr

#logistic_regression
fit = glm(cancer~prs_eGFR_scr+BMI+age+sex+smoking+alcohol+MET+education+TDI,data = data,family = 'binomial')
data$prob = predict(fit,type="response")

#delete 1%
quantile = quantile(data$prob,probs=seq(0,1,0.01))
index = which(data$prob<quantile[2])
data_reduced = data[-index,]
data=data_reduced
data_eGFR_scr_cancer=data[which(data$cancer==1),]
data_eGFR_scr_nocancer=data[which(data$cancer==0),]


#All_infection
#data_eGFR_scr_cancer
data_eGFR_scr_cancer$scr_hat=lm(eGFR_scr_10~prs_eGFR_scr,data=data_eGFR_scr_cancer)$fit
data_eGFR_scr_cancer$scr_residual=data_eGFR_scr_cancer$eGFR_scr_10-data_eGFR_scr_cancer$scr_hat
N=10
data_eGFR_scr_cancer$group=cut(data_eGFR_scr_cancer$scr_residual, labels = F,quantile(data_eGFR_scr_cancer$scr_residual, prob = 0:N / N, names = FALSE), include = TRUE)

data_eGFR_scr_cancer$follow_up_time = data_eGFR_scr_cancer$all_infection_time
data_eGFR_scr_cancer$event.type=data_eGFR_scr_cancer$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scr_cancer[which(data_eGFR_scr_cancer$group==i),]
  temp = summary(lm(eGFR_scr_10~prs_eGFR_scr+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scr_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scr+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeana04 = summary_stat$mean_expos # mean exposure in each quantile
xmeana04
xmean1 = xmeana04-5  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

aa14 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=5, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cr", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
aa14

aa14$coefficients
aa14$p_tests
aa14$p_heterogeneity  
fp_aa14=ifelse(aa14$p_tests[1,2]<0.001," < 0.001",ifelse(aa14$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa14$p_tests[1,2],3)))) 
Q_aa14=ifelse(aa14$p_tests[1,4]<0.001," < 0.001",ifelse(aa14$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa14$p_tests[1,4],3))))

fa114=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scr_cancer$eGFR_scr, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scr_cancer$eGFR_scr), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),aa14$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), aa14$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),aa14$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeana04)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cr]))) + 
  scale_x_continuous(breaks=seq(130,40,by = -10)) + 
  coord_cartesian(xlim =c(  130,40)) + ggtitle(expression(paste(eGFR["cr"]," and All Infections (With Cancer)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_aa14) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_aa14)),  size = 3)

fa114
infa114=data.frame(subgroup="Cancer",exposure="eGFR_scr",n=dim(data_eGFR_scr_cancer)[1],event_n=sum(data_eGFR_scr_cancer$event.type),
                   P_Cochrans_Q=aa14$p_tests[1,4],p_non_linear=aa14$p_tests[1,2])

#data_eGFR_scr_nocancer
data_eGFR_scr_nocancer$scr_hat=lm(eGFR_scr_10~prs_eGFR_scr,data=data_eGFR_scr_nocancer)$fit
data_eGFR_scr_nocancer$scr_residual=data_eGFR_scr_nocancer$eGFR_scr_10-data_eGFR_scr_nocancer$scr_hat
N=10
data_eGFR_scr_nocancer$group=cut(data_eGFR_scr_nocancer$scr_residual, labels = F,quantile(data_eGFR_scr_nocancer$scr_residual, prob = 0:N / N, names = FALSE), include = TRUE)

data_eGFR_scr_nocancer$follow_up_time = data_eGFR_scr_nocancer$all_infection_time
data_eGFR_scr_nocancer$event.type=data_eGFR_scr_nocancer$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scr_nocancer[which(data_eGFR_scr_nocancer$group==i),]
  temp = summary(lm(eGFR_scr_10~prs_eGFR_scr+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scr_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scr+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeanb04 = summary_stat$mean_expos # mean exposure in each quantile
xmeanb04
xmean1 = xmeanb04-6  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

ba14 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=6, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cr", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
ba14

ba14$coefficients
ba14$p_tests
ba14$p_heterogeneity  
fp_ba14=ifelse(ba14$p_tests[1,2]<0.001," < 0.001",ifelse(ba14$p_tests[1,2]>0.999," > 0.999",paste(" =",round(ba14$p_tests[1,2],3)))) 
Q_ba14=ifelse(ba14$p_tests[1,4]<0.001," < 0.001",ifelse(ba14$p_tests[1,4]>0.999," > 0.999",paste(" =",round(ba14$p_tests[1,4],3))))

fb114=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scr_nocancer$eGFR_scr, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scr_nocancer$eGFR_scr), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),ba14$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), ba14$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),ba14$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeanb04)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cr]))) + 
  scale_x_continuous(breaks=seq(130,40,by = -10)) + 
  coord_cartesian(xlim =c(  130,40)) + ggtitle(expression(paste(eGFR["cr"]," and All Infections (Without Cancer)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_ba14) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_ba14)),  size = 3)

fb114
infb114=data.frame(subgroup="noCancer",exposure="eGFR_scr",n=dim(data_eGFR_scr_nocancer)[1],event_n=sum(data_eGFR_scr_nocancer$event.type),
                   P_Cochrans_Q=ba14$p_tests[1,4],p_non_linear=ba14$p_tests[1,2])

#eGFR_scys

data=all_data_scys
#logistic_regression
fit = glm(cancer~prs_eGFR_scys+BMI+age+smoking+alcohol+MET+education+TDI,data = data,family = 'binomial')
data$prob = predict(fit,type="response")

#delete 1%
quantile = quantile(data$prob,probs=seq(0,1,0.01))
index = which(data$prob<quantile[2])
data_reduced = data[-index,]
data=data_reduced
data_eGFR_scys_cancer=data[which(data$cancer==1),]
data_eGFR_scys_nocancer=data[which(data$cancer==0),]


#All_infection
#data_eGFR_scr_cancer
data_eGFR_scys_cancer$scys_hat=lm(eGFR_scys_10~prs_eGFR_scys,data=data_eGFR_scys_cancer)$fit
data_eGFR_scys_cancer$scys_residual=data_eGFR_scys_cancer$eGFR_scys_10-data_eGFR_scys_cancer$scys_hat
N=10
data_eGFR_scys_cancer$group=cut(data_eGFR_scys_cancer$scys_residual, labels = F,quantile(data_eGFR_scys_cancer$scys_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_eGFR_scys_cancer$follow_up_time = data_eGFR_scys_cancer$all_infection_time
data_eGFR_scys_cancer$event.type=data_eGFR_scys_cancer$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scys_cancer[which(data_eGFR_scys_cancer$group==i),]
  temp = summary(lm(eGFR_scys_10~prs_eGFR_scys+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scys_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scys+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeana24 = summary_stat$mean_expos # mean exposure in each quantile
xmeana24
xmean1 = xmeana24-5  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

aa24 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=5, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cys", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
aa24

aa24$coefficients
aa24$p_tests
aa24$p_heterogeneity  
fp_aa24=ifelse(aa24$p_tests[1,2]<0.001," < 0.001",ifelse(aa24$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa24$p_tests[1,2],3)))) 
Q_aa24=ifelse(aa24$p_tests[1,4]<0.001," < 0.001",ifelse(aa24$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa24$p_tests[1,4],3))))

fa214=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scys_cancer$eGFR_scys, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scys_cancer$eGFR_scys), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),aa24$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), aa24$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),aa24$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeana24)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cys]))) + 
  scale_x_continuous(breaks=seq(130,40,by = -10)) + 
  coord_cartesian(xlim =c(  130,40)) + ggtitle(expression(paste(eGFR["cys"]," and All Infections (With Cancer)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_aa24) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_aa24)),  size = 3)

fa214
infa214=data.frame(subgroup="Cancer",exposure="eGFR_scys",n=dim(data_eGFR_scys_cancer)[1],event_n=sum(data_eGFR_scys_cancer$event.type),
                   P_Cochrans_Q=aa24$p_tests[1,4],p_non_linear=aa24$p_tests[1,2])

#data_eGFR_scr_nocancer
data_eGFR_scys_nocancer$scys_hat=lm(eGFR_scys_10~prs_eGFR_scys,data=data_eGFR_scys_nocancer)$fit
data_eGFR_scys_nocancer$scys_residual=data_eGFR_scys_nocancer$eGFR_scys_10-data_eGFR_scys_nocancer$scys_hat
N=10
data_eGFR_scys_nocancer$group=cut(data_eGFR_scys_nocancer$scys_residual, labels = F,quantile(data_eGFR_scys_nocancer$scys_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_eGFR_scys_nocancer$follow_up_time = data_eGFR_scys_nocancer$all_infection_time
data_eGFR_scys_nocancer$event.type=data_eGFR_scys_nocancer$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scys_nocancer[which(data_eGFR_scys_nocancer$group==i),]
  temp = summary(lm(eGFR_scys_10~prs_eGFR_scys+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scys_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scys+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeanb24 = summary_stat$mean_expos # mean exposure in each quantile
xmeanb24
xmean1 = xmeanb24-5  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

ba24 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=5, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cys", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
ba24

ba24$coefficients
ba24$p_tests
ba24$p_heterogeneity  
fp_ba24=ifelse(ba24$p_tests[1,2]<0.001," < 0.001",ifelse(ba24$p_tests[1,2]>0.999," > 0.999",paste(" =",round(ba24$p_tests[1,2],3)))) 
Q_ba24=ifelse(ba24$p_tests[1,4]<0.001," < 0.001",ifelse(ba24$p_tests[1,4]>0.999," > 0.999",paste(" =",round(ba24$p_tests[1,4],3))))

fb214=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scys_nocancer$eGFR_scys, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scys_nocancer$eGFR_scys), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),ba24$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), ba24$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),ba24$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeanb24)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cys]))) + 
  scale_x_continuous(breaks=seq(130,40,by = -10)) + 
  coord_cartesian(xlim =c(  130,40)) + ggtitle(expression(paste(eGFR["cys"]," and All Infections (Without Cancer)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_ba24) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_ba24)),  size = 3)

fb214
infb214=data.frame(subgroup="noCancer",exposure="eGFR_scys",n=dim(data_eGFR_scys_nocancer)[1],event_n=sum(data_eGFR_scys_nocancer$event.type),
                   P_Cochrans_Q=ba24$p_tests[1,4],p_non_linear=ba24$p_tests[1,2])

#BUN
data=all_data_BUN
#logistic_regression
fit = glm(cancer~prs_BUN+BMI+age+smoking+alcohol+MET+education+TDI,data = data,family = 'binomial')
data$prob = predict(fit,type="response")

#delete 1%
quantile = quantile(data$prob,probs=seq(0,1,0.01))
index = which(data$prob<quantile[2])
data_reduced = data[-index,]
data=data_reduced
data_BUN_cancer=data[which(data$cancer==1),]
data_BUN_nocancer=data[which(data$cancer==0),]


#All_infection
#data_BUN_cancer
data_BUN_cancer$BUN_hat=lm(BUN~prs_BUN,data=data_BUN_cancer)$fit
data_BUN_cancer$BUN_residual=data_BUN_cancer$BUN-data_BUN_cancer$BUN_hat
N=10
data_BUN_cancer$group=cut(data_BUN_cancer$BUN_residual, labels = F,quantile(data_BUN_cancer$BUN_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_BUN_cancer$follow_up_time = data_BUN_cancer$all_infection_time
data_BUN_cancer$event.type=data_BUN_cancer$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_BUN_cancer[which(data_BUN_cancer$group==i),]
  temp = summary(lm(BUN~prs_BUN+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$BUN)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_BUN+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeana34 = summary_stat$mean_expos # mean exposure in each quantile
xmeana34
xmean1 = xmeana34-3  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

aa34 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=3, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="BUN", pref_y="Hazard ratio of all infections", breaks=c(0.25,0.5,1,2,4,8))
aa34

aa34$coefficients
aa34$p_tests
aa34$p_heterogeneity  
fp_aa34=ifelse(aa34$p_tests[1,2]<0.001," < 0.001",ifelse(aa34$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa34$p_tests[1,2],3)))) 
Q_aa34=ifelse(aa34$p_tests[1,4]<0.001," < 0.001",ifelse(aa34$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa34$p_tests[1,4],3))))

fa314=ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), data_BUN_cancer,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,12,1))+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x, y=log2(yest)),aa34$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x, y=log2(lci)), aa34$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x, y=log2(uci)),aa34$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeana34), y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 4),  labels =c(0.25,0.5,1,2,4,8,16),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~BUN))) + 
  scale_x_continuous(breaks=seq(0,12,by = 1)) + 
  coord_cartesian(xlim =c(0,12)) + ggtitle("BUN and All Infections (With Cancer)")+
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
  annotate("text", x = 10, y = 3.7, label = bquote(P['Cochran\'s Q test'] ~ .(Q_aa34) ), size = 3) + 
  annotate("text", x = 10, y = 3.4, label = bquote(P['nonlinearity'] ~.(fp_aa34)),  size = 3)

fa314
infa314=data.frame(subgroup="Cancer",exposure="BUN",n=dim(data_BUN_cancer)[1],event_n=sum(data_BUN_cancer$event.type),
                   P_Cochrans_Q=aa34$p_tests[1,4],p_non_linear=aa34$p_tests[1,2])

#data_BUN_nocancer
data_BUN_nocancer$BUN_hat=lm(BUN~prs_BUN,data=data_BUN_nocancer)$fit
data_BUN_nocancer$BUN_residual=data_BUN_nocancer$BUN-data_BUN_nocancer$BUN_hat
N=10
data_BUN_nocancer$group=cut(data_BUN_nocancer$BUN_residual, labels = F,quantile(data_BUN_nocancer$BUN_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_BUN_nocancer$follow_up_time = data_BUN_nocancer$all_infection_time
data_BUN_nocancer$event.type=data_BUN_nocancer$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_BUN_nocancer[which(data_BUN_nocancer$group==i),]
  temp = summary(lm(BUN~prs_BUN+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$BUN)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_BUN+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeanb34 = summary_stat$mean_expos # mean exposure in each quantile
xmeanb34
xmean1 = xmeanb34-3  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

ba34 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=3, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="BUN", pref_y="Hazard ratio of all infections", breaks=c(0.25,0.5,1,2,4,8))
ba34

ba34$coefficients
ba34$p_tests
ba34$p_heterogeneity  
fp_ba34=ifelse(ba34$p_tests[1,2]<0.001," < 0.001",ifelse(ba34$p_tests[1,2]>0.999," > 0.999",paste(" =",round(ba34$p_tests[1,2],3)))) 
Q_ba34=ifelse(ba34$p_tests[1,4]<0.001," < 0.001",ifelse(ba34$p_tests[1,4]>0.999," > 0.999",paste(" =",round(ba34$p_tests[1,4],3))))

fb314=ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), data_BUN_nocancer,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,12,1))+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x, y=log2(yest)),ba34$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x, y=log2(lci)), ba34$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x, y=log2(uci)),ba34$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeanb34), y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 4),  labels =c(0.25,0.5,1,2,4,8,16),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~BUN))) + 
  scale_x_continuous(breaks=seq(0,12,by = 1)) + 
  coord_cartesian(xlim =c(0,12)) + ggtitle("BUN and All Infections (Without Cancer)")+
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
  annotate("text", x = 10, y = 3.7, label = bquote(P['Cochran\'s Q test'] ~ .(Q_ba34) ), size = 3) + 
  annotate("text", x = 10, y = 3.4, label = bquote(P['nonlinearity'] ~.(fp_ba34)),  size = 3)

fb314
infb314=data.frame(subgroup="noCancer",exposure="BUN",n=dim(data_BUN_nocancer)[1],event_n=sum(data_BUN_nocancer$event.type),
                   P_Cochrans_Q=ba34$p_tests[1,4],p_non_linear=ba34$p_tests[1,2])
inf_sub3=rbind(infa114,infb114,infa214,infb214,infa314,infb314)
##############################################subgroup6####################################################################
#subgroup6 CVD

#eGFR_scr
data=all_data_scr

#logistic_regression
fit = glm(CVD~prs_eGFR_scr+BMI+age+sex+smoking+alcohol+MET+education+TDI+cancer,data = data,family = 'binomial')
data$prob = predict(fit,type="response")

data_eGFR_scr_CVD=data[which(data$CVD==1),]
data_eGFR_scr_noCVD=data[which(data$CVD==0),]


#All_infection
#data_eGFR_scr_CVD
data_eGFR_scr_CVD$scr_hat=lm(eGFR_scr_10~prs_eGFR_scr,data=data_eGFR_scr_CVD)$fit
data_eGFR_scr_CVD$scr_residual=data_eGFR_scr_CVD$eGFR_scr_10-data_eGFR_scr_CVD$scr_hat
N=10
data_eGFR_scr_CVD$group=cut(data_eGFR_scr_CVD$scr_residual, labels = F,quantile(data_eGFR_scr_CVD$scr_residual, prob = 0:N / N, names = FALSE), include = TRUE)

data_eGFR_scr_CVD$follow_up_time = data_eGFR_scr_CVD$all_infection_time
data_eGFR_scr_CVD$event.type=data_eGFR_scr_CVD$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scr_CVD[which(data_eGFR_scr_CVD$group==i),]
  temp = summary(lm(eGFR_scr_10~prs_eGFR_scr+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scr_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scr+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeana05 = summary_stat$mean_expos # mean exposure in each quantile
xmeana05
xmean1 = xmeana05-6  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

aa15 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=6, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cr", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
aa15

aa15$coefficients
aa15$p_tests
aa15$p_heterogeneity  

fp_aa15=ifelse(aa15$p_tests[1,2]<0.001," < 0.001",ifelse(aa15$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa15$p_tests[1,2],3))))
Q_aa15=ifelse(aa15$p_tests[1,4]<0.001," < 0.001",ifelse(aa15$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa15$p_tests[1,4],3))))

fa115=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scr_CVD$eGFR_scr, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scr_CVD$eGFR_scr), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),aa15$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), aa15$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),aa15$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeana05)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cr]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cr"]," and All Infections (With CVD)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_aa15) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_aa15)),  size = 3)

fa115
infa115=data.frame(subgroup="CVD",exposure="eGFR_cr",n=dim(data_eGFR_scr_CVD)[1],event_n=sum(data_eGFR_scr_CVD$event.type),
                   P_Cochrans_Q=aa15$p_tests[1,4],p_non_linear=aa15$p_tests[1,2])

#data_eGFR_scr_noCVD
data_eGFR_scr_noCVD$scr_hat=lm(eGFR_scr_10~prs_eGFR_scr,data=data_eGFR_scr_noCVD)$fit
data_eGFR_scr_noCVD$scr_residual=data_eGFR_scr_noCVD$eGFR_scr_10-data_eGFR_scr_noCVD$scr_hat
N=10
data_eGFR_scr_noCVD$group=cut(data_eGFR_scr_noCVD$scr_residual, labels = F,quantile(data_eGFR_scr_noCVD$scr_residual, prob = 0:N / N, names = FALSE), include = TRUE)

data_eGFR_scr_noCVD$follow_up_time = data_eGFR_scr_noCVD$all_infection_time
data_eGFR_scr_noCVD$event.type=data_eGFR_scr_noCVD$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scr_noCVD[which(data_eGFR_scr_noCVD$group==i),]
  temp = summary(lm(eGFR_scr_10~prs_eGFR_scr+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scr_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scr+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeanb05 = summary_stat$mean_expos # mean exposure in each quantile
xmeanb05
xmean1 = xmeanb05-6  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

ba15 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=6, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cr", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
ba15

ba15$coefficients
ba15$p_tests
ba15$p_heterogeneity

fp_ba15=ifelse(ba15$p_tests[1,2]<0.001," < 0.001",ifelse(ba15$p_tests[1,2]>0.999," > 0.999",paste(" =",round(ba15$p_tests[1,2],3))))
Q_ba15=ifelse(ba15$p_tests[1,4]<0.001," < 0.001",ifelse(ba15$p_tests[1,4]>0.999," > 0.999",paste(" =",round(ba15$p_tests[1,4],3))))

fb115=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scr_noCVD$eGFR_scr, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scr_noCVD$eGFR_scr), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),ba15$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), ba15$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),ba15$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeanb05)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cr]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cr"]," and All Infections (Without CVD)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_ba15) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_ba15)),  size = 3)

fb115
infb115=data.frame(subgroup="noCVD",exposure="eGFR_cr",n=dim(data_eGFR_scr_noCVD)[1],event_n=sum(data_eGFR_scr_noCVD$event.type),
                   P_Cochrans_Q=ba15$p_tests[1,4],p_non_linear=ba15$p_tests[1,2])

#eGFR_scys

data=all_data_scys
#logistic_regression
fit = glm(CVD~prs_eGFR_scys+age+smoking+alcohol+MET+education+TDI+cancer,data = data,family = 'binomial')
data$prob = predict(fit,type="response")

#delete 1%
data_eGFR_scys_CVD=data[which(data$CVD==1),]
data_eGFR_scys_noCVD=data[which(data$CVD==0),]


#All_infection
#data_eGFR_scys_CVD
data_eGFR_scys_CVD$scys_hat=lm(eGFR_scys_10~prs_eGFR_scys,data=data_eGFR_scys_CVD)$fit
data_eGFR_scys_CVD$scys_residual=data_eGFR_scys_CVD$eGFR_scys_10-data_eGFR_scys_CVD$scys_hat
N=10
data_eGFR_scys_CVD$group=cut(data_eGFR_scys_CVD$scys_residual, labels = F,quantile(data_eGFR_scys_CVD$scys_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_eGFR_scys_CVD$follow_up_time = data_eGFR_scys_CVD$all_infection_time
data_eGFR_scys_CVD$event.type=data_eGFR_scys_CVD$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scys_CVD[which(data_eGFR_scys_CVD$group==i),]
  temp = summary(lm(eGFR_scys_10~prs_eGFR_scys+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scys_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scys+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeana25 = summary_stat$mean_expos # mean exposure in each quantile
xmeana25
xmean1 = xmeana25-6  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

aa25 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=6, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cys", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
aa25

aa25$coefficients
aa25$p_tests
aa25$p_heterogeneity


fp_aa25=ifelse(aa25$p_tests[1,2]<0.001," < 0.001",ifelse(aa25$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa25$p_tests[1,2],3))))
Q_aa25=ifelse(aa25$p_tests[1,4]<0.001," < 0.001",ifelse(aa25$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa25$p_tests[1,4],3))))

fa215=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scys_CVD$eGFR_scys, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scys_CVD$eGFR_scys), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),aa25$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), aa25$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),aa25$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeana25)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cys]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cys"]," and All Infections (With CVD)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_aa25) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_aa25)),  size = 3)

fa215
infa215=data.frame(subgroup="CVD",exposure="eGFR_cys",n=dim(data_eGFR_scys_CVD)[1],event_n=sum(data_eGFR_scys_CVD$event.type),
                   P_Cochrans_Q=aa25$p_tests[1,4],p_non_linear=aa25$p_tests[1,2])

#data_eGFR_scys_noCVD
data_eGFR_scys_noCVD$scys_hat=lm(eGFR_scys_10~prs_eGFR_scys,data=data_eGFR_scys_noCVD)$fit
data_eGFR_scys_noCVD$scys_residual=data_eGFR_scys_noCVD$eGFR_scys_10-data_eGFR_scys_noCVD$scys_hat
N=10
data_eGFR_scys_noCVD$group=cut(data_eGFR_scys_noCVD$scys_residual, labels = F,quantile(data_eGFR_scys_noCVD$scys_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_eGFR_scys_noCVD$follow_up_time = data_eGFR_scys_noCVD$all_infection_time
data_eGFR_scys_noCVD$event.type=data_eGFR_scys_noCVD$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scys_noCVD[which(data_eGFR_scys_noCVD$group==i),]
  temp = summary(lm(eGFR_scys_10~prs_eGFR_scys+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scys_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scys+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeanb25 = summary_stat$mean_expos # mean exposure in each quantile
xmeanb25
xmean1 = xmeanb25-5  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

ba25 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=5, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cys", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
ba25

ba25$coefficients
ba25$p_tests
ba25$p_heterogeneity

fp_ba25=ifelse(ba25$p_tests[1,2]<0.001," < 0.001",ifelse(ba25$p_tests[1,2]>0.999," > 0.999",paste(" =",round(ba25$p_tests[1,2],3))))
Q_ba25=ifelse(ba25$p_tests[1,4]<0.001," < 0.001",ifelse(ba25$p_tests[1,4]>0.999," > 0.999",paste(" =",round(ba25$p_tests[1,4],3))))

fb215=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scys_noCVD$eGFR_scys, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scys_noCVD$eGFR_scys), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),ba25$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), ba25$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),ba25$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeanb25)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cys]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cys"]," and All Infections (Without CVD)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_ba25) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_ba25)),  size = 3)

fb215
infb215=data.frame(subgroup="noCVD",exposure="eGFR_cys",n=dim(data_eGFR_scys_noCVD)[1],event_n=sum(data_eGFR_scys_noCVD$event.type),
                   P_Cochrans_Q=ba25$p_tests[1,4],p_non_linear=ba25$p_tests[1,2])

#BUN
data=all_data_BUN
#logistic_regression
fit = glm(CVD~prs_BUN+age+BMI+smoking+alcohol+MET+education+TDI+cancer,data = data,family = 'binomial')
data$prob = predict(fit,type="response")

#delete 1%
data_BUN_CVD=data[which(data$CVD==1),]
data_BUN_noCVD=data[which(data$CVD==0),]


#All_infection
#data_BUN_CVD
data_BUN_CVD$BUN_hat=lm(BUN~prs_BUN,data=data_BUN_CVD)$fit
data_BUN_CVD$BUN_residual=data_BUN_CVD$BUN-data_BUN_CVD$BUN_hat
N=10
data_BUN_CVD$group=cut(data_BUN_CVD$BUN_residual, labels = F,quantile(data_BUN_CVD$BUN_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_BUN_CVD$follow_up_time = data_BUN_CVD$all_infection_time
data_BUN_CVD$event.type=data_BUN_CVD$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_BUN_CVD[which(data_BUN_CVD$group==i),]
  temp = summary(lm(BUN~prs_BUN+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$BUN)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_BUN+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeana35 = summary_stat$mean_expos # mean exposure in each quantile
xmeana35
xmean1 = xmeana35-3  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

aa35 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=3, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="BUN", pref_y="Hazard ratio of all infections", breaks=c(0.25,0.5,1,2,4,8))
aa35

aa35$coefficients
aa35$p_tests
aa35$p_heterogeneity

fp_aa35=ifelse(aa35$p_tests[1,2]<0.001," < 0.001",ifelse(aa35$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa35$p_tests[1,2],3))))
Q_aa35=ifelse(aa35$p_tests[1,4]<0.001," < 0.001",ifelse(aa35$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa35$p_tests[1,4],3))))

fa315=ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), data_BUN_CVD,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,12,1))+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x, y=log2(yest)),aa35$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x, y=log2(lci)), aa35$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x, y=log2(uci)),aa35$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeana35), y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 4),  labels =c(0.25,0.5,1,2,4,8,16),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~BUN))) + 
  scale_x_continuous(breaks=seq(0,12,by = 1)) + 
  coord_cartesian(xlim =c(0,12)) + ggtitle("BUN and All Infections (With CVD)")+
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
  annotate("text", x = 10, y = 3.7, label = bquote(P['Cochran\'s Q test'] ~ .(Q_aa35) ), size = 3) + 
  annotate("text", x = 10, y = 3.4, label = bquote(P['nonlinearity'] ~.(fp_aa35)),  size = 3)

fa315
infa315=data.frame(subgroup="CVD",exposure="BUN",n=dim(data_BUN_CVD)[1],event_n=sum(data_BUN_CVD$event.type),
                   P_Cochrans_Q=aa35$p_tests[1,4],p_non_linear=aa35$p_tests[1,2])

#data_BUN_noCVD
data_BUN_noCVD$BUN_hat=lm(BUN~prs_BUN,data=data_BUN_noCVD)$fit
data_BUN_noCVD$BUN_residual=data_BUN_noCVD$BUN-data_BUN_noCVD$BUN_hat
N=10
data_BUN_noCVD$group=cut(data_BUN_noCVD$BUN_residual, labels = F,quantile(data_BUN_noCVD$BUN_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_BUN_noCVD$follow_up_time = data_BUN_noCVD$all_infection_time
data_BUN_noCVD$event.type=data_BUN_noCVD$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_BUN_noCVD[which(data_BUN_noCVD$group==i),]
  temp = summary(lm(BUN~prs_BUN+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$BUN)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_BUN+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeanb35 = summary_stat$mean_expos # mean exposure in each quantile
xmeanb35
xmean1 = xmeanb35-3  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

ba35 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=3, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="BUN", pref_y="Hazard ratio of all infections", breaks=c(0.25,0.5,1,2,4,8))
ba35

ba35$coefficients
ba35$p_tests
ba35$p_heterogeneity

fp_ba35=ifelse(ba35$p_tests[1,2]<0.001," < 0.001",ifelse(ba35$p_tests[1,2]>0.999," > 0.999",paste(" =",round(ba35$p_tests[1,2],3))))
Q_ba35=ifelse(ba35$p_tests[1,4]<0.001," < 0.001",ifelse(ba35$p_tests[1,4]>0.999," > 0.999",paste(" =",round(ba35$p_tests[1,4],3))))

fb315=ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), data_BUN_noCVD,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,12,1))+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x, y=log2(yest)),ba35$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x, y=log2(lci)), ba35$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x, y=log2(uci)),ba35$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeanb35), y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 4),  labels =c(0.25,0.5,1,2,4,8,16),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~BUN))) + 
  scale_x_continuous(breaks=seq(0,12,by = 1)) + 
  coord_cartesian(xlim =c(0,12)) + ggtitle("BUN and All Infections (Without CVD)")+
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
  annotate("text", x = 10, y = 3.7, label = bquote(P['Cochran\'s Q test'] ~ .(Q_ba35) ), size = 3) +
  annotate("text", x = 10, y = 3.4, label = bquote(P['nonlinearity'] ~.(fp_ba35)),  size = 3)

fb315
infb315=data.frame(subgroup="noCVD",exposure="BUN",n=dim(data_BUN_noCVD)[1],event_n=sum(data_BUN_noCVD$event.type),
                   P_Cochrans_Q=ba35$p_tests[1,4],p_non_linear=ba35$p_tests[1,2])
inf_sub4=rbind(infa115,infb115,infa215,infb215,infa315,infb315)
##############################################subgroup7####################################################################
#subgroup7 UACR

#eGFR_scr
data=all_data_scr
data$UACR_2=1
data$UACR_2[which(data$UACR_cat==0)]=0


#logistic_regression
fit = glm(UACR_2~prs_eGFR_scr+BMI+age+sex+smoking+alcohol+MET+education+TDI+cancer,data = data,family = 'binomial')
data$prob = predict(fit,type="response")

data_eGFR_scr_UACR_u30=data[which(data$UACR_2==1),]
data_eGFR_scr_UACR_l30=data[which(data$UACR_2==0),]


#All_infection
#data_eGFR_scr_UACR_u30
data_eGFR_scr_UACR_u30$scr_hat=lm(eGFR_scr_10~prs_eGFR_scr,data=data_eGFR_scr_UACR_u30)$fit
data_eGFR_scr_UACR_u30$scr_residual=data_eGFR_scr_UACR_u30$eGFR_scr_10-data_eGFR_scr_UACR_u30$scr_hat
N=10
data_eGFR_scr_UACR_u30$group=cut(data_eGFR_scr_UACR_u30$scr_residual, labels = F,quantile(data_eGFR_scr_UACR_u30$scr_residual, prob = 0:N / N, names = FALSE), include = TRUE)

data_eGFR_scr_UACR_u30$follow_up_time = data_eGFR_scr_UACR_u30$all_infection_time
data_eGFR_scr_UACR_u30$event.type=data_eGFR_scr_UACR_u30$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scr_UACR_u30[which(data_eGFR_scr_UACR_u30$group==i),]
  temp = summary(lm(eGFR_scr_10~prs_eGFR_scr+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scr_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scr+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeana06 = summary_stat$mean_expos # mean exposure in each quantile
xmeana06
xmean1 = xmeana06-6  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

aa16 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=6, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cr", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
aa16

aa16$coefficients
aa16$p_tests
aa16$p_heterogeneity


fp_aa16=ifelse(aa16$p_tests[1,2]<0.001," < 0.001",ifelse(aa16$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa16$p_tests[1,2],3))))
Q_aa16=ifelse(aa16$p_tests[1,4]<0.001," < 0.001",ifelse(aa16$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa16$p_tests[1,4],3))))

fa116=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scr_UACR_u30$eGFR_scr, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scr_UACR_u30$eGFR_scr), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),aa16$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), aa16$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),aa16$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeana05)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cr]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cr"]," and All Infections (UACR > 30)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_aa16) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_aa16)),  size = 3)

fa116
infa116=data.frame(subgroup="UACR>30",exposure="eGFR_cr",n=dim(data_eGFR_scr_UACR_u30)[1],event_n=sum(data_eGFR_scr_UACR_u30$event.type),
                   P_Cochrans_Q=aa16$p_tests[1,4],p_non_linear=aa16$p_tests[1,2])


#data_eGFR_scr_UACR_l30
data_eGFR_scr_UACR_l30$scr_hat=lm(eGFR_scr_10~prs_eGFR_scr,data=data_eGFR_scr_UACR_l30)$fit
data_eGFR_scr_UACR_l30$scr_residual=data_eGFR_scr_UACR_l30$eGFR_scr_10-data_eGFR_scr_UACR_l30$scr_hat
N=10
data_eGFR_scr_UACR_l30$group=cut(data_eGFR_scr_UACR_l30$scr_residual, labels = F,quantile(data_eGFR_scr_UACR_l30$scr_residual, prob = 0:N / N, names = FALSE), include = TRUE)

data_eGFR_scr_UACR_l30$follow_up_time = data_eGFR_scr_UACR_l30$all_infection_time
data_eGFR_scr_UACR_l30$event.type=data_eGFR_scr_UACR_l30$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scr_UACR_l30[which(data_eGFR_scr_UACR_l30$group==i),]
  temp = summary(lm(eGFR_scr_10~prs_eGFR_scr+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scr_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scr+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeanb05 = summary_stat$mean_expos # mean exposure in each quantile
xmeanb05
xmean1 = xmeanb05-6  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

ba16 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=6, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cr", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
ba16

ba16$coefficients
ba16$p_tests
ba16$p_heterogeneity

fp_ba16=ifelse(ba16$p_tests[1,2]<0.001," < 0.001",ifelse(ba16$p_tests[1,2]>0.999," > 0.999",paste(" =",round(ba16$p_tests[1,2],3))))
Q_ba16=ifelse(ba16$p_tests[1,4]<0.001," < 0.001",ifelse(ba16$p_tests[1,4]>0.999," > 0.999",paste(" =",round(ba16$p_tests[1,4],3))))

fb116=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scr_UACR_l30$eGFR_scr, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scr_UACR_l30$eGFR_scr), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),ba16$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), ba16$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),ba16$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeanb05)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cr]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cr"]," and All Infections (UACR < 30)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_ba16) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_ba16)),  size = 3)

fb116
infb116=data.frame(subgroup="UACR<30",exposure="eGFR_cr",n=dim(data_eGFR_scr_UACR_l30)[1],event_n=sum(data_eGFR_scr_UACR_l30$event.type),
                   P_Cochrans_Q=ba16$p_tests[1,4],p_non_linear=ba16$p_tests[1,2])


#eGFR_scys

data=all_data_scys
data$UACR_2=1
data$UACR_2[which(data$UACR_cat==0)]=0


#logistic_regression
fit = glm(UACR_2~prs_eGFR_scys+BMI+age+sex+smoking+alcohol+MET+education+TDI+cancer,data = data,family = 'binomial')
data$prob = predict(fit,type="response")

data_eGFR_scys_UACR_u30=data[which(data$UACR_2==1),]
data_eGFR_scys_UACR_l30=data[which(data$UACR_2==0),]


#All_infection
#data_eGFR_scys_UACR_u30
data_eGFR_scys_UACR_u30$scys_hat=lm(eGFR_scys_10~prs_eGFR_scys,data=data_eGFR_scys_UACR_u30)$fit
data_eGFR_scys_UACR_u30$scys_residual=data_eGFR_scys_UACR_u30$eGFR_scys_10-data_eGFR_scys_UACR_u30$scys_hat
N=10
data_eGFR_scys_UACR_u30$group=cut(data_eGFR_scys_UACR_u30$scys_residual, labels = F,quantile(data_eGFR_scys_UACR_u30$scys_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_eGFR_scys_UACR_u30$follow_up_time = data_eGFR_scys_UACR_u30$all_infection_time
data_eGFR_scys_UACR_u30$event.type=data_eGFR_scys_UACR_u30$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scys_UACR_u30[which(data_eGFR_scys_UACR_u30$group==i),]
  temp = summary(lm(eGFR_scys_10~prs_eGFR_scys+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scys_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scys+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeana25 = summary_stat$mean_expos # mean exposure in each quantile
xmeana25
xmean1 = xmeana25-6  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

aa26 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=6, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cys", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
aa26

aa26$coefficients
aa26$p_tests
aa26$p_heterogeneity


fp_aa26=ifelse(aa26$p_tests[1,2]<0.001," < 0.001",ifelse(aa26$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa26$p_tests[1,2],3))))
Q_aa26=ifelse(aa26$p_tests[1,4]<0.001," < 0.001",ifelse(aa26$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa26$p_tests[1,4],3))))

fa216=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scys_UACR_u30$eGFR_scys, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scys_UACR_u30$eGFR_scys), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),aa26$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), aa26$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),aa26$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeana25)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cys]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cys"]," and All Infections (UACR > 30)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_aa26) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_aa26)),  size = 3)

fa216
infa216=data.frame(subgroup="UACR>30",exposure="eGFR_cys",n=dim(data_eGFR_scys_UACR_u30)[1],event_n=sum(data_eGFR_scys_UACR_u30$event.type),
                   P_Cochrans_Q=aa26$p_tests[1,4],p_non_linear=aa26$p_tests[1,2])

#data_eGFR_scys_UACR_l30
data_eGFR_scys_UACR_l30$scys_hat=lm(eGFR_scys_10~prs_eGFR_scys,data=data_eGFR_scys_UACR_l30)$fit
data_eGFR_scys_UACR_l30$scys_residual=data_eGFR_scys_UACR_l30$eGFR_scys_10-data_eGFR_scys_UACR_l30$scys_hat
N=10
data_eGFR_scys_UACR_l30$group=cut(data_eGFR_scys_UACR_l30$scys_residual, labels = F,quantile(data_eGFR_scys_UACR_l30$scys_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_eGFR_scys_UACR_l30$follow_up_time = data_eGFR_scys_UACR_l30$all_infection_time
data_eGFR_scys_UACR_l30$event.type=data_eGFR_scys_UACR_l30$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_eGFR_scys_UACR_l30[which(data_eGFR_scys_UACR_l30$group==i),]
  temp = summary(lm(eGFR_scys_10~prs_eGFR_scys+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$eGFR_scys_10)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_eGFR_scys+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeanb26 = summary_stat$mean_expos # mean exposure in each quantile
xmeanb26
xmean1 = xmeanb26-5  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

ba26 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=5, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="eGFR_cys", pref_y="Hazard ratio of All_infections", breaks=c(0.25,0.5,1,2,4,8))
ba26


ba26$coefficients
ba26$p_tests
ba26$p_heterogeneity

fp_ba26=ifelse(ba26$p_tests[1,2]<0.001," < 0.001",ifelse(ba26$p_tests[1,2]>0.999," > 0.999",paste(" =",round(ba26$p_tests[1,2],3))))
Q_ba26=ifelse(ba26$p_tests[1,4]<0.001," < 0.001",ifelse(ba26$p_tests[1,4]>0.999," > 0.999",paste(" =",round(ba26$p_tests[1,4],3))))

fb216=ggplot() + 
  geom_histogram(aes(x = data_eGFR_scys_UACR_l30$eGFR_scys, y = ((..count..) / sum(..count..)*10) ), data.frame(data_eGFR_scys_UACR_l30$eGFR_scys), binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1)+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x*10, y=log2(yest)),ba26$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x*10, y=log2(lci)), ba26$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x*10, y=log2(uci)),ba26$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeanb25)*10, y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3),  labels =c(0.25,0.5,1,2,4,8),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~eGFR[cys]))) + 
  scale_x_continuous(breaks=seq(130,50,by = -10)) + 
  coord_cartesian(xlim =c(  130,50)) + ggtitle(expression(paste(eGFR["cys"]," and All Infections (UACR < 30)")))+
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
  annotate("text", x = 65, y = 2.8, label = bquote(P['Cochran\'s Q test'] ~ .(Q_ba26) ), size = 3) + 
  annotate("text", x = 64.5, y = 2.5, label = bquote(P['nonlinearity'] ~.(fp_ba26)),  size = 3)

fb216
infb216=data.frame(subgroup="UACR<30",exposure="eGFR_cys",n=dim(data_eGFR_scys_UACR_l30)[1],event_n=sum(data_eGFR_scys_UACR_l30$event.type),
                   P_Cochrans_Q=ba26$p_tests[1,4],p_non_linear=ba26$p_tests[1,2])

#BUN
data=all_data_BUN
data$UACR_2=1
data$UACR_2[which(data$UACR_cat==0)]=0


#logistic_regression
fit = glm(UACR_2~prs_BUN+BMI+age+sex+smoking+alcohol+MET+education+TDI+cancer,data = data,family = 'binomial')
data$prob = predict(fit,type="response")

data_BUN_UACR_u30=data[which(data$UACR_2==1),]
data_BUN_UACR_l30=data[which(data$UACR_2==0),]


#All_infection
#data_BUN_UACR_u30
data_BUN_UACR_u30$BUN_hat=lm(BUN~prs_BUN,data=data_BUN_UACR_u30)$fit
data_BUN_UACR_u30$BUN_residual=data_BUN_UACR_u30$BUN-data_BUN_UACR_u30$BUN_hat
N=10
data_BUN_UACR_u30$group=cut(data_BUN_UACR_u30$BUN_residual, labels = F,quantile(data_BUN_UACR_u30$BUN_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_BUN_UACR_u30$follow_up_time = data_BUN_UACR_u30$all_infection_time
data_BUN_UACR_u30$event.type=data_BUN_UACR_u30$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_BUN_UACR_u30[which(data_BUN_UACR_u30$group==i),]
  temp = summary(lm(BUN~prs_BUN+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$BUN)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_BUN+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeana36 = summary_stat$mean_expos # mean exposure in each quantile
xmeana36
xmean1 = xmeana36-3  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

aa36 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=3, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="BUN", pref_y="Hazard ratio of all infections", breaks=c(0.25,0.5,1,2,4,8))
aa36

aa36$coefficients
aa36$p_tests
aa36$p_heterogeneity

fp_aa36=ifelse(aa36$p_tests[1,2]<0.001," < 0.001",ifelse(aa36$p_tests[1,2]>0.999," > 0.999",paste(" =",round(aa36$p_tests[1,2],3))))
Q_aa36=ifelse(aa36$p_tests[1,4]<0.001," < 0.001",ifelse(aa36$p_tests[1,4]>0.999," > 0.999",paste(" =",round(aa36$p_tests[1,4],3))))

fa316=ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), data_BUN_UACR_u30,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,12,1))+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x, y=log2(yest)),aa36$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x, y=log2(lci)), aa36$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x, y=log2(uci)),aa36$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeana36), y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 4),  labels =c(0.25,0.5,1,2,4,8,16),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~BUN))) + 
  scale_x_continuous(breaks=seq(0,12,by = 1)) + 
  coord_cartesian(xlim =c(0,12)) + ggtitle("BUN and All Infections (UACR > 30)")+
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
  annotate("text", x = 10, y = 3.7, label = bquote(P['Cochran\'s Q test'] ~ .(Q_aa36) ), size = 3) + 
  annotate("text", x = 10, y = 3.4, label = bquote(P['nonlinearity'] ~.(fp_aa36)),  size = 3)

fa316
infa316=data.frame(subgroup="UACR>30",exposure="BUN",n=dim(data_BUN_UACR_u30)[1],event_n=sum(data_BUN_UACR_u30$event.type),
                   P_Cochrans_Q=aa36$p_tests[1,4],p_non_linear=aa36$p_tests[1,2])

#data_BUN_UACR_l30
data_BUN_UACR_l30$BUN_hat=lm(BUN~prs_BUN,data=data_BUN_UACR_l30)$fit
data_BUN_UACR_l30$BUN_residual=data_BUN_UACR_l30$BUN-data_BUN_UACR_l30$BUN_hat
N=10
data_BUN_UACR_l30$group=cut(data_BUN_UACR_l30$BUN_residual, labels = F,quantile(data_BUN_UACR_l30$BUN_residual, prob = 0:N / N, names = FALSE), include = TRUE)


data_BUN_UACR_l30$follow_up_time = data_BUN_UACR_l30$all_infection_time
data_BUN_UACR_l30$event.type=data_BUN_UACR_l30$all_infection_event
summary_stat = data.frame(beta_expos=rep(NA,N),se_expos=rep(NA,N),beta_outcome=rep(NA,N),se_outcome=rep(NA,N),mean_expos=rep(NA,N))
for(i in 1:N){
  pheno = data_BUN_UACR_l30[which(data_BUN_UACR_l30$group==i),]
  temp = summary(lm(BUN~prs_BUN+age+sex+pc1+pc2,data=pheno, weights=1/prob))
  summary_stat$beta_expos[i] = temp$coefficients[2,1]
  summary_stat$se_expos[i] = temp$coefficients[2,2]
  summary_stat$mean_expos[i] = mean(pheno$BUN)
  temp2 = summary(coxph(Surv(follow_up_time, event.type) ~ prs_BUN+age+sex+pc1+pc2,data = pheno, weights = 1/prob))
  summary_stat$beta_outcome[i]=temp2$coefficients[1,1]
  summary_stat$se_outcome[i]=temp2$coefficients[1,3]
  
}
by   = summary_stat$beta_outcome  # genetic associations with outcome per quantile
byse = summary_stat$se_outcome    # standard errors of genetic associations with outcome
bx   = summary_stat$beta_expos  # genetic associations with exposure per quantile
bxse = summary_stat$se_expos    # standard errors of genetic associations with exposure
xmeanb36 = summary_stat$mean_expos # mean exposure in each quantile
xmeanb36
xmean1 = xmeanb36-3  # in our example we reduced BMI by 18 (see offset in function) - modelling is better when values start close to zero
# if your exposure values start near zero, this isn't needed

ba36 = frac_poly_summ_mr(by, bx, byse, bxse, xmean1, family="binomial",pd=0.05,ref=NA, d="both", offset=3, xlim_upper=NA, ylim_lower=NA, fig=TRUE, pref_x="BUN", pref_y="Hazard ratio of all infections", breaks=c(0.25,0.5,1,2,4,8))
ba36

ba36$coefficients
ba36$p_tests
ba36$p_heterogeneity

fp_ba36=ifelse(ba36$p_tests[1,2]<0.001," < 0.001",ifelse(ba36$p_tests[1,2]>0.999," > 0.999",paste(" =",round(ba36$p_tests[1,2],3))))
Q_ba36=ifelse(ba36$p_tests[1,4]<0.001," < 0.001",ifelse(ba36$p_tests[1,4]>0.999," > 0.999",paste(" =",round(ba36$p_tests[1,4],3))))

fb316=ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), data_BUN_UACR_l30,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,12,1))+  
  geom_hline(aes(yintercept=0), colour="grey") + 
  geom_line(aes(x=x, y=log2(yest)),ba36$plot.data, color=mycolors[9]) +
  #geom_ribbon(aes(x = x*10, y = log2(yest), ymin = log2(lci), ymax = log2(uci)),a$plot.data,  fill = mycolors[7], alpha = 0.1) +
  geom_line(aes(x=x, y=log2(lci)), ba36$plot.data, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = x, y=log2(uci)),ba36$plot.data, color=mycolors[8],linetype='dashed') +
  geom_point(aes(x = mean(xmeanb36), y = 0), colour="red", size=3) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 4),  labels =c(0.25,0.5,1,2,4,8,16),
                     sec.axis = sec_axis(~. / 10, name = expression(Density~of~BUN))) + 
  scale_x_continuous(breaks=seq(0,12,by = 1)) + 
  coord_cartesian(xlim =c(0,12)) + ggtitle("BUN and All Infections (UACR < 30)")+
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
  annotate("text", x = 10, y = 3.7, label = bquote(P['Cochran\'s Q test'] ~ .(Q_ba36) ), size = 3) +
  annotate("text", x = 10, y = 3.4, label = bquote(P['nonlinearity'] ~.(fp_ba36)),  size = 3)

fb316
infb316=data.frame(subgroup="UACR<30",exposure="BUN",n=dim(data_BUN_UACR_l30)[1],event_n=sum(data_BUN_UACR_l30$event.type),
                   P_Cochrans_Q=ba36$p_tests[1,4],p_non_linear=ba36$p_tests[1,2])
inf_sub5=rbind(infa116,infb116,infa216,infb216,infa316,infb316)

inf=rbind(inf_sub1,inf_sub2,inf_sub3,inf_sub4,inf_sub5)
write.csv(inf,"information_subgroup_all_infection.csv")
#################################################################################################
library(cowplot)


pdf("Nonlinear_MR_subgroup_all_infection.pdf",family = "Times",height = 16,width = 30)
plot_grid(fa112,fa212,fa312,fb112,fb212,fb312,
          fa113,fa213,fa313,fb113,fb213,fb313,
          fa114,fa214,fa314,fb114,fb214,fb314,
          fa115,fa215,fa315,fb115,fb215,fb315,
          fa116,fa216,fa316,fb116,fb216,fb316,
          nrow = 5,ncol =6)
dev.off()
