library(ggplot2)
library(rms)
library(RColorBrewer)
library(cowplot)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
all_data=read.table("all_data_imputation_cohort.txt", header=T)
#all_data1=read.table("all_data_imputation.txt", header=T)
#all_data_add=read.csv("medication_data_clean.csv")
#all_data <- merge(all_data1, all_data_add, by = "eid", all.x = TRUE)


all_data$education=factor(all_data$education,order=T)
all_data$smoking=factor(all_data$smoking,order=F)
all_data$alcohol=factor(all_data$alcohol,order=F)
all_data$UACR_cat=factor(all_data$UACR_cat,order=T)
all_data$eGFR_scr_cat=factor(all_data$eGFR_scr_cat,order=F)
all_data$eGFR_scys_cat=factor(all_data$eGFR_scys_cat,order=F)
all_data$BUN_cat=factor(all_data$BUN_cat,order=F)
all_data$medicine_g_i=factor(all_data$medicine_g_i,order=F)

#all_infection

all_data$event.type <- all_data$all_infection_event
all_data$follow_up_time<-all_data$all_infection_time

lambda=1
mycolors<-brewer.pal(9, "Blues")
covar = c("age","sex","education","TDI","smoking","alcohol","CVD","chronic_pulmonary_disease",
            "cancer","diabetes","hypertension","UACR_cat","BMI","medicine_g_i")

#d
####################################################
##eGFR_scr
need_term <- c(covar,"eGFR_scr")
data_c2_nona <- all_data[complete.cases(all_data[,c(need_term)]),c(need_term,"event.type","follow_up_time")]
dd <- datadist(data_c2_nona) #为后续程序设定数据环境
options(datadist='dd')

sm1.1 = cph(Surv(follow_up_time, event.type) ~ age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
              cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i+
              rcs(eGFR_scr,quantile(eGFR_scr, c( .05,0.275, .5,0.725,  .95))),
            #rcs(eGFR_10,c( 4.5, 6.0, 7.5, 9.0, 10.5 )),
            data = data_c2_nona)
fit=update(sm1.1)
#
HR<-Predict(fit,eGFR_scr,fun=exp,ref.zero = TRUE)

a <- as.data.frame(anova(fit))
a
p_overall <- ifelse(a$P[which(rownames(a)=="eGFR_scr")]<0.001,"<0.001",round(a$P[which(rownames(a)=="eGFR_scr")],digits = 3))
p_overall  #<0.001
p_nonlinear <- ifelse(a$P[which(rownames(a)=="X.Nonlinear")]<0.001,"<0.001",round(a$P[which(rownames(a)=="X.Nonlinear")],digits = 3))
p_nonlinear  #<0.001

HR$yhat = exp(log(HR$yhat)/lambda)
HR$lower = exp(log(HR$lower)/lambda)
HR$upper = exp(log(HR$upper)/lambda)

p_major_eGFRscr1 <-ggplot() + 
  geom_histogram(aes(x = eGFR_scr, y = (..count..) / sum(..count..) * 10), data_c2_nona,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(5,155,5)) +  
  geom_line(aes(x = eGFR_scr, y = yhat), HR, linetype = 1, col = mycolors[9], size = 0.8) + 
  geom_line(aes(x = eGFR_scr, y=lower), HR, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = eGFR_scr, y=upper),HR, color=mycolors[8],linetype='dashed') +
  geom_hline(yintercept = 1, linetype = 8, size = 0.5,color="grey40") + 
  # geom_point(aes(x = PM2.5_knots, y = PM2.5_yhat)) + 
  scale_x_reverse(breaks = seq(10,160,10))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,6,1),  
                     sec.axis = sec_axis(~. / 10, name = expression(paste("Density of ",eGFR["cr"])))) + 
  coord_cartesian(xlim = c(160,10), ylim = c(0,6)) +
  labs(x = expression(paste(eGFR["cr"]," (","ml/min/1.73 ",m^2,")")), y = "Hazard Ratio",title = expression(paste(eGFR["cr"]," and All Infections"))) + 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text( size = 10,vjust = 0.5), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")),    
        plot.tag = element_text( size = 14, face = "bold"),
        plot.title = element_text(size = 15, face = "bold")) + 
  annotate("text", x = 40, y = 5.3, label = paste("P(overall)", p_overall), size = 3.5)+
  annotate("text", x = 40, y = 5.7, label = paste("P(nonlinearity)", p_nonlinear),  size = 3.5)
p_major_eGFRscr1


#d
####################################################
##eGFR_scys
need_term <- c(covar,"eGFR_scys")
data_c2_nona <- all_data[complete.cases(all_data[,c(need_term)]),c(need_term,"event.type","follow_up_time")]
dd <- datadist(data_c2_nona) #为后续程序设定数据环境
options(datadist='dd')

sm1.2 = cph(Surv(follow_up_time, event.type) ~ age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
              cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i+
              rcs(eGFR_scys,quantile(eGFR_scys, c( .05,0.275, .5,0.725,  .95))),
            #rcs(eGFR_10,c( 4.5, 6.0, 7.5, 9.0, 10.5 )),
            data = data_c2_nona)
fit=update(sm1.2)
#
HR<-Predict(fit,eGFR_scys,fun=exp,ref.zero = TRUE)

a <- as.data.frame(anova(fit))
a
p_overall <- ifelse(a$P[which(rownames(a)=="eGFR_scys")]<0.001,"<0.001",round(a$P[which(rownames(a)=="eGFR_scys")],digits = 3))
p_overall  #<0.001
p_nonlinear <- ifelse(a$P[which(rownames(a)=="X.Nonlinear")]<0.001,"<0.001",round(a$P[which(rownames(a)=="X.Nonlinear")],digits = 3))
p_nonlinear  #<0.001

HR$yhat = exp(log(HR$yhat)/lambda)
HR$lower = exp(log(HR$lower)/lambda)
HR$upper = exp(log(HR$upper)/lambda)

p_major_eGFRscys1 <-ggplot() + 
  geom_histogram(aes(x = eGFR_scys, y = (..count..) / sum(..count..) * 10), data_c2_nona,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(5,155,5)) +  
  geom_line(aes(x = eGFR_scys, y = yhat), HR, linetype = 1, col = mycolors[9], size = 0.8) + 
  geom_line(aes(x = eGFR_scys, y=lower), HR, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = eGFR_scys, y=upper),HR, color=mycolors[8],linetype='dashed') +
  geom_hline(yintercept = 1, linetype = 8, size = 0.5,color="grey40") + 
  # geom_point(aes(x = PM2.5_knots, y = PM2.5_yhat)) + 
  scale_x_reverse(breaks = seq(10,160,10))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,6,1),  
                     sec.axis = sec_axis(~. / 10, name = expression(paste("Density of ",eGFR["cys"])))) + 
  coord_cartesian(xlim = c(160,10), ylim = c(0,6)) +
  labs(x = expression(paste(eGFR["cys"]," (","ml/min/1.73 ",m^2,")")), y = "Hazard Ratio",title = expression(paste(eGFR["cys"]," and All Infections"))) + 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text( size = 10,vjust = 0.5), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")),    
        plot.tag = element_text( size = 14, face = "bold"),
        plot.title = element_text(size = 15, face = "bold")) + 
  annotate("text", x = 120, y = 5.3, label = paste("P(overall)", p_overall), size = 3.5)+
  annotate("text", x = 120, y = 5.7, label = paste("P(nonlinearity)", p_nonlinear),  size = 3.5)
p_major_eGFRscys1


####################################################
##BUN
need_term <- c(covar,"BUN")
data_c2_nona <- all_data[complete.cases(all_data[,c(need_term)]),c(need_term,"event.type","follow_up_time")]
dd <- datadist(data_c2_nona) #为后续程序设定数据环境
options(datadist='dd')

sm1.3 = cph(Surv(follow_up_time, event.type) ~ age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
              cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i+
              rcs(BUN,quantile(BUN, c( .05,0.275, .5,0.725,  .95))),
            #rcs(eGFR_10,c( 4.5, 6.0, 7.5, 9.0, 10.5 )),
            data = data_c2_nona)
fit=update(sm1.3)
#
HR<-Predict(fit,BUN,fun=exp,ref.zero = TRUE)

a <- as.data.frame(anova(fit))
a
p_overall <- ifelse(a$P[which(rownames(a)=="BUN")]<0.001,"<0.001",round(a$P[which(rownames(a)=="BUN")],digits = 3))
p_overall  #<0.001
p_nonlinear <- ifelse(a$P[which(rownames(a)=="X.Nonlinear")]<0.001,"<0.001",round(a$P[which(rownames(a)=="X.Nonlinear")],digits = 3))
p_nonlinear  #<0.001

HR$yhat = exp(log(HR$yhat)/lambda)
HR$lower = exp(log(HR$lower)/lambda)
HR$upper = exp(log(HR$upper)/lambda)

p_major_BUN1 <-ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), data_c2_nona,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,35,1)) +  
  geom_line(aes(x = BUN, y = yhat), HR, linetype = 1, col = mycolors[9], size = 0.8) + 
  geom_line(aes(x = BUN, y=lower), HR, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = BUN, y=upper),HR, color=mycolors[8],linetype='dashed') +
  geom_hline(yintercept = 1, linetype = 8, size = 0.5,color="grey40") + 
  # geom_point(aes(x = PM2.5_knots, y = PM2.5_yhat)) + 
  scale_x_continuous(breaks = seq(0,35,2))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,6,1),  
                     sec.axis = sec_axis(~. / 10, name = expression(paste("Density of ","BUN")))) + 
  coord_cartesian(xlim = c(0,35), ylim = c(0,6)) +
  labs(x = "BUN", y = "Hazard Ratio",title = "BUN and All Infections") + 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text( size = 10,vjust = 0.5), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")),    
        plot.tag = element_text( size = 14, face = "bold"),
        plot.title = element_text(size = 15, face = "bold")) + 
  annotate("text", x = 8, y = 5.3, label = paste("P(overall)", p_overall), size = 3.5)+
  annotate("text", x = 8, y = 5.7, label = paste("P(nonlinearity)", p_nonlinear),  size = 3.5)
p_major_BUN1




#pneumonia

all_data$event.type <- all_data$pneumonia_event
all_data$follow_up_time<-all_data$pneumonia_time

lambda=1
mycolors<-brewer.pal(9, "Blues")
covar = c("age","sex","education","TDI","smoking","alcohol","CVD","chronic_pulmonary_disease",
          "cancer","diabetes","hypertension","UACR_cat","BMI","medicine_g_i")

#d
####################################################
##eGFR_scr
need_term <- c(covar,"eGFR_scr")
data_c2_nona <- all_data[complete.cases(all_data[,c(need_term)]),c(need_term,"event.type","follow_up_time")]
dd <- datadist(data_c2_nona) #为后续程序设定数据环境
options(datadist='dd')

sm1.1 = cph(Surv(follow_up_time, event.type) ~ age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
              cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i+
              rcs(eGFR_scr,quantile(eGFR_scr, c( .05,0.275, .5,0.725,  .95))),
            #rcs(eGFR_10,c( 4.5, 6.0, 7.5, 9.0, 10.5 )),
            data = data_c2_nona)
fit=update(sm1.1)
#
HR<-Predict(fit,eGFR_scr,fun=exp,ref.zero = TRUE)

a <- as.data.frame(anova(fit))
a
p_overall <- ifelse(a$P[which(rownames(a)=="eGFR_scr")]<0.001,"<0.001",round(a$P[which(rownames(a)=="eGFR_scr")],digits = 3))
p_overall  #<0.001
p_nonlinear <- ifelse(a$P[which(rownames(a)=="X.Nonlinear")]<0.001,"<0.001",round(a$P[which(rownames(a)=="X.Nonlinear")],digits = 3))
p_nonlinear  #<0.001

HR$yhat = exp(log(HR$yhat)/lambda)
HR$lower = exp(log(HR$lower)/lambda)
HR$upper = exp(log(HR$upper)/lambda)

p_major_eGFRscr2 <-ggplot() + 
  geom_histogram(aes(x = eGFR_scr, y = (..count..) / sum(..count..) * 10), data_c2_nona,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(5,155,5)) +  
  geom_line(aes(x = eGFR_scr, y = yhat), HR, linetype = 1, col = mycolors[9], size = 0.8) + 
  geom_line(aes(x = eGFR_scr, y=lower), HR, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = eGFR_scr, y=upper),HR, color=mycolors[8],linetype='dashed') +
  geom_hline(yintercept = 1, linetype = 8, size = 0.5,color="grey40") + 
  # geom_point(aes(x = PM2.5_knots, y = PM2.5_yhat)) + 
  scale_x_reverse(breaks = seq(10,160,10))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,6,1),  
                     sec.axis = sec_axis(~. / 10, name = expression(paste("Density of ",eGFR["cr"])))) + 
  coord_cartesian(xlim = c(160,10), ylim = c(0,6)) +
  labs(x = expression(paste(eGFR["cr"]," (","ml/min/1.73 ",m^2,")")), y = "Hazard Ratio",title = expression(paste(eGFR["cr"]," and Pneumonia"))) + 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text( size = 10,vjust = 0.5), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")),    
        plot.tag = element_text( size = 14, face = "bold"),
        plot.title = element_text(size = 15, face = "bold")) + 
  annotate("text", x = 60, y = 5.3, label = paste("P(overall)", p_overall), size = 3.5)+
  annotate("text", x = 60, y = 5.7, label = paste("P(nonlinearity)", p_nonlinear),  size = 3.5)
p_major_eGFRscr2


#d
####################################################
##eGFR_scys
need_term <- c(covar,"eGFR_scys")
data_c2_nona <- all_data[complete.cases(all_data[,c(need_term)]),c(need_term,"event.type","follow_up_time")]
dd <- datadist(data_c2_nona) #为后续程序设定数据环境
options(datadist='dd')

sm1.2 = cph(Surv(follow_up_time, event.type) ~ age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
              cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i+
              rcs(eGFR_scys,quantile(eGFR_scys, c( .05,0.275, .5,0.725,  .95))),
            #rcs(eGFR_10,c( 4.5, 6.0, 7.5, 9.0, 10.5 )),
            data = data_c2_nona)
fit=update(sm1.2)
#
HR<-Predict(fit,eGFR_scys,fun=exp,ref.zero = TRUE)

a <- as.data.frame(anova(fit))
a
p_overall <- ifelse(a$P[which(rownames(a)=="eGFR_scys")]<0.001,"<0.001",round(a$P[which(rownames(a)=="eGFR_scys")],digits = 3))
p_overall  #<0.001
p_nonlinear <- ifelse(a$P[which(rownames(a)=="X.Nonlinear")]<0.001,"<0.001",round(a$P[which(rownames(a)=="X.Nonlinear")],digits = 3))
p_nonlinear  #<0.001

HR$yhat = exp(log(HR$yhat)/lambda)
HR$lower = exp(log(HR$lower)/lambda)
HR$upper = exp(log(HR$upper)/lambda)

p_major_eGFRscys2 <-ggplot() + 
  geom_histogram(aes(x = eGFR_scys, y = (..count..) / sum(..count..) * 10), data_c2_nona,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(5,155,5)) +  
  geom_line(aes(x = eGFR_scys, y = yhat), HR, linetype = 1, col = mycolors[9], size = 0.8) + 
  geom_line(aes(x = eGFR_scys, y=lower), HR, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = eGFR_scys, y=upper),HR, color=mycolors[8],linetype='dashed') +
  geom_hline(yintercept = 1, linetype = 8, size = 0.5,color="grey40") + 
  # geom_point(aes(x = PM2.5_knots, y = PM2.5_yhat)) + 
  scale_x_reverse(breaks = seq(10,160,10))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,6,1),  
                     sec.axis = sec_axis(~. / 10, name = expression(paste("Density of ",eGFR["cys"])))) + 
  coord_cartesian(xlim = c(160,10), ylim = c(0,6)) +
  labs(x = expression(paste(eGFR["cys"]," (","ml/min/1.73 ",m^2,")")), y = "Hazard Ratio",title = expression(paste(eGFR["cys"]," and Pneumonia"))) + 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text( size = 10,vjust = 0.5), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")),    
        plot.tag = element_text( size = 14, face = "bold"),
        plot.title = element_text(size = 15, face = "bold")) + 
  annotate("text", x = 120, y = 5.3, label = paste("P(overall)", p_overall), size = 3.5)+
  annotate("text", x = 120, y = 5.7, label = paste("P(nonlinearity)", p_nonlinear),  size = 3.5)
p_major_eGFRscys2


####################################################
##BUN
need_term <- c(covar,"BUN")
data_c2_nona <- all_data[complete.cases(all_data[,c(need_term)]),c(need_term,"event.type","follow_up_time")]
dd <- datadist(data_c2_nona) #为后续程序设定数据环境
options(datadist='dd')

sm1.3 = cph(Surv(follow_up_time, event.type) ~ age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
              cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i+
              rcs(BUN,quantile(BUN, c( .05,0.275, .5,0.725,  .95))),
            #rcs(eGFR_10,c( 4.5, 6.0, 7.5, 9.0, 10.5 )),
            data = data_c2_nona)
fit=update(sm1.3)
#
HR<-Predict(fit,BUN,fun=exp,ref.zero = TRUE)

a <- as.data.frame(anova(fit))
a
p_overall <- ifelse(a$P[which(rownames(a)=="BUN")]<0.001,"<0.001",round(a$P[which(rownames(a)=="BUN")],digits = 3))
p_overall  #<0.001
p_nonlinear <- ifelse(a$P[which(rownames(a)=="X.Nonlinear")]<0.001,"<0.001",round(a$P[which(rownames(a)=="X.Nonlinear")],digits = 3))
p_nonlinear  #<0.001

HR$yhat = exp(log(HR$yhat)/lambda)
HR$lower = exp(log(HR$lower)/lambda)
HR$upper = exp(log(HR$upper)/lambda)

p_major_BUN2 <-ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), data_c2_nona,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,35,1)) +  
  geom_line(aes(x = BUN, y = yhat), HR, linetype = 1, col = mycolors[9], size = 0.8) + 
  geom_line(aes(x = BUN, y=lower), HR, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = BUN, y=upper),HR, color=mycolors[8],linetype='dashed') +
  geom_hline(yintercept = 1, linetype = 8, size = 0.5,color="grey40") + 
  # geom_point(aes(x = PM2.5_knots, y = PM2.5_yhat)) + 
  scale_x_continuous(breaks = seq(0,35,2))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,6,1),  
                     sec.axis = sec_axis(~. / 10, name = expression(paste("Density of ","BUN")))) + 
  coord_cartesian(xlim = c(0,35), ylim = c(0,6)) +
  labs(x = "BUN", y = "Hazard Ratio",title = "BUN and Pneumonia") + 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text( size = 10,vjust = 0.5), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")),    
        plot.tag = element_text( size = 14, face = "bold"),
        plot.title = element_text(size = 15, face = "bold")) + 
  annotate("text", x = 8, y = 5.3, label = paste("P(overall)", p_overall), size = 3.5)+
  annotate("text", x = 8, y = 5.7, label = paste("P(nonlinearity)", p_nonlinear),  size = 3.5)
p_major_BUN2



#Sepsis

all_data$event.type <- all_data$sepsis_event
all_data$follow_up_time<-all_data$sepsis_time

lambda=1
mycolors<-brewer.pal(9, "Blues")
covar = c("age","sex","education","TDI","smoking","alcohol","CVD","chronic_pulmonary_disease",
          "cancer","diabetes","hypertension","UACR_cat","BMI","medicine_g_i")

#d
####################################################
##eGFR_scr
need_term <- c(covar,"eGFR_scr")
data_c2_nona <- all_data[complete.cases(all_data[,c(need_term)]),c(need_term,"event.type","follow_up_time")]
dd <- datadist(data_c2_nona) #为后续程序设定数据环境
options(datadist='dd')

sm1.1 = cph(Surv(follow_up_time, event.type) ~ age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
              cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i+
              rcs(eGFR_scr,quantile(eGFR_scr, c( .05,0.275, .5,0.725,  .95))),
            #rcs(eGFR_10,c( 4.5, 6.0, 7.5, 9.0, 10.5 )),
            data = data_c2_nona)
fit=update(sm1.1)
#
HR<-Predict(fit,eGFR_scr,fun=exp,ref.zero = TRUE)

a <- as.data.frame(anova(fit))
a
p_overall <- ifelse(a$P[which(rownames(a)=="eGFR_scr")]<0.001,"<0.001",round(a$P[which(rownames(a)=="eGFR_scr")],digits = 3))
p_overall  #<0.001
p_nonlinear <- ifelse(a$P[which(rownames(a)=="X.Nonlinear")]<0.001,"<0.001",round(a$P[which(rownames(a)=="X.Nonlinear")],digits = 3))
p_nonlinear  #<0.001

HR$yhat = exp(log(HR$yhat)/lambda)
HR$lower = exp(log(HR$lower)/lambda)
HR$upper = exp(log(HR$upper)/lambda)

p_major_eGFRscr3 <-ggplot() + 
  geom_histogram(aes(x = eGFR_scr, y = (..count..) / sum(..count..) * 10), data_c2_nona,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(5,155,5)) +  
  geom_line(aes(x = eGFR_scr, y = yhat), HR, linetype = 1, col = mycolors[9], size = 0.8) + 
  geom_line(aes(x = eGFR_scr, y=lower), HR, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = eGFR_scr, y=upper),HR, color=mycolors[8],linetype='dashed') + 
  geom_hline(yintercept = 1, linetype = 8, size = 0.5,color="grey40") + 
  # geom_point(aes(x = PM2.5_knots, y = PM2.5_yhat)) + 
  scale_x_reverse(breaks = seq(10,160,10))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,6,1),  
                     sec.axis = sec_axis(~. / 10, name = expression(paste("Density of ",eGFR["cr"])))) + 
  coord_cartesian(xlim = c(160,10), ylim = c(0,6)) +
  labs(x = expression(paste(eGFR["cr"]," (","ml/min/1.73 ",m^2,")")), y = "Hazard Ratio",title = expression(paste(eGFR["cr"]," and Sepsis"))) + 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text( size = 10,vjust = 0.5), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")),    
        plot.tag = element_text( size = 14, face = "bold"),
        plot.title = element_text(size = 15, face = "bold")) + 
  annotate("text", x = 55, y = 5.3, label = paste("P(overall)", p_overall), size = 3.5)+
  annotate("text", x = 55, y = 5.7, label = paste("P(nonlinearity)", p_nonlinear),  size = 3.5)
p_major_eGFRscr3


#d
####################################################
##eGFR_scys
need_term <- c(covar,"eGFR_scys")
data_c2_nona <- all_data[complete.cases(all_data[,c(need_term)]),c(need_term,"event.type","follow_up_time")]
dd <- datadist(data_c2_nona) #为后续程序设定数据环境
options(datadist='dd')

sm1.2 = cph(Surv(follow_up_time, event.type) ~ age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
              cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i+
              rcs(eGFR_scys,quantile(eGFR_scys, c( .05,0.275, .5,0.725,  .95))),
            #rcs(eGFR_10,c( 4.5, 6.0, 7.5, 9.0, 10.5 )),
            data = data_c2_nona)
fit=update(sm1.2)
#
HR<-Predict(fit,eGFR_scys,fun=exp,ref.zero = TRUE)

a <- as.data.frame(anova(fit))
a
p_overall <- ifelse(a$P[which(rownames(a)=="eGFR_scys")]<0.001,"<0.001",round(a$P[which(rownames(a)=="eGFR_scys")],digits = 3))
p_overall  #<0.001
p_nonlinear <- ifelse(a$P[which(rownames(a)=="X.Nonlinear")]<0.001,"<0.001",round(a$P[which(rownames(a)=="X.Nonlinear")],digits = 3))
p_nonlinear  #<0.001

HR$yhat = exp(log(HR$yhat)/lambda)
HR$lower = exp(log(HR$lower)/lambda)
HR$upper = exp(log(HR$upper)/lambda)

p_major_eGFRscys3 <-ggplot() + 
  geom_histogram(aes(x = eGFR_scys, y = (..count..) / sum(..count..) * 10), data_c2_nona,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(5,155,5)) +  
  geom_line(aes(x = eGFR_scys, y = yhat), HR, linetype = 1, col = mycolors[9], size = 0.8) + 
  geom_line(aes(x = eGFR_scys, y=lower), HR, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = eGFR_scys, y=upper),HR, color=mycolors[8],linetype='dashed') + 
  geom_hline(yintercept = 1, linetype = 8, size = 0.5,color="grey40") + 
  # geom_point(aes(x = PM2.5_knots, y = PM2.5_yhat)) + 
  scale_x_reverse(breaks = seq(10,160,10))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,6,1),  
                     sec.axis = sec_axis(~. / 10, name = expression(paste("Density of ",eGFR["cys"])))) + 
  coord_cartesian(xlim = c(160,10), ylim = c(0,6)) +
  labs(x = expression(paste(eGFR["cys"]," (","ml/min/1.73 ",m^2,")")), y = "Hazard Ratio",title = expression(paste(eGFR["cys"]," and Sepsis"))) + 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text( size = 10,vjust = 0.5), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")),    
        plot.tag = element_text( size = 14, face = "bold"),
        plot.title = element_text(size = 15, face = "bold")) + 
  annotate("text", x = 120, y = 5.3, label = paste("P(overall)", p_overall), size = 3.5)+
  annotate("text", x = 120, y = 5.7, label = paste("P(nonlinearity)", p_nonlinear),  size = 3.5)
p_major_eGFRscys3


####################################################
##BUN
need_term <- c(covar,"BUN")
data_c2_nona <- all_data[complete.cases(all_data[,c(need_term)]),c(need_term,"event.type","follow_up_time")]
dd <- datadist(data_c2_nona) #为后续程序设定数据环境
options(datadist='dd')

sm1.3 = cph(Surv(follow_up_time, event.type) ~ age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
              cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i+
              rcs(BUN,quantile(BUN, c( .05,0.275, .5,0.725,  .95))),
            #rcs(eGFR_10,c( 4.5, 6.0, 7.5, 9.0, 10.5 )),
            data = data_c2_nona)
fit=update(sm1.3)
#
HR<-Predict(fit,BUN,fun=exp,ref.zero = TRUE)

a <- as.data.frame(anova(fit))
a
p_overall <- ifelse(a$P[which(rownames(a)=="BUN")]<0.001,"<0.001",round(a$P[which(rownames(a)=="BUN")],digits = 3))
p_overall  #<0.001
p_nonlinear <- ifelse(a$P[which(rownames(a)=="X.Nonlinear")]<0.001,"<0.001",round(a$P[which(rownames(a)=="X.Nonlinear")],digits = 3))
p_nonlinear  #<0.001

HR$yhat = exp(log(HR$yhat)/lambda)
HR$lower = exp(log(HR$lower)/lambda)
HR$upper = exp(log(HR$upper)/lambda)

p_major_BUN3 <-ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), data_c2_nona,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,35,1)) +  
  geom_line(aes(x = BUN, y = yhat), HR, linetype = 1, col = mycolors[9], size = 0.8) + 
  geom_line(aes(x = BUN, y=lower), HR, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = BUN, y=upper),HR, color=mycolors[8],linetype='dashed') +
  geom_hline(yintercept = 1, linetype = 8, size = 0.5,color="grey40") + 
  # geom_point(aes(x = PM2.5_knots, y = PM2.5_yhat)) + 
  scale_x_continuous(breaks = seq(0,35,2))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,6,1),  
                     sec.axis = sec_axis(~. / 10, name = expression(paste("Density of ","BUN")))) + 
  coord_cartesian(xlim = c(0,35), ylim = c(0,6)) +
  labs(x = "BUN", y = "Hazard Ratio",title = "BUN and Sepsis") + 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text( size = 10,vjust = 0.5), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")),    
        plot.tag = element_text( size = 14, face = "bold"),
        plot.title = element_text(size = 15, face = "bold")) + 
  annotate("text", x = 8, y = 5.3, label = paste("P(overall)", p_overall), size = 3.5)+
  annotate("text", x = 8, y = 5.7, label = paste("P(nonlinearity)", p_nonlinear),  size = 3.5)
p_major_BUN3



#skin_soft_tissue

all_data$event.type <- all_data$skin_soft_tissue_event
all_data$follow_up_time<-all_data$skin_soft_tissue_time

lambda=1
mycolors<-brewer.pal(9, "Blues")
covar = c("age","sex","education","TDI","smoking","alcohol","CVD","chronic_pulmonary_disease",
          "cancer","diabetes","hypertension","UACR_cat","BMI","medicine_g_i")

#d
####################################################
##eGFR_scr
need_term <- c(covar,"eGFR_scr")
data_c2_nona <- all_data[complete.cases(all_data[,c(need_term)]),c(need_term,"event.type","follow_up_time")]
dd <- datadist(data_c2_nona) #为后续程序设定数据环境
options(datadist='dd')

sm1.1 = cph(Surv(follow_up_time, event.type) ~ age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
              cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i+
              rcs(eGFR_scr,quantile(eGFR_scr, c( .05,0.275, .5,0.725,  .95))),
            #rcs(eGFR_10,c( 4.5, 6.0, 7.5, 9.0, 10.5 )),
            data = data_c2_nona)
fit=update(sm1.1)
#
HR<-Predict(fit,eGFR_scr,fun=exp,ref.zero = TRUE)

a <- as.data.frame(anova(fit))
a
p_overall <- ifelse(a$P[which(rownames(a)=="eGFR_scr")]<0.001,"<0.001",round(a$P[which(rownames(a)=="eGFR_scr")],digits = 3))
p_overall  #<0.001
p_nonlinear <- ifelse(a$P[which(rownames(a)=="X.Nonlinear")]<0.001,"<0.001",round(a$P[which(rownames(a)=="X.Nonlinear")],digits = 3))
p_nonlinear  #<0.001

HR$yhat = exp(log(HR$yhat)/lambda)
HR$lower = exp(log(HR$lower)/lambda)
HR$upper = exp(log(HR$upper)/lambda)

p_major_eGFRscr4 <-ggplot() + 
  geom_histogram(aes(x = eGFR_scr, y = (..count..) / sum(..count..) * 10), data_c2_nona,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(5,155,5)) +  
  geom_line(aes(x = eGFR_scr, y = yhat), HR, linetype = 1, col = mycolors[9], size = 0.8) + 
  geom_line(aes(x = eGFR_scr, y=lower), HR, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = eGFR_scr, y=upper),HR, color=mycolors[8],linetype='dashed') +
  geom_hline(yintercept = 1, linetype = 8, size = 0.5,color="grey40") + 
  # geom_point(aes(x = PM2.5_knots, y = PM2.5_yhat)) + 
  scale_x_reverse(breaks = seq(10,160,10))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,6,1),  
                     sec.axis = sec_axis(~. / 10, name = expression(paste("Density of ",eGFR["cr"])))) + 
  coord_cartesian(xlim = c(160,10), ylim = c(0,6)) +
  labs(x = expression(paste(eGFR["cr"]," (","ml/min/1.73 ",m^2,")")), y = "Hazard Ratio",title = expression(paste(eGFR["cr"]," and Skin and Soft Tissue Infections"))) + 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text( size = 10,vjust = 0.5), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")),    
        plot.tag = element_text( size = 14, face = "bold"),
        plot.title = element_text(size = 13, face = "bold")) + 
  annotate("text", x = 40, y = 5.3, label = paste("P(overall)", p_overall), size = 3.5)+
  annotate("text", x = 40, y = 5.7, label = paste("P(nonlinearity)", p_nonlinear),  size = 3.5)
p_major_eGFRscr4


#d
####################################################
##eGFR_scys
need_term <- c(covar,"eGFR_scys")
data_c2_nona <- all_data[complete.cases(all_data[,c(need_term)]),c(need_term,"event.type","follow_up_time")]
dd <- datadist(data_c2_nona) #为后续程序设定数据环境
options(datadist='dd')

sm1.2 = cph(Surv(follow_up_time, event.type) ~ age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
              cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i+
              rcs(eGFR_scys,quantile(eGFR_scys, c( .05,0.275, .5,0.725,  .95))),
            #rcs(eGFR_10,c( 4.5, 6.0, 7.5, 9.0, 10.5 )),
            data = data_c2_nona)
fit=update(sm1.2)
#
HR<-Predict(fit,eGFR_scys,fun=exp,ref.zero = TRUE)

a <- as.data.frame(anova(fit))
a
p_overall <- ifelse(a$P[which(rownames(a)=="eGFR_scys")]<0.001,"<0.001",round(a$P[which(rownames(a)=="eGFR_scys")],digits = 3))
p_overall  #<0.001
p_nonlinear <- ifelse(a$P[which(rownames(a)=="X.Nonlinear")]<0.001,"<0.001",round(a$P[which(rownames(a)=="X.Nonlinear")],digits = 3))
p_nonlinear  #<0.001

HR$yhat = exp(log(HR$yhat)/lambda)
HR$lower = exp(log(HR$lower)/lambda)
HR$upper = exp(log(HR$upper)/lambda)

p_major_eGFRscys4 <-ggplot() + 
  geom_histogram(aes(x = eGFR_scys, y = (..count..) / sum(..count..) * 10), data_c2_nona,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(5,155,5)) +  
  geom_line(aes(x = eGFR_scys, y = yhat), HR, linetype = 1, col = mycolors[9], size = 0.8) + 
  geom_line(aes(x = eGFR_scys, y=lower), HR, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = eGFR_scys, y=upper),HR, color=mycolors[8],linetype='dashed') +
  geom_hline(yintercept = 1, linetype = 8, size = 0.5,color="grey40") + 
  # geom_point(aes(x = PM2.5_knots, y = PM2.5_yhat)) + 
  scale_x_reverse(breaks = seq(10,160,10))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,6,1),  
                     sec.axis = sec_axis(~. / 10, name = expression(paste("Density of ",eGFR["cys"])))) + 
  coord_cartesian(xlim = c(160,10), ylim = c(0,6)) +
  labs(x = expression(paste(eGFR["cys"]," (","ml/min/1.73 ",m^2,")")), y = "Hazard Ratio",title = expression(paste(eGFR["cys"]," and Skin and Soft Tissue Infections"))) + 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text( size = 10,vjust = 0.5), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")),    
        plot.tag = element_text( size = 14, face = "bold"),
        plot.title = element_text(size = 13, face = "bold")) + 
  annotate("text", x = 120, y = 5.3, label = paste("P(overall)", p_overall), size = 3.5)+
  annotate("text", x = 120, y = 5.7, label = paste("P(nonlinearity)", p_nonlinear),  size = 3.5)
p_major_eGFRscys4


####################################################
##BUN
need_term <- c(covar,"BUN")
data_c2_nona <- all_data[complete.cases(all_data[,c(need_term)]),c(need_term,"event.type","follow_up_time")]
dd <- datadist(data_c2_nona) #为后续程序设定数据环境
options(datadist='dd')

sm1.3 = cph(Surv(follow_up_time, event.type) ~ age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
              cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i+
              rcs(BUN,quantile(BUN, c( .05,0.275, .5,0.725,  .95))),
            #rcs(eGFR_10,c( 4.5, 6.0, 7.5, 9.0, 10.5 )),
            data = data_c2_nona)
fit=update(sm1.3)
#
HR<-Predict(fit,BUN,fun=exp,ref.zero = TRUE)

a <- as.data.frame(anova(fit))
a
p_overall <- ifelse(a$P[which(rownames(a)=="BUN")]<0.001,"<0.001",round(a$P[which(rownames(a)=="BUN")],digits = 3))
p_overall  #<0.001
p_nonlinear <- ifelse(a$P[which(rownames(a)=="X.Nonlinear")]<0.001,"<0.001",round(a$P[which(rownames(a)=="X.Nonlinear")],digits = 3))
p_nonlinear  #<0.001

HR$yhat = exp(log(HR$yhat)/lambda)
HR$lower = exp(log(HR$lower)/lambda)
HR$upper = exp(log(HR$upper)/lambda)

p_major_BUN4 <-ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), data_c2_nona,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,35,1)) +  
  geom_line(aes(x = BUN, y = yhat), HR, linetype = 1, col = mycolors[9], size = 0.8) + 
  geom_line(aes(x = BUN, y=lower), HR, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = BUN, y=upper),HR, color=mycolors[8],linetype='dashed') +
  geom_hline(yintercept = 1, linetype = 8, size = 0.5,color="grey40") + 
  # geom_point(aes(x = PM2.5_knots, y = PM2.5_yhat)) + 
  scale_x_reverse(breaks = seq(0,35,2))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,6,1),  
                     sec.axis = sec_axis(~. / 10, name = expression(paste("Density of ","BUN")))) + 
  coord_cartesian(xlim = c(0,35), ylim = c(0,6)) +
  labs(x = "BUN", y = "Hazard Ratio",title = "BUN and Skin and Soft Tissue Infections") + 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text( size = 10,vjust = 0.5), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")),    
        plot.tag = element_text( size = 14, face = "bold"),
        plot.title = element_text(size = 13, face = "bold")) + 
  annotate("text", x = 8, y = 4.3, label = paste("P(overall)", p_overall), size = 3.5)+
  annotate("text", x = 8, y = 4.7, label = paste("P(nonlinearity)", p_nonlinear),  size = 3.5)
p_major_BUN4



#urinary_tract

all_data$event.type <- all_data$urinary_tract_event
all_data$follow_up_time<-all_data$urinary_tract_time

lambda=1
mycolors<-brewer.pal(9, "Blues")
covar = c("age","sex","education","TDI","smoking","alcohol","CVD","chronic_pulmonary_disease",
          "cancer","diabetes","hypertension","UACR_cat","BMI","medicine_g_i")

#d
####################################################
##eGFR_scr
need_term <- c(covar,"eGFR_scr")
data_c2_nona <- all_data[complete.cases(all_data[,c(need_term)]),c(need_term,"event.type","follow_up_time")]
dd <- datadist(data_c2_nona) #为后续程序设定数据环境
options(datadist='dd')

sm1.1 = cph(Surv(follow_up_time, event.type) ~ age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
              cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i+
              rcs(eGFR_scr,quantile(eGFR_scr, c( .05,0.275, .5,0.725,  .95))),
            #rcs(eGFR_10,c( 4.5, 6.0, 7.5, 9.0, 10.5 )),
            data = data_c2_nona)
fit=update(sm1.1)
#
HR<-Predict(fit,eGFR_scr,fun=exp,ref.zero = TRUE)

a <- as.data.frame(anova(fit))
a
p_overall <- ifelse(a$P[which(rownames(a)=="eGFR_scr")]<0.001,"<0.001",round(a$P[which(rownames(a)=="eGFR_scr")],digits = 3))
p_overall  #<0.001
p_nonlinear <- ifelse(a$P[which(rownames(a)=="X.Nonlinear")]<0.001,"<0.001",round(a$P[which(rownames(a)=="X.Nonlinear")],digits = 3))
p_nonlinear  #<0.001

HR$yhat = exp(log(HR$yhat)/lambda)
HR$lower = exp(log(HR$lower)/lambda)
HR$upper = exp(log(HR$upper)/lambda)

p_major_eGFRscr5 <-ggplot() + 
  geom_histogram(aes(x = eGFR_scr, y = (..count..) / sum(..count..) * 10), data_c2_nona,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(5,155,5)) +  
  geom_line(aes(x = eGFR_scr, y = yhat), HR, linetype = 1, col = mycolors[9], size = 0.8) + 
  geom_line(aes(x = eGFR_scr, y=lower), HR, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = eGFR_scr, y=upper),HR, color=mycolors[8],linetype='dashed') + 
  geom_hline(yintercept = 1, linetype = 8, size = 0.5,color="grey40") + 
  # geom_point(aes(x = PM2.5_knots, y = PM2.5_yhat)) + 
  scale_x_reverse(breaks = seq(10,160,10))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,6,1),  
                     sec.axis = sec_axis(~. / 10, name = expression(paste("Density of ",eGFR["cr"])))) + 
  coord_cartesian(xlim = c(160,10), ylim = c(0,6)) +
  labs(x = expression(paste(eGFR["cr"]," (","ml/min/1.73 ",m^2,")")), y = "Hazard Ratio",title = expression(paste(eGFR["cr"]," and Urinary Tract Infections"))) + 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text( size = 10,vjust = 0.5), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")),    
        plot.tag = element_text( size = 14, face = "bold"),
        plot.title = element_text(size = 15, face = "bold")) + 
  annotate("text", x = 40, y = 5.3, label = paste("P(overall)", p_overall), size = 3.5)+
  annotate("text", x = 40, y = 5.7, label = paste("P(nonlinearity)", p_nonlinear),  size = 3.5)
p_major_eGFRscr5


#d
####################################################
##eGFR_scys
need_term <- c(covar,"eGFR_scys")
data_c2_nona <- all_data[complete.cases(all_data[,c(need_term)]),c(need_term,"event.type","follow_up_time")]
dd <- datadist(data_c2_nona) #为后续程序设定数据环境
options(datadist='dd')

sm1.2 = cph(Surv(follow_up_time, event.type) ~ age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
              cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i+
              rcs(eGFR_scys,quantile(eGFR_scys, c( .05,0.275, .5,0.725,  .95))),
            #rcs(eGFR_10,c( 4.5, 6.0, 7.5, 9.0, 10.5 )),
            data = data_c2_nona)
fit=update(sm1.2)
#
HR<-Predict(fit,eGFR_scys,fun=exp,ref.zero = TRUE)

a <- as.data.frame(anova(fit))
a
p_overall <- ifelse(a$P[which(rownames(a)=="eGFR_scys")]<0.001,"<0.001",round(a$P[which(rownames(a)=="eGFR_scys")],digits = 3))
p_overall  #<0.001
p_nonlinear <- ifelse(a$P[which(rownames(a)=="X.Nonlinear")]<0.001,"<0.001",round(a$P[which(rownames(a)=="X.Nonlinear")],digits = 3))
p_nonlinear  #<0.001

HR$yhat = exp(log(HR$yhat)/lambda)
HR$lower = exp(log(HR$lower)/lambda)
HR$upper = exp(log(HR$upper)/lambda)

p_major_eGFRscys5 <-ggplot() + 
  geom_histogram(aes(x = eGFR_scys, y = (..count..) / sum(..count..) * 10), data_c2_nona,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(5,155,5)) +  
  geom_line(aes(x = eGFR_scys, y = yhat), HR, linetype = 1, col = mycolors[9], size = 0.8) + 
  geom_line(aes(x = eGFR_scys, y=lower), HR, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = eGFR_scys, y=upper),HR, color=mycolors[8],linetype='dashed') +
  geom_hline(yintercept = 1, linetype = 8, size = 0.5,color="grey40") + 
  # geom_point(aes(x = PM2.5_knots, y = PM2.5_yhat)) + 
  scale_x_reverse(breaks = seq(10,160,10))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,6,1),  
                     sec.axis = sec_axis(~. / 10, name = expression(paste("Density of ",eGFR["cys"])))) + 
  coord_cartesian(xlim = c(160,10), ylim = c(0,6)) +
  labs(x = expression(paste(eGFR["cys"]," (","ml/min/1.73 ",m^2,")")), y = "Hazard Ratio",title = expression(paste(eGFR["cys"]," and Urinary Tract Infections"))) + 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text( size = 10,vjust = 0.5), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")),    
        plot.tag = element_text( size = 14, face = "bold"),
        plot.title = element_text(size = 15, face = "bold")) + 
  annotate("text", x = 120, y = 5.3, label = paste("P(overall)", p_overall), size = 3.5)+
  annotate("text", x = 120, y = 5.7, label = paste("P(nonlinearity)", p_nonlinear),  size = 3.5)
p_major_eGFRscys5


####################################################
##BUN
need_term <- c(covar,"BUN")
data_c2_nona <- all_data[complete.cases(all_data[,c(need_term)]),c(need_term,"event.type","follow_up_time")]
dd <- datadist(data_c2_nona) #为后续程序设定数据环境
options(datadist='dd')

sm1.3 = cph(Surv(follow_up_time, event.type) ~ age+sex+education+TDI+smoking+alcohol+CVD+chronic_pulmonary_disease+
              cancer+diabetes+hypertension+UACR_cat+BMI+medicine_g_i+
              rcs(BUN,quantile(BUN, c( .05,0.275, .5,0.725,  .95))),
            #rcs(eGFR_10,c( 4.5, 6.0, 7.5, 9.0, 10.5 )),
            data = data_c2_nona)
fit=update(sm1.3)
#
HR<-Predict(fit,BUN,fun=exp,ref.zero = TRUE)

a <- as.data.frame(anova(fit))
a
p_overall <- ifelse(a$P[which(rownames(a)=="BUN")]<0.001,"<0.001",round(a$P[which(rownames(a)=="BUN")],digits = 3))
p_overall  #<0.001
p_nonlinear <- ifelse(a$P[which(rownames(a)=="X.Nonlinear")]<0.001,"<0.001",round(a$P[which(rownames(a)=="X.Nonlinear")],digits = 3))
p_nonlinear  #<0.001

HR$yhat = exp(log(HR$yhat)/lambda)
HR$lower = exp(log(HR$lower)/lambda)
HR$upper = exp(log(HR$upper)/lambda)

p_major_BUN5 <-ggplot() + 
  geom_histogram(aes(x = BUN, y = (..count..) / sum(..count..) * 10), data_c2_nona,# binwidth = 5, 
                 fill = mycolors[4], color = "grey80", size = 0.5, alpha = 0.1,
                 breaks=seq(0,35,1)) +  
  geom_line(aes(x = BUN, y = yhat), HR, linetype = 1, col = mycolors[9], size = 0.8) + 
  geom_line(aes(x = BUN, y=lower), HR, color=mycolors[8],linetype='dashed') +
  geom_line(aes(x = BUN, y=upper),HR, color=mycolors[8],linetype='dashed') +
  geom_hline(yintercept = 1, linetype = 8, size = 0.5,color="grey40") + 
  # geom_point(aes(x = PM2.5_knots, y = PM2.5_yhat)) + 
  scale_x_continuous(breaks = seq(0,35,2))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,6,1),  
                     sec.axis = sec_axis(~. / 10, name = expression(paste("Density of ","BUN")))) + 
  coord_cartesian(xlim = c(0,35), ylim = c(0,6)) +
  labs(x = "BUN", y = "Hazard Ratio",title = "BUN and Urinary Tract Infections") + 
  theme(panel.background = element_rect(fill = "transparent"), 
        axis.line.x.bottom = element_line(color = "black"), 
        axis.line.y.left = element_line(color = "black"), 
        axis.line.y.right = element_line(color = "black"), 
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text( size = 10,vjust = 0.5), 
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 3), "mm")),    
        plot.tag = element_text( size = 14, face = "bold"),
        plot.title = element_text(size = 15, face = "bold")) + 
  annotate("text", x = 8, y = 4.3, label = paste("P(overall)", p_overall), size = 3.5)+
  annotate("text", x = 8, y = 4.7, label = paste("P(nonlinearity)", p_nonlinear),  size = 3.5)
p_major_BUN5


pdf("Restricted_cubic_spline_new.pdf",family = "Times",height = 15,width = 15)

plot_grid(p_major_eGFRscr1,p_major_eGFRscys1,p_major_BUN1,
          p_major_eGFRscr2,p_major_eGFRscys2,p_major_BUN2,
          p_major_eGFRscr3,p_major_eGFRscys3,p_major_BUN3,
          p_major_eGFRscr4,p_major_eGFRscys4,p_major_BUN4,
          p_major_eGFRscr5,p_major_eGFRscys5,p_major_BUN5,
          nrow = 5,ncol =3)
dev.off()



