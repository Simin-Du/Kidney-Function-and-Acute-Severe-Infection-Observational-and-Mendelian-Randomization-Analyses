setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
SNP_data=read.csv("SNP_table10.csv")
#SNP_trait=read.csv("traits.csv")
#library(devtools)
#install_github("phenoscanner/phenoscanner")
library(phenoscanner)

trait=c()

colnames(SNP_data)

snplist = SNP_data[,"RS.number"]
res1 <- phenoscanner(snpquery=snplist[1:100],pvalue = 1e-08,proxies = "EUR",r2 = 0.8)
res2 <- phenoscanner(snpquery=snplist[101:200],pvalue = 1e-08,proxies = "EUR",r2 = 0.8)
res3 <- phenoscanner(snpquery=snplist[201:256],pvalue = 1e-08,proxies = "EUR",r2 = 0.8)
related_traits1 = unique(sort(c(res1$results$trait,res2$results$trait,res3$results$trait)))
#write.csv(related_traits1,"256SNP_traits1.csv")

related_traits1=read.csv("256SNP_traits1_20221114.csv",header = T)
traits1=related_traits1[which(related_traits1$related==1),"x"]
#如果related_traits多的话，把这个存到变量存到csv中，让黄主任选

#以下针对黄主任选出的变量：
selected_traits = traits1
removeSNP = c()
for(i in 1:length(selected_traits)){
  temp = grep(selected_traits[i],res1$results$trait)
  if(sum(temp,na.rm = T)!=0){
    removeSNP = c(removeSNP,res1$results$snp[temp])
  }
}
for(i in 1:length(selected_traits)){
  temp = grep(selected_traits[i],res2$results$trait)
  if(sum(temp,na.rm = T)!=0){
    removeSNP = c(removeSNP,res2$results$snp[temp])
  }
}
for(i in 1:length(selected_traits)){
  temp = grep(selected_traits[i],res3$results$trait)
  if(sum(temp,na.rm = T)!=0){
    removeSNP = c(removeSNP,res3$results$snp[temp])
  }
}
removeSNP1=unique(sort(removeSNP))#11


#######################################################################################

traits_data2=read.csv("traits2_all.csv")
related_traits2=unique(sort(traits_data2$GWAS.Trait))
#write.csv(related_traits2,"256SNP_traits2.csv")

related_traits2=read.csv("256SNP_traits2_20221114.csv")

traits2=related_traits2[which(related_traits2$related==1),"x"]
#如果related_traits多的话，把这个存到变量存到csv中，让黄主任选

#以下针对黄主任选出的变量：
selected_traits = traits2
removeSNP = c()
for(i in 1:length(selected_traits)){
  temp = grep(selected_traits[i],traits_data2$GWAS.Trait)
  if(sum(temp,na.rm = T)!=0){
    removeSNP = c(removeSNP,traits_data2$Query[temp])
  }
}
removeSNP2=unique(sort(removeSNP))

removeSNP=unique(sort(c(removeSNP1,removeSNP2)))
removeSNP=data.frame("removeSNP"=removeSNP)

write.csv(removeSNP,"remove_SNP.csv",quote = F,row.names = F)







