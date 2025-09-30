setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


SNP_list=read.table("eGFR_SNP_significant.txt")

########################################################################################################################
# clumping
#library(TwoSampleMR)
SNP_data=read.csv("SNP_table10.csv")
SNP_data$Chr <- as.integer(sub(":.*", "", SNP_data$Chr.Pos..b37.))
SNP_data$Pos <- as.integer(sub(".*:", "", SNP_data$Chr.Pos..b37.))

SNP_significant_data=merge(SNP_data,SNP_list,by.x = "RS.number",by.y="V1")

for (i in 1:length(SNP_significant_data$RS.number)){
  for (j in i:length(SNP_significant_data$RS.number)){
    if(i!=j & SNP_significant_data$Chr[i]==SNP_significant_data$Chr[j] &abs(SNP_significant_data$Pos[i]-SNP_significant_data$Pos[j])<=500000){
      print("pairs:")
      print(SNP_significant_data$RS.number[i])
      print(SNP_significant_data$RS.number[j])
    }
  }
}



write.table(SNP_list,"SNP_after_clumping.txt",quote=F,row.names = F)


########################################################################################################

#test
SNP_list=read.table("SNP_after_clumping.txt",header = T)
SNP_list_new=read.table("27SNP_1Proxy.txt",header=T)
all_SNP=read.table("ALLSNPrsID.txt")

for (i in SNP_list$V1){
  if (!(i %in% all_SNP$V1)){
    print(i)
  }
}









