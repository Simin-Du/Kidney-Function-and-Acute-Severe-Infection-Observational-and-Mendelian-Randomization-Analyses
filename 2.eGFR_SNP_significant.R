library(readxl)
library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

sTable6 <- read_excel("NIHMS1038587-supplement-Tables.xlsx", sheet = 7, range = cell_rows(2:266))

sTable10 <- read_excel("NIHMS1038587-supplement-Tables.xlsx", sheet = 11, range = cell_rows(2:258))
sTable10 <- sTable10[sTable10$`Gray font` == 0, ]

sTable6$Chr <- as.integer(sub(":.*", "", sTable6$`Chr/Pos (b37)`))
sTable6$Pos <- as.integer(sub(".*:", "", sTable6$`Chr/Pos (b37)`))

sTable10$Chr <- as.integer(sub(":.*", "", sTable10$`Chr/Pos (b37)`))
sTable10$Pos <- as.integer(sub(".*:", "", sTable10$`Chr/Pos (b37)`))

sTable6$Pos_50kb_downstream <- sTable6$Pos - 50000
sTable6$Pos_50kb_upstream <- sTable6$Pos + 50000

range(unique(sTable6$Chr))
range(unique(sTable10$Chr))

# -------------------------------------------------------------------------
match_in_sTable6 <- c()
snp1_in_sTable6 <- c()
snp2_in_sTable6 <- c()

for (Chr in 1:22) {
  
  temp1 <- lapply(sTable10$Pos[sTable10$Chr == Chr], function(x) {
    x >= sTable6$Pos_50kb_downstream[sTable6$Chr == Chr] & 
      x <= sTable6$Pos_50kb_upstream[sTable6$Chr == Chr]})
  
  temp2 <- as.data.frame(matrix(unlist(temp1), byrow = T, nrow = length(sTable10$Pos[sTable10$Chr == Chr])))
  row.names(temp2) <- sTable10$`RS number`[sTable10$Chr == Chr]
  colnames(temp2) <- sTable6$`RS number`[sTable6$Chr == Chr]
  
  match_in_sTable6 <- c(match_in_sTable6, apply(temp2, 1, sum))

  snp1_in_sTable6 <- c(snp1_in_sTable6, sTable6$`RS number`[sTable6$Chr == Chr][unlist(lapply(apply(temp2, 1, function(x) {which(x == T)}), `[`, 1))])
  snp2_in_sTable6 <- c(snp2_in_sTable6, sTable6$`RS number`[sTable6$Chr == Chr][unlist(lapply(apply(temp2, 1, function(x) {which(x == T)}), `[`, 2))])
}

sTable10$match_in_sTable6 <- match_in_sTable6
sum(sTable10$match_in_sTable6 > 0)  ## n=185

sTable10$snp1_in_sTable6 <- snp1_in_sTable6
sTable10$snp2_in_sTable6 <- snp2_in_sTable6

sTable10$`BUN p-value 1` <- sTable6$`BUN p-value`[match(sTable10$snp1_in_sTable6, sTable6$`RS number`)]
sTable10$`BUN p-value 2` <- sTable6$`BUN p-value`[match(sTable10$snp2_in_sTable6, sTable6$`RS number`)]

sTable10$`BUN p-value` <- apply(sTable10[, c("BUN p-value 1", "BUN p-value 2")], 1, function(x) {min(x, na.rm = T)})
sTable10$snp_in_sTable6 <-ifelse(sTable10$`BUN p-value 1` == sTable10$`BUN p-value`, sTable10$snp1_in_sTable6, sTable10$snp2_in_sTable6)

# -------------------------------------------------------------------------
SNP_rawdata=read.csv("SNP_table10.csv")
remove_SNP=read.csv("remove_SNP.csv",header=T)
SNPlist=SNP_rawdata[-which(SNP_rawdata$RS.number %in% remove_SNP$removeSNP),"RS.number"]
k=length(SNPlist)

eGFR_SNP_significant <- intersect(sTable10$`RS number`[sTable10$match_in_sTable6 > 0 & sTable10$`BUN p-value` < 0.05 / k & sTable10$`Support for kidney function relevance`!=3], SNPlist)
write.table(eGFR_SNP_significant, "eGFR_SNP_significant.txt", row.names = F, col.names = F, quote = F)


##################

