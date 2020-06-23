#!/usr/bin/env Rscript
options(stringsAsFactors = F)

args <- commandArgs(T)

# args <- c("../tmp/OnX4IRTGlT/out_pool_merge.txt", "../tmp/OnX4IRTGlT/pc_criteria.txt")

info_df <- read.table(args[1], header = T, sep = "\t", stringsAsFactors = F)
rownames(info_df) <- info_df[,2]

#tmp_info_df <- as.character(info_df[, c(9,12,17,23,25,26)])
as.character(info_df[1,])
tmp_info_df <- info_df[, c(9,12,13,17,23,25,26)]
tmp_info_df[tmp_info_df[,5]=="-", 5] <- 0
tmp_info_df[1:5,]
pc_criteria <- apply(tmp_info_df, 1, function(x){
  if(x[3]=="3p"){
    len_criterion <- as.numeric(x[2])>=20&as.numeric(x[2])<=24
  } else{
    len_criterion <- as.numeric(x[1])>=20&as.numeric(x[1])<=24
  }
  return(as.numeric(x[4])<300&as.numeric(x[5])>0.75&
           x[6]==TRUE&x[7]==TRUE&len_criterion) #&as.numeric(x[21])>0.9
})

pc_criteria_name <- rownames(info_df)[pc_criteria]
write.table(pc_criteria, args[2], quote = F, sep = "\t", col.names = F)
