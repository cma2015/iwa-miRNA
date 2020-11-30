#!/usr/bin/env Rscript
options(stringsAsFactors = F)

args <- commandArgs(T)
# args <- c("out_pool_merge.txt", "pc_criteria.txt")

info_df <- read.table(args[1], header = T, sep = "\t", stringsAsFactors = F)
rownames(info_df) <- info_df[,2]

#tmp_info_df <- as.character(info_df[, c(9,12,17,23,25,26)])
as.character(info_df[1,])
tmp_info_df <- info_df[, c(9,12,13,17,23,25,26)]
tmp_info_df <- info_df[, c("Len5p", "Len3p", "Mature_arm", "Stem_loop_len", 
                           "Abundance_bias", "RNAfold", "Centroidfold")]
tmp_info_df[tmp_info_df[,5]=="-", 5] <- 0

stemloop <- args[3]
structure <- args[4]
abias <- args[5]
sbias <- args[6]
minlen <- args[7]
maxlen <- args[8]

cat(stemloop, structure, abias, sbias, minlen, maxlen)

pc_criteria <- apply(tmp_info_df, 1, function(x){
  if(x[3]=="3p"){
    len_criterion <- as.numeric(x[2])>=as.numeric(minlen)&as.numeric(x[2])<=as.numeric(maxlen)
  } else{
    len_criterion <- as.numeric(x[1])>=as.numeric(minlen)&as.numeric(x[1])<=as.numeric(maxlen)
  }
  if(structure == "no"){
    return(as.numeric(x[4])<as.numeric(stemloop)&as.numeric(x[5])>as.numeric(abias)&len_criterion) #&as.numeric(x[21])>0.9
  }else{
    return(as.numeric(x[4])<as.numeric(stemloop)&as.numeric(x[5])>as.numeric(abias)&
             x[6]=='True'&x[7]=='True'&len_criterion) #&as.numeric(x[21])>0.9
  }
})

pc_criteria_name <- rownames(info_df)[pc_criteria]
write.table(pc_criteria, args[2], quote = F, sep = "\t", col.names = F)
