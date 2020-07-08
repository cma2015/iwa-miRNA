#!/usr/bin/env Rscript
options(stringsAsFactors = F)

args <- commandArgs(T)
# args <- c("Translate_result.txt",
#           "mir.txt")

info_df <- read.table(args[1], header = T, sep = "\t", stringsAsFactors = F)
if("Extended_stem_loop_loc"%in%colnames(info_df)){
  tmp_df <- do.call("rbind", strsplit(info_df$Extended_stem_loop_loc, ":"))
}else{
  tmp_df <- do.call("rbind", strsplit(info_df$pLoc, ":"))
}
mature_len <- apply(info_df, 1, function(x){
  arm_index <- colnames(info_df)%in%'Mature_arm'
  if(x[arm_index]=="3p"){
    arm_len_index <- colnames(info_df)%in%'Len3p'
    len_res <- as.numeric(x[arm_len_index])
  } else{
    arm_len_index <- colnames(info_df)%in%'Len5p'
    len_res <- as.numeric(x[arm_len_index])
  }
  return(len_res)
})

if("ID"%in%colnames(info_df)){
  tmp_df <- cbind(tmp_df[,1:3], info_df$ID, mature_len, tmp_df[,4], info_df$Mean)
}else{
  tmp_df <- cbind(tmp_df[,1:3], info_df$Precursors, mature_len, tmp_df[,4], "-")
}

write.table(tmp_df, args[2], row.names = F, quote = F, sep = "\t", col.names = F)
