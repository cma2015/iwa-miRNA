#!/usr/bin/env Rscript
options(stringsAsFactors = F)

args <- commandArgs(T)

args <- c("final_table.txt", "mir.gff", "mir.fasta")

info_df <- read.table(args[1], header = T, sep = "\t", stringsAsFactors = F)
if("Extended_stem_loop_loc" %in% colnames(info_df)){
  tmp_df <- do.call("rbind", strsplit(info_df$Extended_stem_loop_loc, ":"))
  tmp_df <- cbind(tmp_df[,1], info_df$Extended_stem_loop_loc, tmp_df[,2:3])
  
  write.table(tmp_df, args[2], row.names = F, quote = F, sep = "\t", col.names = F)
  
  tmp_vec <- tmp_id <- vector()
  for(i in 1:nrow(info_df)){
    if( ! info_df[i,]$Extended_stem_loop_loc %in% tmp_id){
      tmp_vec <- c(tmp_vec, paste0(">", info_df[i,]$Extended_stem_loop_loc), info_df[i,]$Extended_stem_loop_seq)
      tmp_id <- c(tmp_id, info_df[i,]$Extended_stem_loop_loc)
    }
  }
  write.table(tmp_vec, args[3], row.names = F, quote = F, sep = "\t", col.names = F)
}else{
  tmp_df <- do.call("rbind", strsplit(info_df$pLoc, ":"))
  tmp_df <- cbind(tmp_df[,1], info_df$pLoc, tmp_df[,2:3])
  
  write.table(tmp_df, args[2], row.names = F, quote = F, sep = "\t", col.names = F)
  
  tmp_vec <- tmp_id <- vector()
  for(i in 1:nrow(info_df)){
    if( ! info_df[i,]$pLoc %in% tmp_id){
      tmp_vec <- c(tmp_vec, paste0(">", info_df[i,]$pLoc), info_df[i,]$pSeq)
      tmp_id <- c(tmp_id, info_df[i,]$pLoc)
    }
  }
  write.table(tmp_vec, args[3], row.names = F, quote = F, sep = "\t", col.names = F) 
}
