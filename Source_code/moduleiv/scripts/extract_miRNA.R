## R for miRNA bed and mature location
options(stringsAsFactors = F)

data_input <- read.table("final_table.txt", sep = "\t", header = T, stringsAsFactors = F)

if("Extended_stem_loop_loc" %in% colnames(data_input)){
  pre_miRNAs <- do.call('rbind', strsplit(data_input$Extended_stem_loop_loc, ":"))
  pre_miRNAs <- cbind(pre_miRNAs[, 1:3], data_input$ID, ".", pre_miRNAs[,4])
  write.table(pre_miRNAs, "pre_miRNA.bed", quote = F, 
              row.names = F, col.names = F, sep = "\t")
  
  mature_loc <- t(apply(data_input, 1, function(x){
    preseq <- strsplit(x[5], ":")[[1]]
    mature <- strsplit(x[10], ":")[[1]]
    star <- strsplit(x[13], ":")[[1]]
    if(preseq[4] == "+"){
      locval <- as.numeric(c(mature[2:3],star[2:3]))-as.numeric(preseq[2])
    }else{
      locval <- as.numeric(preseq[3])-as.numeric(c(mature[2:3],star[2:3]))
    }
    c(x[2], sort(locval))
  }))
  write.table(mature_loc, "mature_location.txt", quote = F, 
              row.names = F, col.names = F, sep = "\t")
} else if ("Precursors" %in% colnames(data_input) ){
  pre_miRNAs <- do.call('rbind', strsplit(data_input$Precursors, ":"))
  pre_miRNAs <- cbind(pre_miRNAs[, 1:3], data_input$Precursors, ".", pre_miRNAs[,4])
  write.table(pre_miRNAs, "pre_miRNA.bed", quote = F, 
              row.names = F, col.names = F, sep = "\t")
  
  mature_loc <- t(apply(data_input, 1, function(x){
    preseq <- strsplit(x[2], ":")[[1]]
    mature <- strsplit(x[5], ":")[[1]]
    star <- strsplit(x[8], ":")[[1]]
    if(preseq[4] == "+"){
      locval <- as.numeric(c(mature[2:3],star[2:3]))-as.numeric(preseq[2])
    }else{
      locval <- as.numeric(preseq[3])-as.numeric(c(mature[2:3],star[2:3]))
    }
    c(x[2], sort(locval))
  }))
  write.table(mature_loc, "mature_location.txt", quote = F, 
              row.names = F, col.names = F, sep = "\t")
}

snp_input <- read.table("SNPs.txt", sep = "\t", header = T, stringsAsFactors = F)
snp_input <- cbind(snp_input[, c(1,2,2,3)], ".", ".", snp_input[,4:5])
write.table(snp_input, "SNP.bed", quote = F, 
            row.names = F, col.names = F, sep = "\t")
  