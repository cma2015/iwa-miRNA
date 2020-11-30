options(stringsAsFactors = F)
table_data <- read.table("final_table.txt", header=T, sep="\t")
miRNAs_list <- read.table("mature_miRNAs.txt", header=F, sep="\t", row.names = 1)
psRNAtarget_output <- read.table("psRNAtarget_MIT.out", header=T, sep="\t")

for(ii in unique(table_data[,5])){
    tmp_vec <- c(paste0(ii, "_5p"), paste0(ii, "_3p"))
    names(tmp_vec) <- miRNAs_list[tmp_vec,1]
    tmp_target <- psRNAtarget_output[psRNAtarget_output[,1]%in%names(tmp_vec), ]
    tmp_target[,1] <- tmp_vec[tmp_target[,1]]
    write.table(tmp_target, paste0("../miRNASelection/data/", ii, ".mti"), quote=F, sep="\t", row.names = F)
}
