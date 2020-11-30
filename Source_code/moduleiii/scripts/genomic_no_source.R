
options(stringsAsFactors = F)
table_data <- read.table("final_table.txt", header=T, sep="\t")

table_data_add <- cbind(table_data[,1:4], "-", table_data[,5:ncol(table_data)])
colnames(table_data_add)[5] <- "Genomic_source"

write.table(table_data_add, "final_table_source.txt", quote=F, sep="\t", row.names = F)
