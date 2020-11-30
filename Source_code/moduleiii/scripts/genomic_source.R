
options(stringsAsFactors = F)
TE_inter <- read.csv("miRNA_TE.txt", sep = "\t", header = F)
TE_inter_names <- unique(TE_inter[,4])
non_TE_inter <- read.csv("miRNA_non_TE.txt", sep = "\t", header = F)
non_TE_inter <- non_TE_inter[!non_TE_inter[,4]%in%TE_inter_names, ]
non_TE_PCG_names <- non_TE_inter[!grepl("[Pp][CcEe][Gg]", non_TE_inter[,14]), 4]
PCG_index <- grepl("[Pp][CcEe][Gg]", non_TE_inter[,14])
non_TE_PCG_index <-  !(PCG_index&non_TE_inter[,4]%in%non_TE_PCG_names)
non_TE_inter <- non_TE_inter[non_TE_PCG_index, ]

new_table <- rbind(TE_inter, non_TE_inter)
all_mirs <- read.csv("miRNA_list.txt", sep = "\t", header = F)
all_mirs <- all_mirs[!all_mirs[,4]%in%new_table[,4], ]
all_mirs[,(ncol(all_mirs)+1):ncol(new_table)] <- "."
data <- rbind(new_table, all_mirs)
tmp_names <- apply(data, 1, function(x){paste0(c(x[1:3], x[6]), collapse = ":")})
uni_index <- vector()
for(i in unique(tmp_names)){
  uni_index <- c(uni_index, which(tmp_names%in%i)[1])
}
data <- data[uni_index,]
data[data[,14]==".", 14] <- "Intergenic"

rownames(data) <- apply(data, 1, function(x){paste0(c(x[1:3], x[6]), collapse = ":")}) #4:11795179:11795270:+
rownames(data) <- gsub(" ", "", rownames(data))
table_data <- read.table("final_table.txt", header=T, sep="\t")

table_data_add <- cbind(table_data[,1:4], data[table_data[,5], 14], table_data[,5:ncol(table_data)])
colnames(table_data_add)[5] <- "Genomic_source"

write.table(table_data_add, "final_table_source.txt", quote=F, sep="\t", row.names = F)
