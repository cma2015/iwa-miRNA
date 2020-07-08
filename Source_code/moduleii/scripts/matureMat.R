
options(stringsAsFactors = F)
library(reshape2)

# sample_information <- read.table("sample_info.txt", row.names = 1, sep = "\t")

expre_mat <- read.table("seq_mat.txt", sep = "\t")
colnames(expre_mat) <- c("SRR", "Seq", "abundance")
abundance_mat <- dcast(expre_mat, Seq~SRR, value.var = "abundance", )
rownames(abundance_mat) <- abundance_mat[,1]
abundance_mat <- abundance_mat[,-1]

abundance_mat[is.na(abundance_mat)] <- 0

exp_mat <- abundance_mat

write.table(exp_mat, "miRNA_in_sample.txt", quote = F, sep = "\t")

# sam_name <- unique(sample_information[,1])
# sam_name <- sam_name[!sam_name%in%"-"]
# #sam_name <- intersect(c("root", "stem", "leaf", "flower", "seed"), sam_name)
# sam_name
# merge_mat <- matrix(nrow = nrow(exp_mat), ncol = length(sam_name))
# rownames(merge_mat) <- rownames(exp_mat)
# colnames(merge_mat) <- sort(sam_name)
# 
# for(ii in rownames(merge_mat)){
#   for(jj in colnames(merge_mat)){
#     srr_index <- rownames(sample_information)[sample_information[,1]==jj]
#     srr_index <- intersect(srr_index, colnames(exp_mat))
#     merge_mat[ii, jj] <- round(mean(as.numeric(exp_mat[ii, srr_index])), 2)
#   }
# }

#tmp_order <- read.table("sequences.txt")

#merge_mat <- as.data.frame(merge_mat, stringsAsFactors = F)
#merge_mat <- merge_mat[tmp_order[,1], ]
#rownames(merge_mat) <- tmp_order[,1]
#merge_mat[is.na(merge_mat)] <- 0

# write.table(merge_mat, "./miRNA_in_tissue.txt", quote = F, sep = "\t")

