
library(reshape2)

data_input <- read.table("corr_pre-mirna_seq_mat.txt", sep = "\t", stringsAsFactors = F)


all_abundance <- data_input[data_input[,1]=="abundance"&data_input[,2]=="all", 3:5]
colnames(all_abundance) <- c("SRR", "Loc", "abundance")
all_abundance_mat <- dcast(all_abundance, Loc~SRR, value.var = "abundance", )

all_abundance_mat[is.na(all_abundance_mat)] <- 0
all_abundance_mat[1:5,]

write.table(all_abundance_mat, "expressionMat.txt", quote = F, sep = "\t", row.names = F)

# range(apply(all_abundance_mat, 1, max))
# 
# limit_abundance <- data_input[data_input[,1]=="abundance"&data_input[,2]=="limit", 3:5]
# colnames(limit_abundance) <- c("SRR", "Loc", "abundance")
# limit_abundance_mat <- dcast(limit_abundance, Loc~SRR, value.var = "abundance", )
# rownames(limit_abundance_mat) <- limit_abundance_mat[,1]
# limit_abundance_mat <- limit_abundance_mat[, -1]
# limit_abundance_mat <- as.matrix(limit_abundance_mat)
# limit_abundance_mat[is.na(limit_abundance_mat)] <- 0
# limit_abundance_mat[1:5,1:5]
# range(apply(limit_abundance_mat, 1, max))
# 
# save(limit_abundance_mat, all_abundance_mat, file = "expMats.RData")
#   
