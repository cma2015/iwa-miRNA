##
library(dplyr)
options(stringsAsFactors = F)
args <- commandArgs(T)

resultdf <- read.table(file = args[3], sep = "\t", header = T)
rownames(resultdf) <- resultdf[,2]

pc_result <- read.table(file = args[1], sep = "\t", header = F, row.names = 1)
ml_result <- read.table(file = args[2], sep = "\t", header = F, row.names = 1)
rename_result <- read.table(file = args[4], sep = "\t", header = T)
exp_result <- read.table(file = args[5], sep = "\t", header = T)
exp_summary <- data.frame("Mean" = apply(exp_result, 1, mean),
                          "Max" = apply(exp_result, 1, max),
                          "Sample" = apply(exp_result, 1, function(x){sum(x>1)}))
rownames(exp_summary) <- rownames(exp_result)

resultdf$"HTcriteria" <- pc_result[rownames(resultdf), 1]
resultdf$"One_class_SVM" <- ml_result[rownames(resultdf), 1]
resultdf$ID <- rename_result[,2]
resultdf$ID[is.na(resultdf$ID)] <- resultdf$Pre.miRNAs[is.na(resultdf$ID)]
resultdf$Mean <- apply(exp_result, 1, mean)[rownames(resultdf)]
resultdf$Max <- apply(exp_result, 1, max)[rownames(resultdf)]
resultdf$samples <- apply(exp_result, 1, function(x){sum(x>1)})[rownames(resultdf)]
resultdf[is.na(resultdf)] <- 0

resultdf <- resultdf[, c(1, 29, 27:28, 2:26, 30:ncol(resultdf))]
write.table(resultdf, args[6], quote = F, row.names = F, sep = "\t")
