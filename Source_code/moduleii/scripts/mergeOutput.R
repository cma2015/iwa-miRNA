##
library(dplyr)
options(stringsAsFactors = F)
args <- commandArgs(T)
# args <- c("pc_criteria.txt", "ML_result.txt", "out_pool_merge.txt",
#           "NameChange_out.txt", "expressionMat.txt", "final_table.txt")

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
resultdf$"ID" <- rename_result[,2]
resultdf$"ID"[is.na(resultdf$"ID")] <- resultdf$Precursors[is.na(resultdf$"ID")]
resultdf$Mean <- apply(exp_result, 1, mean)[rownames(resultdf)]
resultdf$Max <- apply(exp_result, 1, max)[rownames(resultdf)]
resultdf$Samples <- apply(exp_result, 1, function(x){sum(x>1)})[rownames(resultdf)]
resultdf[is.na(resultdf)] <- 0

output_name <- c("Precursors", "ID", "HTcriteria", "One_class_SVM", "Extended_stem_loop_loc", "Extended_stem_loop_seq",
                 "Extended_stem_loop_len", "Stem_loop_loc", "Stem_loop_seq", "Loc5p", "Seq5p", "Len5p", "Loc3p", "Seq3p",
                 "Len3p", "Mature_arm", "Source", "RPM5p", "RPM3p", "Stem_loop_len", "Stem_loop_MFE", "Stem_loop_AMFE", 
                 "The_total_abundance","The_number_of_sequences_in_miRNA.miRNA._and_3nt_variant_region",
                 "The_number_of_sequences_in_pre.miRNAs", "Abundance_bias", "Strand_bias", "RNAfold", "Centroidfold",
                 "Mean", "Max", "Samples")
resultdf <- resultdf[, output_name]
write.table(resultdf, args[6], quote = F, row.names = F, sep = "\t")
