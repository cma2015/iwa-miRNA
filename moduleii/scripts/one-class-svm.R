library(e1071)
options(stringsAsFactors = F)
args <- commandArgs(T)

# args <- c("../tmp/OnX4IRTGlT/00feature_out.txt","../tmp/OnX4IRTGlT/out_pool_merge.txt", 0.3 ,"../ML_result.txt")

feature.dir <- args[1] # directory of feature matrix, the first columns are sample names
pos.dir <- args[2] # directory of positive samples
nu <- args[3] # parameter nu, default is 0.3
out_res <- args[4]

featureMat <- as.matrix(read.table(file = feature.dir, sep = "\t", row.names = 1, header = T))
featureMat[1:5,1:5]
colnames(featureMat) <- paste0("feature-", 1:ncol(featureMat))
class(featureMat) <- "numeric"

positivesDf <- read.table(file = pos.dir,  sep = "\t", header = T)
positives <- positivesDf[positivesDf$Source!="+",]$Extended_stem_loop_loc
# unlabels <- setdiff(rownames(featureMat), positives)
# scale 
# featureMatScale <- apply(featureMat, 2, scaleFeature)

posMat <- featureMat[positives, ] #featureMatScale[positives, ]
sdRes <- apply(posMat, 2, sd)
idx <- which(sdRes != 0)
length(idx)
posMat <- posMat[,idx]
dim(posMat)

predictMat <- featureMat[!rownames(featureMat)%in%positives, idx]
dim(predictMat)

svm.model <- svm(posMat,
                 y = NULL,
                 type = 'one-classification',
                 nu = nu,
                 kernel = "radial",
                 scale = TRUE)

positive.label <- predict(svm.model, posMat)
table(positive.label)

# sum(names(predict.label[predict.label])%in%rownames(pc_res)[pc_res[,2]])
# sum(names(predict.label[!predict.label])%in%rownames(pc_res)[pc_res[,2]])
# a <- sum(positive.label)
# b <- sum(!positive.label)
# c <- sum(names(positive.label[positive.label])%in%rownames(pc_res)[pc_res[,2]])
# d <- sum(names(positive.label[!positive.label])%in%rownames(pc_res)[pc_res[,2]])
# x <- matrix(c(a, c, b, d), ncol = 2)
# cat("nu: ", nu, " p-value:", chisq.test(x)$p.value, "\n")

class_result1 <- rep("removed_positive", length(positive.label))
class_result1[positive.label] <- "remained_positive"

if(sum(positivesDf$Source=="+")>0){
  predict.label <- predict(svm.model, predictMat)
  table(predict.label)
  class_result2 <- rep("others", length(predict.label))
  class_result2[predict.label] <- "novel_prediction"
  class_result <- c(class_result1, class_result2)
  names(class_result) <- c(names(positive.label), names(predict.label))
} else{
  class_result <- class_result1
  names(class_result) <- names(positive.label)
}

write.table(as.data.frame(class_result), out_res, quote = F, sep = "\t", col.names = F)
