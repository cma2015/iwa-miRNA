##
library(igraph)
options(stringsAsFactors = F)

##
mydata <- read.table("Conserved.txt")
g=graph.data.frame(mydata)
miRNA_graph <- cluster_infomap(g)
names(g[1])
species_name <- names(rev(sort(table(sapply(strsplit(mydata[,1], "-"), 
                                            function(x)x[1])))))[1]

##
novel_miRNA <- read.table("Non_conserved.txt")
novel_mat <- matrix(nrow = length(names(g[1]))+nrow(novel_miRNA), ncol = 1)
rownames(novel_mat) <- c(names(g[1]), novel_miRNA[,1])

miRNA_name <- sort(unique(miRNA_graph$membership))
ll <- vector()
for(x in letters){for(y in letters){ll <- c(ll, paste0(x, y, collapse = ""))}}
rawNum =1
for(ii in miRNA_name){
  name_index <- miRNA_graph$membership%in%ii
  mirna_name <- names(g[1])[name_index]
  cat(mirna_name,"\n")
  know_number <- any(!grepl(":", mirna_name))
  not_all_novel <- any(grepl(":", mirna_name))
  if(know_number&not_all_novel){
    cat("one_known_one_novel", ii, "\n")
    novel_mat[mirna_name,1] <- mirna_name
    novel_name <- mirna_name[grepl(":", mirna_name)]
    len_nov <- length(novel_name)
    known_name <- mirna_name[!grepl(":", mirna_name)]
    known_name <- gsub("-Known", "", known_name)
    tmp_inedx <- grepl("-Novel", known_name)
    known_name[!tmp_inedx] <- gsub("[a-z]$", "", known_name[!tmp_inedx])
    known_name <- names(sort(table(known_name), decreasing = T))[1]
    novel_mat[novel_name,1] <- unique(paste0(known_name, "_N", 1:len_nov))
  } else if(all(grepl(":", mirna_name))){
    cat("two_novel", ii, "\n")
    len_nov <- length(mirna_name)
    novel_mat[mirna_name,1] <- unique(paste0(species_name, 
                                             "-MIR_N", rawNum, 
                                             c(letters,ll)[1:len_nov]))
    rawNum <- rawNum+1
  } else {
    cat("all_known", ii, "\n")
    novel_mat[mirna_name,1] <- mirna_name
  }
}

know_gene <- novel_miRNA[!grepl(":", novel_miRNA[,1]),1]
novel_gene <- novel_miRNA[grepl(":", novel_miRNA[,1]),1]
novel_mat[know_gene,1] <- know_gene
novel_mat[novel_gene,1] <- paste0(species_name, "-MIR_N", 
                                  rawNum:(rawNum+length(novel_gene)-1), "a")

write.table(novel_mat, "NameChange.txt", quote = F, sep = "\t", col.names = F)
