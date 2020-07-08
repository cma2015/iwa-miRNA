library(seqinr)
library(tibble)

args <- commandArgs(trailingOnly = T)
cur_path <- args[1]
database_name <- list.files(path = args[1])
species_path <- args[2]

raw_data <- matrix(ncol = 12)
colnames(raw_data) <- 1:12
da_vec <- c("miRBase", "PmiREN", "sRNAanno", "Psgenes", "Novel")
da_abr <- c("1","2","3","4", "+")

database_list = vector()
if ("miRBase.txt" %in% database_name){
  tmp_data <- read.table(paste0(cur_path, "/miRBase.txt"), header = T)[,c(1:4,6:11,5)]
  tmp_data <- cbind(tmp_data, "miRBase")
  cat(ncol(tmp_data), "miRBase\n")
    colnames(tmp_data) <- colnames(raw_data)
  raw_data <- rbind(raw_data, tmp_data)
  database_list <- c(database_list, "miRBase")
}
if  ("PmiREN.txt" %in% database_name) {
  tmp_data <- read.table(paste0(cur_path, "/PmiREN.txt"), header = T)[,c(1:4,6:11,5)]
  tmp_data <- cbind(tmp_data, "PmiREN")
  cat(ncol(tmp_data), "PmiREN\n")
  colnames(tmp_data) <- colnames(raw_data)
  raw_data <- rbind(raw_data, tmp_data) 
  database_list <- c(database_list, "PmiREN")
}
if  ("sRNAanno.txt" %in% database_name) {
  tmp_data <- read.table(paste0(cur_path, "/sRNAanno.txt"), header = T)[,c(1:4,6:11,5)]
  tmp_data <- cbind(tmp_data, "sRNAanno")
  cat(ncol(tmp_data), "sRNAanno\n")
  colnames(tmp_data) <- colnames(raw_data)
  raw_data <- rbind(raw_data, tmp_data)
  database_list <- c(database_list, "sRNAanno")
}
if  ("PlantsmallRNAgenes.txt" %in% database_name) {
  tmp_data <- read.table(paste0(cur_path, "/PlantsmallRNAgenes.txt"), header = T)[,c(1:4,6:11,5)]
  tmp_data <- cbind(tmp_data, "Psgenes")
  cat(ncol(tmp_data), "Psgenes\n")
  colnames(tmp_data) <- colnames(raw_data)
  raw_data <- rbind(raw_data, tmp_data)
  database_list <- c(database_list, "Psgenes")
}

if  ("prediction.txt" %in% database_name) {
  tmp_data <- read.table(paste0(cur_path, "/prediction.txt"), header = T)[,c(1:7,9:11,13)]
  tmp_data <- cbind(tmp_data, "Novel")
  cat(ncol(tmp_data), "Novel\n")
  colnames(tmp_data) <- colnames(raw_data)
  raw_data <- rbind(raw_data, tmp_data)
  database_list <- c(database_list, "Novel")
}

rm_dup <- function(loctest, loclist){
  input_loc <- unlist(strsplit(loctest, ":"))
  loc_line <- do.call("rbind",strsplit(loclist, ":"))
  tmp_index <- loc_line[,1]%in%input_loc[1]
  if(sum(tmp_index) == 1){
    return_res <- TRUE
  } else {
    con1 <- as.numeric(input_loc[2])>=as.numeric(loc_line[,2]) & as.numeric(input_loc[2])<=as.numeric(loc_line[,3])
    con2 <- as.numeric(input_loc[3])>=as.numeric(loc_line[,2]) & as.numeric(input_loc[3])<=as.numeric(loc_line[,3])
    con3 <- as.numeric(input_loc[2])<=as.numeric(loc_line[,2]) & as.numeric(input_loc[3])>=as.numeric(loc_line[,3])
    all_related <- which(con1|con2|con3)
    order_index <- which(loc_line[all_related,2]%in%input_loc[2] & loc_line[all_related,3]%in%input_loc[3])
    if(order_index==1){
      return_res <- TRUE
    }else{
      return_res <- FALSE
    }
  }
  return_res
}

like_loc <- function(loctest, loclist){
  input_loc <- unlist(strsplit(loctest, ":"))
  loc_line <- do.call("rbind",strsplit(loclist, ":"))
  chr <- loc_line[,1]%in%input_loc[1]
  con1 <- as.numeric(input_loc[2])>=as.numeric(loc_line[,2]) & as.numeric(input_loc[2])<=as.numeric(loc_line[,3])
  con2 <- as.numeric(input_loc[3])>=as.numeric(loc_line[,2]) & as.numeric(input_loc[3])<=as.numeric(loc_line[,3])
  con3 <- as.numeric(input_loc[2])<=as.numeric(loc_line[,2]) & as.numeric(input_loc[3])>=as.numeric(loc_line[,3])
  return_res <- which(chr&(con1|con2|con3))
  return_res
}

raw_data <- raw_data[-1,]
raw_data <- raw_data[raw_data[,5]!="-"|raw_data[,8]!="-", ]
raw_data[,2] <- gsub("[cC]hr", "", raw_data[,2])
raw_data[,5] <- gsub("[cC]hr", "", raw_data[,5])
raw_data[,8] <- gsub("[cC]hr", "", raw_data[,8])
uniq_mature <- unique(raw_data[,5])
mature_index <- vector()
for( i in uniq_mature){
  mature_index <- c(mature_index, rm_dup(loctest = i, loclist = uniq_mature))
}

loc_number <- uniq_mature[mature_index]
final_mat <- matrix(ncol = 12 + length(database_list), nrow = length(loc_number))
colnames(final_mat) <- c("Precursors", "pLoc", "pSeq", "pLen", 
                         "Loc5p", "Seq5p", "Len5p", "Loc3p", 
                         "Seq3p", "Len3p", 'Mature_arm',database_list, "db_num")
rownames(final_mat) <- loc_number
final_mat[,12:(ncol(final_mat)-1)] <- "<i class='incorrect'></i>"

for(i in loc_number){
  mat_loc <- like_loc(i, raw_data[,5])
  mirline <- grep("MIR", raw_data[mat_loc,1])
  if(length(mirline)>0){
    final_line <- mat_loc[mirline][1]
  } else {
    final_line <- mat_loc[which(raw_data[mat_loc,5]%in%names(table(raw_data[mat_loc,5])[1]))][1]
  }
  final_mat[i, 1:10] <- as.character(raw_data[final_line,1:10])
  tmp_arm <- table(raw_data[mat_loc,11])
  final_mat[i, 11] <- names(sort(tmp_arm)[length(tmp_arm)])
  final_mat[i, colnames(final_mat)%in%raw_data[mat_loc,12]] <- "<i class='correct'></i>"
  # final_mat[i, ncol(final_mat)-1] <- sum(colnames(final_mat)%in%raw_data[mat_loc,11])
  final_mat[i, ncol(final_mat)] <- paste0(da_abr[which(da_vec%in%raw_data[mat_loc,12])], collapse = "")
}

final_mat[,3] <- gsub("U", "T", final_mat[,3])
final_mat[,6] <- gsub("U", "T", final_mat[,6])
final_mat[,9] <- gsub("U", "T", final_mat[,9])

new_loc <- t(apply(final_mat, 1, function(x){
  aa <- unlist(strsplit(x[5], ":"))
  bb <- unlist(strsplit(x[8], ":"))
  tmp_min <- min(as.numeric(c(aa[2],aa[3],bb[2],bb[3])))-20
  tmp_max <- max(as.numeric(c(aa[2],aa[3],bb[2],bb[3])))+20
  tmp_name <-paste0(c(aa[1], tmp_min, tmp_max, aa[4]), collapse = ":")
  out_tmp <- c(aa[1], tmp_min-1, tmp_max, tmp_name, ".", aa[4])
  out_tmp
}))

write.table(new_loc, paste0(cur_path, "/tmp.bed"), quote = F, col.names = F, row.names = F, sep = "\t")

aa <- system(paste0("bedtools getfasta -s -name -fi ", species_path, "/Genome/Genome.fa",
              " -bed ", cur_path, "/tmp.bed -fo ", cur_path, "/tmp_seq.fa"), intern = T)

seq_input <- read.fasta(paste0(cur_path, "/tmp_seq.fa"), as.string = T)
seq_input <- unlist(seq_input)
seq_input <- toupper(seq_input)

tmp_name <- gsub("\\(.+", "", names(seq_input))
final_index <- which(new_loc[,4]%in%tmp_name)
final_mat <- final_mat[final_index, ]
final_mat[,2] <- tmp_name
final_mat[,3] <- seq_input
final_mat[,4] <- nchar(seq_input)

final_mat <- final_mat[final_mat[,11]%in%c("5p","3p","-"), ]
write.table(final_mat, paste0(cur_path, "/Translate_result.txt"), quote = F, sep = "\t", row.names = F)
  