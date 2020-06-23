args <- commandArgs(trailingOnly = T)

options(stringsAsFactors = F)
tmp_data1 <- read.table(args[1], header = T, sep = "\t")
tmp_data1 <- cbind(tmp_data1, "+")
tmp_data2 <- read.table(args[2], header = T, sep="\t")
tmp_data2 <- tmp_data2[c(1:11,ncol(tmp_data2))]
colnames(tmp_data2) <- colnames(tmp_data1)
raw_data <- rbind(tmp_data1, tmp_data2)

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
  return_res <- which(chr&con1|con2|con3)
  return_res
}

raw_data[,2] <- gsub("[cC]hr", "", raw_data[,2])
raw_data[,5] <- gsub("[cC]hr", "", raw_data[,5])
raw_data[,8] <- gsub("[cC]hr", "", raw_data[,8])
uniq_mature <- unique(raw_data[,5])
mature_index <- vector()
for( i in uniq_mature){
  mature_index <- c(mature_index, rm_dup(loctest = i, loclist = uniq_mature))
}

loc_number <- uniq_mature[mature_index]
final_mat <- matrix(ncol = 12, nrow = length(loc_number))
colnames(final_mat) <- c("Pre-miRNAs", "pLoc", "pSeq", "pLen", 
                         "Loc5p", "Seq5p", "Len5p", "Loc3p", 
                         "Seq3p", "Len3p", "Mature_arm", "Count")
rownames(final_mat) <- loc_number

rearr_db <- function(ii){
  atmp <- unlist(strsplit(ii, ""))
  if("+" %in% atmp){
    zatmp <- atmp[!atmp%in%"+"]
    btmp <- paste0(c(sort(zatmp), "+"), collapse = "")
  }else{
    btmp <- paste0(sort(atmp), collapse = "")
  }
  return(btmp)
}

for(i in loc_number){
  mat_loc <- like_loc(i, raw_data[,5])
  occur_count <- rearr_db(ii = raw_data[mat_loc,12])
  mirline <- which(!grepl(":",raw_data[mat_loc,1]))
  if(length(mirline)>0){
    tmp_line <- mat_loc[mirline]
    tmp_arm <- table(raw_data[tmp_line,11])
    final_arm <- names(sort(tmp_arm)[length(tmp_arm)])      
    final_line <- tmp_line[1]
    final_mat[i, ] <- c(as.character(raw_data[final_line,1:10]), final_arm, occur_count)
    } else {
    final_line <- mat_loc[which(raw_data[mat_loc,5]%in%names(table(raw_data[mat_loc,5])[1]))][1]
    final_mat[i, ] <- c(as.character(raw_data[final_line,1:11]), occur_count)
  }
 
}

final_mat <- final_mat[order(final_mat[,1], decreasing = T), ]
write.table(final_mat, args[3], quote = F, sep = "\t", row.names = F)
