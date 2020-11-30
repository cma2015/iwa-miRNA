#!/usr/bin/env Rscript
options(stringsAsFactors = F)
library(ggplot2)
library(ggthemes)
library(ggsci)

##
theme_Publication <- function(base_size=14, base_family="sans") {
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(angle=30,vjust =1, hjust = 1),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            # legend.position = "bottom",
            # legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            # legend.margin = unit(0, "cm"),
            # plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

pair_genes <- read.table("input.collinearity", sep = "\t")
pair_genes[,1] <- apply(pair_genes, 1, function(x){unlist(strsplit(gsub(" ", "", x[1]), "-"))[1]})

miRNA_gff <- read.table("miRNAs.gff", stringsAsFactors = F, sep = "\t")

library(igraph)
miRNA_pair <- pair_genes[grepl(":", pair_genes[,2])&grepl(":", pair_genes[,3]), ]

g=graph.data.frame(miRNA_pair[,2:3])
miRNA_graph <- cluster_infomap(g)
names(g[1])
miRNA_name <- sort(unique(miRNA_graph$membership))
miRNA_typelist <- list()
miRNA_typevec <- vector()

mir_table <- read.table("final_table.txt", sep = "\t", header = T, stringsAsFactors = F)
if(!"Extended_stem_loop_loc"%in%colnames(mir_table)){
  colnames(mir_table)[1:2] <- c("ID", "Extended_stem_loop_loc")
  colnames(mir_table)[ncol(mir_table)] <- "Source"
}
rownames(mir_table) <- mir_table$Extended_stem_loop_loc

for(ii in miRNA_name){
  name_index <- miRNA_graph$membership%in%ii
  mi_name <- names(g[1])[name_index]
  mi_name <- mir_table[mi_name, ]$ID
  write.table(paste(paste0("Group",ii), paste0(mi_name, collapse = "\t"), sep = "\t"), "miRNA_duplication_in_blocks.txt", 
              append = T, row.names = F, col.names = F, sep = "\t", quote = F)
}

## duplication type plot
gene_type <- read.table("input.gene_type", sep = "\t")
rownames(gene_type) <- gene_type[,1]
gene_type[gene_type[,2]==0, 2] <- "singleton"
gene_type[gene_type[,2]==1, 2] <- "dispersed"
gene_type[gene_type[,2]==2, 2] <- "proximal"
gene_type[gene_type[,2]==3, 2] <- "tandem"
gene_type[gene_type[,2]==4, 2] <- "WGD/segmental"

type_out <- gene_type[mir_table$Extended_stem_loop_loc,]
type_out[,1] <- mir_table[type_out[,1], ]$ID
colnames(type_out) <- c("ID", "Type")

write.table(type_out, "summary_of_miRNAs_duplication.txt", row.names = F, sep = "\t", quote = F)

tmp_a <- table(gene_type[mir_table[mir_table$Source=="+", ]$Extended_stem_loop_loc, 2])
tmp_b <- table(gene_type[mir_table[mir_table$Source!="+", ]$Extended_stem_loop_loc, 2])
tmp_df <- data.frame("Type" = c(rep("Newly annotated",length(tmp_a)), rep("Alreadly annotated", length(tmp_b))),
                     "Gene_type" = c(names(tmp_a), names(tmp_b)),
                     "Value"=c(tmp_a, tmp_b), stringsAsFactors = F)
tmp_df[,2] <- factor(tmp_df[,2], levels = c("WGD/segmental", "tandem", "proximal", "dispersed", "singleton"))

pdf("Genomic_duplication_and_miRNA_expansion.pdf", 4, 5)
ggplot(data = tmp_df, aes(x=Type, y=Value, fill=Gene_type))+
  geom_bar(stat="identity") +
  ylab("Number of miRNAs") +
  xlab("") +
  theme(axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(colour = "black", size = 0.5),
        legend.title = element_blank())+
  scale_fill_lancet() + theme_Publication()
dev.off()
