library(ggplot2)
library(rtracklayer)
options(stringsAsFactors = F)

args = commandArgs(trailingOnly=TRUE)

barplotForMyself <- function(dataset, fillcol, yname, perchange, title){
  minVal <- min(dataset$count)
  maxVal <- max(dataset$count)
  p <- ggplot(data = dataset, aes(x = type, y=count))+geom_bar(stat = "identity", fill=fillcol)+
    theme_classic(base_size = 16)+ theme(axis.text = element_text(colour = "black"), 
                                         axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(x="", y=yname, title = title)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 6),limits=c(0,maxVal*1.1))
  if(is.na(perchange[1])){
    p <- p + geom_text(mapping = aes(label= count), size = 5, colour = 'black', vjust = -0.5)
  } else {
    p <- p + geom_text(mapping = aes(label=perchange), size = 5, colour = 'black', vjust = -0)
  }
  p
} 

bedread <- read.table("GFF3/Annotation.bed", sep = "\t", stringsAsFactors = F)
featurecount <- as.data.frame(table(bedread[,5]))
featurecount <- featurecount[order(featurecount$Freq ,decreasing = T), ]

featurecount[,1] <- factor(featurecount[,1], levels = featurecount[,1])

pdf(args[1], width = 8, height = 8)
par(mfrow=c(2,1))

maxVal <- max(featurecount$Freq)
ggplot(data = featurecount, aes(x = Var1, y=Freq))+geom_bar(stat = "identity", fill="#C0504D")+
  theme_classic(base_size = 16)+ theme(axis.text = element_text(colour = "black"), 
                                       axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="", y="Count", title = "The number of different annotation categories")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6),limits=c(0,maxVal*1.1)) + 
  geom_text(mapping = aes(label= Freq), size = 5, colour = 'black', vjust = -0.5)


for(typ in featurecount$Var1){
  gIndexbed <- bedread[bedread$V5%in%typ, ]
  tmpgra <- GRanges(seqnames = gIndexbed$V1,
          ranges = IRanges(start = gIndexbed$V2,
                           end = gIndexbed$V3,
                           names = gIndexbed$V4))
  featurecount[featurecount$Var1==typ, 2] <- round(sum(width(reduce(tmpgra)))/1000000, digits = 2)
}

genomeSize <- read.table("Genome/Genome.fa.fai", sep = "\t", stringsAsFactors = F)
genomeSize <- sum(genomeSize[,2])/1000000

maxVal <- max(featurecount$Freq)
inputLP <- sprintf("%sMb\n(%s%s)" ,featurecount$Freq,round(featurecount$Freq/genomeSize*100, digits = 2),"%")

ggplot(data = featurecount, aes(x = Var1, y=Freq))+geom_bar(stat = "identity", fill="#C0504D")+
  theme_classic(base_size = 16)+ theme(axis.text = element_text(colour = "black"), 
                                       axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="", y="Length (Mb)", title = "The length of different annotation categories")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6),limits=c(0,maxVal*1.1)) + 
  geom_text(mapping = aes(label= inputLP), size = 5, colour = 'black', vjust = -0)

dev.off()

# dataInput <- data.frame("type"=gffName,
#                         "count"=numberCount)
# dataInput[,1] <- factor(dataInput[,1], levels = dataInput[,1])
# barplotForMyself(dataset = dataInput, fillcol="#C0504D", title = "The number of different regions", 
#                                    yname="Counts", perchange = NA)
# 
# genomeSize <- sum(width(makeGRangesFromDataFrame(gffInfo[gffInfo$type%in%"chromosome", ])))/1000000
# lenPer <- data.frame("type"=gffName[-2],
#                         "count"=lenCount[-2])
# lenPer[,1] <- factor(lenPer[,1], levels = lenPer[,1])
# inputLP <- sprintf("%sMb\n(%s%s)" ,lenCount[-2],round(lenCount[-2]/genomeSize*100, digits = 2),"%")
# barplotForMyself(dataset = lenPer, fillcol="#4F81BD", title = "The length and percentage of different regions",
#                                    yname="Length and Percentage", 
#                                    perchange=inputLP)
