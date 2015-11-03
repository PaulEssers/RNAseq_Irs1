####### this script merges multiple output files from CIRCexplorer and does a crude differential expression and sequence retrieval

### merge the files into one dataframe
files <- list.files(pattern="_circ.txt")
circles<-data.frame()
for(i in 1:length(files)){
  data<-read.table(files[i],sep="\t",stringsAsFactors=FALSE)
  data$V1<-paste(data$V1,paste(data$V2,data$V3,sep="-"),sep=":")
  data<-data[,c(1,6,15,13)]
  names(data)<-c("location","strand","host",files[i])
  if(i==1){circles<-data}
  else{circles<-merge(circles,data,by="location",all=T,fill=0)}
}

#rearrange the dataframe, remove NAs
circles<-circles[,c(1,seq(2,ncol(circles),3),seq(3,ncol(circles),3),seq(4,ncol(circles),3))]
circles[is.na(circles)]<-0

#make one column with all the host gene names
condense<-function(names=rep("",10)){
  for(i in c(1:length(names))){
    #print(names[i])
    if(names[i]!=0){return(names[i])}
}}

circles$new_strand<-rep(0,nrow(circles))
circles$host_gene<-rep(0,nrow(circles))
for(i in c(1:nrow(circles))){
  circles$host_gene[i]<-condense(circles[i,c(9:15)])
  circles$new_strand[i]<-condense(circles[i,c(2:8)])
}
safe<-circles
circles<-safe

circles$strand<-circles$new_strand
circles<-circles[,-c(2:7,9:15,23)]
circles$host_gene<-as.character(circles$host_gene)

safe<-circles

# select only the circles which have > 20 counts overall
library("dplyr")
circles<-safe
circles$total_expression<-rep(0,nrow(circles))
#circles$fold_change<-rep(0,nrow(circles))
circles$wt.avg<-rep(0,nrow(circles))
circles$wt.sem<-rep(0,nrow(circles))
circles$irs1ko.avg<-rep(0,nrow(circles))
circles$irs1ko.sem<-rep(0,nrow(circles))

#i=1
#mean(as.numeric(circles[i,c(2:4)])

std <- function(x) sd(x, na.rm=T)/sqrt(length(x)) ## function for standard error of the mean
for(i in c(1:nrow(circles))){
  circles$total_expression[i]<-sum(as.numeric(circles[i,c(3:9)]))
  #circles$fold_change[i]<-mean(as.numeric(circles[i,c(2:4)]))/mean(as.numeric(circles[i,c(5:8)]))
  circles$wt.avg[i]<-mean(as.numeric(circles[i,c(6:9)]))
  circles$wt.sem[i]<-std(as.numeric(circles[i,c(6:9)]))
  circles$irs1ko.avg[i]<-mean(as.numeric(circles[i,c(3:5)]))
  circles$irs1ko.sem[i]<-std(as.numeric(circles[i,c(3:5)]))
}
save(circles,file="allCircles.Rdata")
circlesExpressed<-circles[circles$total_expression>20,]
rm(circles)
rm(safe)
rm(data)
rm(files)
rm(i)
### check which of these are differentially expressed between the samples

circlesExpressed<-circlesExpressed[order(circlesExpressed$total_expression, decreasing=T),]
circlesExpressed$fold_change<-circlesExpressed$irs1ko.avg/circlesExpressed$wt.avg
circlesExpressed$log2FC<-log2(circlesExpressed$fold_change)

for(i in c(1:nrow(circlesExpressed))){
  circlesExpressed$dpois[i]<-dpois(sum(as.numeric(circlesExpressed[i,c(3:5)])), lambda=sum(as.numeric(circlesExpressed[i,c(6:9)])), log = FALSE)
}
circlesExpressed$dpois_corr<-p.adjust(circlesExpressed$dpois, method = "fdr")

#write.table(circlesExpressed,file="circlesExpressed.txt",sep="\t")
#typeof(data.frame(circlesExpressed)) # =list? 

DE<-circlesExpressed[circlesExpressed$dpois_corr<0.05,]
DE<-as.tbl(DE)
DE<-select(DE, host_gene, location, wt.avg, wt.sem, irs1ko.avg, irs1ko.sem,fold_change)

## get into decent format for plotting with ggplot2 and print a file with all circRNAs

library("tidyr")
x<-data.frame(DE)
x$host_gene<-as.character(x$host_gene)
x<-gather(x, genotype, value, c(3:6))
x<-separate(x, genotype, into = c("genotype", "metric"), sep = "\\.")
x<-spread(x, metric, value)
x$host_gene<-factor(x$host_gene,levels=x[order(x$fold_change, decreasing=T),]$host_gene)
x$genotype<-factor(x$genotype,levels=c("wt","irs1ko"))
x<-x[order(x$fold_change),]

library("ggplot2")
dodge <- position_dodge(width=0.9)
pdf("circleRNA_expression.pdf")
for(level in levels(x$host_gene)){
  temp<-x[x$host_gene==level,]
  p=ggplot(data=temp, aes(x=location, y=avg,fill=genotype,guide=genotype)) + 
    geom_bar(stat="identity",position=dodge, colour="black") +
    geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem),width=0.25, position=dodge)+
    scale_fill_manual(values=c("darkslateblue","red"))+
    ggtitle(level)+
    guides(fill=guide_legend(reverse=FALSE))+
    theme(text = element_text(size=12,colour="black"), 
        axis.text.x = element_text(angle=0, vjust=1, colour="black"),
        axis.text.y = element_text(colour="black"))
  print(p)
}
dev.off()

### plot only the highest expressed of these with high log2FC
DEhigh<-circlesExpressed[circlesExpressed$dpois_corr<0.05 & circlesExpressed$total_expression>quantile(circlesExpressed$total_expression,0.8) & abs(circlesExpressed$log2FC)>0.6,]
write.table(DEhigh,file="DEhigh.txt",sep="\t")


x<-DEhigh[,c("host_gene", "location", "strand", "wt.avg", "wt.sem", "irs1ko.avg", "irs1ko.sem","total_expression")]
x$host_gene<-as.character(x$host_gene)
x<-gather(x, genotype, value, c(4:7))
x<-separate(x, genotype, into = c("genotype", "metric"), sep = "\\.")
x<-spread(x, metric, value) ##stopped working for unknown reason
x$host_gene<-factor(x$host_gene,levels=x[order(x$total_expression),]$host_gene)
x$genotype<-factor(x$genotype,levels=c("wt","irs1ko"))
x<-x[order(x$fold_change),]

library("ggplot2")
dodge <- position_dodge(width=0.9)
pdf("circleRNA_expression(high).pdf")
  temp<-x#[x$total_expression<250,]
  p=ggplot(data=temp, aes(x=host_gene, y=avg,fill=genotype,guide=genotype)) + 
    geom_bar(stat="identity",position=dodge, colour="black") +
    geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem),width=0.25, position=dodge)+
    scale_fill_manual(values=c("darkslateblue","red"))+
    ggtitle("circRNA expression")+
    ylab("# junction reads")+
    guides(fill=guide_legend(reverse=FALSE))+
    theme(text = element_text(size=12,colour="black"), 
          axis.text.x = element_text(angle=90, vjust=0, colour="black"),
          axis.text.y = element_text(colour="black"))
  print(p)
dev.off()



###### extract the exon sequences of the circles from their location
library(GenomicFeatures)
library(biomaRt)
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("BSgenome.Mmusculus.UCSC.mm10")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


##convert names to entrez ID, make a table with the name, ID, location and strand of the circle, to extract the sequence with 
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
name2entrez <- getBM(attributes = c("external_gene_name", "entrezgene"), filters = "external_gene_name", values = DEhigh$host_gene, mart = mart)

DEhigh$external_gene_name<-DEhigh$host_gene
myCircles<-merge(DEhigh,name2entrez,by="external_gene_name",all.x=T)
myCircles<-myCircles[,c("external_gene_name","location","strand","entrezgene")]
coordinates<-strsplit(myCircles$location, "[:-]", fixed = FALSE, perl = TRUE, useBytes = FALSE)
coordinates<-t(data.frame(coordinates))
colnames(coordinates)<-c("chromosome","start","end")
myCircles<-cbind(myCircles,coordinates)
row.names(myCircles)<-NULL
myCircles$start<-as.numeric(as.character(myCircles$start))
myCircles$end<-as.numeric(as.character(myCircles$end))
#change location to match the database (0 vs 1 position?)
myCircles$start<-as.numeric(myCircles$start)+1

### find the genes in the USCS database and get the exon sequences
# extract gene from txdb
getGene<-function(myCircles_row){
  mygene<-as.character(myCircles_row$entrezgene)
  mygene_txs <- transcriptsBy(txdb, by="gene")[[mygene]]
  # get the names of the exons of my gene
  mygene_tx_names <- mcols(mygene_txs)$tx_name
  # use these to extract those exons from the complete database
  mygene_exbytx <- exonsBy(txdb, "tx", use.names=TRUE)[mygene_tx_names]
  return(mygene_exbytx)
}

#return the sequence of the two joined exons, switched so the junction is in the middle
getCircleJunction<-function(myCircles_row, gene, junction="[junction]",verbose=FALSE){
  #myCircles_row<-myCircles[1,]
  #gene<-my_transcripts  
  #junction="[junction]"
  # this function gets the sequence of the two exons that make up the detected junction, separated by the string "junction"
  # left and right refer to their orientation in the linear sense here
  if(verbose){print(myCircles_row)}
  start<-myCircles_row$start
  end<-myCircles_row$end
  for(gene_nr in 1:length(gene)){
    # go through all transcripts until I find one with both exons
    left=0
    right=0
    transcript<-gene[[gene_nr]]
    my_exons <- split(transcript,c(1:length(transcript)))  
    #head(my_exons)
    for(i in 1:length(my_exons)){
      ranges<-data.frame(my_exons[[i]]@ranges)
      if(ranges$start==start){
        left<-i
        if(verbose){print("left_exon <")}
      }
      if(verbose){print(ranges)}
      if(ranges$end==end){
        right<-i
        if(verbose){print("> right_exon")}
      }
    }  
    #if I found both exons for any transcript, I'm done
    if(left!=0&right!=0){break}
  }
  circle<-my_exons[c(right,left)]
  circle_seq <- getSeq(Mmusculus, circle)
  #sequence<-paste(toString(circle_seq[[1]]),junction,toString(circle_seq[[2]]))
  print("done!")
  #return(sequence)
  return(circle_seq)
}

getCircleExons<-function(myCircles_row, transcript){
  # this function gets the full sequence of the circle, assuming that it consists of all exons in between
  # !!!!! currently assumes + strand, have to rewrite to take into account the strand (?switch start and end might be enough)
  print(myCircles_row)
  my_exons <- split(transcript,c(1:length(transcript)))  
  started=FALSE
  ended=FALSE
  sel=FALSE*length(my_exons)
  for(i in 1:length(my_exons)){
    ranges<-data.frame(my_exons[[i]]@ranges)
    if(ranges$start==myCircles_row$start){
      started=T
      print("started <")
    }
    print(ranges)
    if(!ended){
      sel[i]<-started
    }
    else{sel[i]<-FALSE}
    if(ranges$end==myCircles_row$end){
      ended=T
      print("> ended")
    }
  }
  circle<-my_exons[as.logical(sel)]
  circle_seq <- getSeq(Mmusculus, circle)
  sequence<-""
  for(i in 1:length(circle_seq)){sequence<-paste0(sequence,toString(circle_seq[[i]]))}
  return(sequence)
}

for(row in 1:nrow(myCircles)){
  if(!is.na(myCircles[row,]$entrezgene)){
    my_transcripts<-getGene(myCircles[row,])
    #length(my_transcripts)
    #myCircles$junction[row]<-getCircleJunction(myCircles[row,],gene=my_transcripts)
    circle_seq<-getCircleJunction(myCircles[row,],gene=my_transcripts)
    myCircles$junction_upstream[row]<-circle_seq[1]
    myCircles$junction_downstream[row]<-circle_seq[2]
  }
}

safe<-myCircles
### make a file for primer3 input
fileConn<-file("primer3input.txt")
writePrimer3<-function(myCircles_row,fileConn){
  writeLines(paste0("SEQUENCE_ID=circ",myCirclesRow$external_gene_name), fileConn)
  
  

close(fileConn)
