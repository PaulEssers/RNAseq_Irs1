####### this script merges multiple output files from CIRCexplorer and does a crude differential expression and sequence retrieval

### merge the files into one dataframe
files <- list.files(pattern="_circ.txt")
circles<-data.frame()
for(i in 1:length(files)){
  data<-read.table(files[i],sep="\t",stringsAsFactors=FALSE)
  data$V1<-paste(data$V1,paste(data$V2,data$V3,sep="-"),sep=":")
  data<-data[,c(1,15,13)]
  names(data)<-c("location","host",files[i])
  if(i==1){circles<-data}
  else{circles<-merge(circles,data,by="location",all=T,fill=0)}
}

#rearrange the dataframe, remove NAs
circles<-circles[,c(1,seq(2,ncol(circles),2),seq(3,ncol(circles),2))]
circles[is.na(circles)]<-0

#make one column with all the host gene names
condense<-function(names=rep("",10)){
  for(i in c(1:length(names))){
    #print(names[i])
    if(names[i]!=0){return(names[i])}
}}

circles$host_gene<-rep(0,nrow(circles))
for(i in c(1:nrow(circles))){
  circles$host_gene[i]<-condense(circles[i,c(2:8)])
}
circles<-circles[,-c(2:8)]
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
  circles$total_expression[i]<-sum(as.numeric(circles[i,c(2:8)]))
  #circles$fold_change[i]<-mean(as.numeric(circles[i,c(2:4)]))/mean(as.numeric(circles[i,c(5:8)]))
  circles$wt.avg[i]<-mean(as.numeric(circles[i,c(5:8)]))
  circles$wt.sem[i]<-std(as.numeric(circles[i,c(5:8)]))
  circles$irs1ko.avg[i]<-mean(as.numeric(circles[i,c(2:4)]))
  circles$irs1ko.sem[i]<-std(as.numeric(circles[i,c(2:4)]))
}
save(circles,file="allCircles.Rdata")
circlesExpressed<-circles[circles$total_expression>20,]
rm(circles)
rm(safe)
rm(data)

### check which of these are differentially expressed between the samples

circlesExpressed<-circlesExpressed[order(circlesExpressed$total_expression, decreasing=T),]
circlesExpressed$fold_change<-circlesExpressed$irs1ko.avg/circlesExpressed$wt.avg
circlesExpressed$log2FC<-log2(circlesExpressed$fold_change)

for(i in c(1:nrow(circlesExpressed))){
  circlesExpressed$dpois[i]<-dpois(sum(as.numeric(circlesExpressed[i,c(2:4)])), lambda=sum(as.numeric(circlesExpressed[i,c(5:8)])), log = FALSE)
}
circlesExpressed$dpois_corr<-circlesExpressed$dpois*nrow(circlesExpressed)
circlesExpressed$dpois_corr2<-p.adjust(circlesExpressed$dpois, method = "fdr")
write.table(circlesExpressed,file="circlesExpressed.txt",sep="\t")

DE<-circlesExpressed[circlesExpressed$dpois_corr2<0.05,]
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
DEhigh<-circlesExpressed[circlesExpressed$dpois_corr2<0.05 & circlesExpressed$total_expression>quantile(circlesExpressed$total_expression,0.8) & abs(circlesExpressed$log2FC)>0.6,]
DEhigh<-as.tbl(DEhigh)
x<-select(DEhigh, host_gene, location, wt.avg, wt.sem, irs1ko.avg, irs1ko.sem,total_expression)

x<-data.frame(x)
x$host_gene<-as.character(x$host_gene)
x<-gather(x, genotype, value, c(3:6))
x<-separate(x, genotype, into = c("genotype", "metric"), sep = "\\.")
x<-spread(x, metric, value)
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
write.table(DEhigh,file="DEhigh.txt",sep="\t")


###### extract the exon sequences of the circles from their location
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(GenomicFeatures)
library(biomaRt)

#convert names to entrez ID 
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
name2entrez <- getBM(attributes = c("external_gene_name", "entrezgene"), filters = "external_gene_name", values = DEhigh$host_gene, mart = mart)


DEhigh$external_gene_name<-DEhigh$host_gene
myCircles<-merge(DEhigh,name2entrez,by="external_gene_name",all.x=T)
myCircles<-myCircles[,c("external_gene_name","location","entrezgene")]
coordinates<-strsplit(myCircles$location, "[:-]", fixed = FALSE, perl = TRUE, useBytes = FALSE)
coordinates<-t(data.frame(coordinates))
colnames(coordinates)<-c("chromosome","start","end")
myCircles<-cbind(myCircles,coordinates)
#change location to match the database
myCircles$start<-as.numeric(myCircles$start)+1
myCircles$end<-as.numeric(myCircles$end)+1

# extract genes from txdb
mygene<-as.character(myCircles[1,]$entrezgene)
mygene_txs <- transcriptsBy(txdb, by="gene")[[mygene]]
test<-cbind(myCircles,t(data.frame(coordinates)))

#get the names of the exons of my gene
mygene_tx_names <- mcols(mygene_txs)$tx_name
#use these to extract those exons from the complete database
mygene_exbytx <- exonsBy(txdb, "tx", use.names=TRUE)[mygene_tx_names]

library("BSgenome.Mmusculus.UCSC.mm10")


elementLengths(mygene_exbytx)
head(mygene_exbytx)
trak2_tx_names

mygene_txs@ranges


### plot the circles on a genome track with their host gene


