load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/allGenes.Rdata")
load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/regulatedGenes.Rdata")

all.ncRNA<-allGenes[allGenes$gene_biotype!="protein_coding",]
regulated.ncRNA<-regulatedGenes[regulatedGenes$gene_biotype!="protein_coding",]

novel<-regulatedGenes[is.na(regulatedGenes$gene_biotype),]
known<-regulatedGenes[regulatedGenes$gene!="-",]
known<-known[known$gene_biotype!="protein_coding",]

#the following types of ncRNAs are present:
  
table(known$gene_biotype)
hist(known$transcript_length, main="transcript lengths", breaks=50)
hist(known$transcript_length, main="transcript lengths", breaks=2000,xlim=c(0,300))

#So a lot of the ncRNAs are smaller than 200nt. How reliable are these?

known<-known[-c(4),]

lncRNA<-known[which(as.integer(known$transcript_length)>200),]
lncRNA<-lncRNA[-1,]#duplicate entry for some reason, this messes up the names on the axis later

library("cummeRbund")
cuff<-readCufflinks()
id<-as.character(lncRNA$gene_id)
cb_genes<-cummeRbund::getGenes(cuff,id)
print(expressionBarplot(cb_genes,labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE)+scale_x_discrete(labels=lncRNA[order(lncRNA$gene_id),]$external_gene_name))

geneNamesBarplot(myGeneNames=id,all_genes=allGenes,field="gene_id",title="")
id

lncRNAfollow_up<-lncRNA[c(1,3,6),]
View(lncRNAfollow_up)
setwd("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffnorm_no_options")
gene.count<-read.table("genes.count_table",header=T)
gene.count<-gene.count[which(gene.count$tracking_id %in% lncRNAfollow_up$test_id),]
View(gene.count)


## novel ncRNAs. The following are multi-exon ncRNAs found to be differentially regulated:
novel_id<-c("XLOC_004849","XLOC_015211","XLOC_024455","XLOC_027797","XLOC_038391","XLOC_042045")
cb_genes<-cummeRbund::getGenes(cuff,novel_id)
print(expressionBarplot(cb_genes,labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE))


### also do some general expression tests
allGenesHigh10<-allGenes[allGenes$value_1>quantile(allGenes$value_1,0.90),]


plot(density(all.novel.over1$value_1),xlim=c(0,10),col="red", main="ncRNA expression density")
lines(density(all.ncRNA$value_1),xlim=c(0,10),col="black")

