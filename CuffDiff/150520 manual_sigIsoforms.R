### get the significantly regulated genes, as well as the ones that are absent in one condition
isoform_exp.diff<-read.table("isoform_exp.diff", sep="\t", header=T)
isoform_exp.diff<-isoform_exp.diff[isoform_exp.diff$status=="OK",]
isoform_exp.diff<-isoform_exp.diff[isoform_exp.diff$value_1>1 | isoform_exp.diff$value_2>1,]
nrow(isoform_exp.diff)
# [1] 25531

sigIsof_table<-isoform_exp.diff[isoform_exp.diff$q_value<0.05,]
ko_only<-isoform_exp.diff[isoform_exp.diff$value_1==0,]
wt_only<-isoform_exp.diff[isoform_exp.diff$value_2==0,]
interest<-rbind(sigIsof_table,wt_only,ko_only)

#sigIDs<-c(as.character(sigIsof_table$test_id),as.character(wt_only$test_id),as.character(ko_only$test_id))
save(interest,file="sigIsofs.Rdata")


### plot a histogram of the log2fold changes
lim.x<-c(-8,8)
p1 <- hist(isoform_exp.diff$log2.fold_change,breaks=50,xlim=lim.x) 
p2 <- hist(isoform_exp.diff[isoform_exp.diff$q_value<0.05,]$log2.fold_change,breaks=50)                  
plot( p1, col=rgb(0,0,1,1/4), xlim=lim.x,main="Log2FC")  
plot( p2, col=rgb(1,0,0,1/4), xlim=lim.x, add=T)  
plot( p1, col=rgb(0,0,1,1/4), xlim=lim.x,ylim=c(0,100),main="Log2FC (zoom)")  
plot( p2, col=rgb(1,0,0,1/4), xlim=lim.x,ylim=c(0,100), add=T)

##### get associated gene names etc
library(biomaRt)
# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name","transcript_biotype","description","transcript_length"),
                 filters = "ensembl_transcript_id", values = substr(as.character(isoform_exp.diff$gene),1,18),
                 mart = mart)

id2ens<-cbind(as.character(isoform_exp.diff$gene_id),substr(as.character(isoform_exp.diff$gene),1,18))
colnames(id2ens)<-c("gene_id","ensembl_transcript_id")
id2name<-merge(id2ens,results,by="ensembl_transcript_id",all.x=T)

allIsoforms<-merge(isoform_exp.diff,id2name,by="gene_id",all.x=T)
regulatedIsoforms<-merge(interest,id2name,by="gene_id",all.x=T)

rm(id2ens)
rm(id2name)

write.table(allGenes,file="allIsoformsAnnotated.txt",sep="\t")
write.table(regulatedIsoforms,file="regulatedIsoformsAnnotated.txt",sep="\t")
save(allIsoforms,file="allIsoforms.Rdata")
save(regulatedIsoforms,file="regulatedIsoforms.Rdata")

regulatedIsoforms<-regulatedIsoforms[order(regulatedIsoforms$log2.fold_change.,decreasing=T),]
regulatedIsoformsHigh<-(regulatedIsoforms[which(regulatedIsoforms$value_1>10 | regulatedIsoforms$value_2>10),])
regulatedGenesHigh<-(regulatedGenes[which(regulatedGenes$value_1>10 | regulatedGenes$value_2>10),])

isoformRegulatedGeneNames<-unique(regulatedIsoforms$external_gene_name)
RegulatedGeneNames<-unique(regulatedGenes$external_gene_name)

area1<-as.numeric(length(isoformRegulatedGeneNames))
area2<-as.numeric(length(RegulatedGeneNames))
cross.area<-sum(as.numeric(RegulatedGeneNames %in% isoformRegulatedGeneNames))
draw.pairwise.venn(area1, area2, cross.area, category = c("isoform level","gene_level"),fill=c("blue","red"),main='overlap between gene and transcript level DE')


isoformSpecific<-regulatedIsoformsHigh[-which(regulatedIsoformsHigh$external_gene_name %in% regulatedGenesHigh$external_gene_name),]
isoformSpecificNames<-unique(isoformSpecific$external_gene_name)
isoformSpecificNames<-isoformSpecificNames[-1]
out<-geneNamesBarplot(myGeneNames=isoformSpecificNames,all_genes=allGenes,field="name")
  
geneIDs<-unique(isoformSpecific$gene_id)

pdf("isoform_specific_regulation.pdf")
for(id in geneIDs){
  cb_genes<-cummeRbund::getGenes(cuff,id)
  print(expressionBarplot(isoforms(cb_genes),labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE))
}
dev.off()

isoformSpecific[grep("46168",isoformSpecific$ensembl_transcript_id),]
