### get the significantly regulated genes, as well as the ones that are absent in one condition
gene_exp.diff<-read.table("gene_exp.diff", sep="\t", header=T)
gene_exp.diff<-gene_exp.diff[gene_exp.diff$status=="OK",]
gene_exp.diff<-gene_exp.diff[gene_exp.diff$value_1>1 | gene_exp.diff$value_2>1,]
nrow(gene_exp.diff)
# [1] 25531

sigGene_table<-gene_exp.diff[gene_exp.diff$q_value<0.05,]
ko_only<-gene_exp.diff[gene_exp.diff$value_1==0,]
wt_only<-gene_exp.diff[gene_exp.diff$value_2==0,]
interest<-rbind(sigGene_table,wt_only,ko_only)

### plot a histogram of the log2fold changes
lim.x<-c(-5,5)
p1 <- hist(gene_exp.diff$log2.fold_change,breaks=50) 
p2 <- hist(gene_exp.diff[gene_exp.diff$q_value<0.05,]$log2.fold_change,breaks=50)                  
plot( p1, col=rgb(0,0,1,1/4), xlim=lim.x,main="Log2FC")  
plot( p2, col=rgb(1,0,0,1/4), xlim=lim.x, add=T)  
plot( p1, col=rgb(0,0,1,1/4), xlim=lim.x,ylim=c(0,100),main="Log2FC (zoom)")  
plot( p2, col=rgb(1,0,0,1/4), xlim=lim.x,ylim=c(0,100), add=T)

##### get associated gene names etc
library(biomaRt)
# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name","gene_biotype","description","transcript_length"),
                 filters = "ensembl_transcript_id", values = substr(as.character(gene_exp.diff$gene),1,18),
                 mart = mart)

id2ens<-cbind(as.character(gene_exp.diff$gene_id),substr(as.character(gene_exp.diff$gene),1,18))
colnames(id2ens)<-c("gene_id","ensembl_transcript_id")
id2name<-merge(id2ens,results,by="ensembl_transcript_id",all.x=T)

allGenes<-merge(gene_exp.diff,id2name,by="gene_id",all.x=T)
regulatedGenes<-merge(interest,id2name,by="gene_id",all.x=T)
write.table(allGenes,file="allGenesAnnotated.txt",sep="\t")
write.table(regulatedGenes,file="regulatedGenesAnnotated.txt",sep="\t")
save(allGenes,file="allGenes.Rdata")
save(regulatedGenes,file="regulatedGenes.Rdata")


alg<-factor(rep(0,times=nrow(allGenes)))
alg<-factor(as.numeric(allGenes$gene_id %in% regulatedGenes$gene_id))
names(alg)<-allGenes$external_gene_name
### make GO plots
# make a vector with gene IDs and whether or not they are significantly regulated



library("cellplot")
library("topGO")
library("plyr")

onts = c( "MF", "BP", "CC" )
tab = as.list(onts)
names(tab) = onts
# get data f
for(i in 1:3){
  ## prepare data
  tgd <- new( "topGOdata", ontology=onts[i], allGenes = alg, nodeSize=5,
              annot=annFUN.org, mapping="org.Mm.eg.db", ID = "symbol"  )
  ## run tests
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  ## look at results
  tab[[i]] <- GenTable( tgd, Fisher.elim = resultTopGO.elim, orderBy = "Fisher.elim" , topNodes = 200)
}

topGOResults <- rbind.fill(tab)
write.table(topGOResults, file = "topGOResults.txt", sep="\t")

#For each significant GO term, get the genes associated with this term combine with their expression information.

write.table(NULL,file="GO_gene_data.txt",append=FALSE, sep="\t")
pdf("GOterms.pdf")
for(iter in 1:3){
  goBP<-tab[[iter]]
  gene_exp_GO<-data.frame()
  
  nr_terms<-30
  for(i in c(1:nr_terms)){
    print(goBP$GO.ID[i])
    #gets gene name and go_id for all genes annotated with GO-term
    GO.members<- getBM(attributes=c('external_gene_name', 'go_id', 'name_1006'), filters = 'go_id', values = goBP$GO.ID[i], mart = mart)
    #select the relevant genes and store their DE information
    sel<-allGenes$external_gene_name %in% GO.members$external_gene_name
    GO.gene_data<-allGenes[sel & allGenes$q_value<0.1,]
    GO.gene_data$go_id<-rep(goBP$GO.ID[i],nrow(GO.gene_data))
    GO.gene_data$go_name<-rep(GO.members$name_1006[1],nrow(GO.gene_data))
    GO.gene_data$Fisher.elim<-rep(goBP$Fisher.elim[i],nrow(GO.gene_data))
    gene_exp_GO<-rbind(gene_exp_GO,GO.gene_data)
    write.table(GO.gene_data,file="GO_gene_data.txt",append=TRUE, sep="\t")
  }
  
  # remove entries which have no genes (how did these get there??)
  x<-rep(TRUE,nrow(goBP))
  for(i in c(1:nrow(goBP))){
    sel<-which(gene_exp_GO$go_id==goBP$GO.ID[i])
    if(length(sel)<3){x[i]=FALSE}
  }
  goBP<-goBP[x,]
  if(nr_terms>nrow(goBP)){nr_terms<-as.numeric(nrow(goBP))}
  
  x <- -log2(as.numeric(goBP$Fisher.elim[1:nr_terms]))
  names(x)<-goBP$Term[1:nr_terms]
  cells<-list()
  for(i in c(1:nr_terms)){
    sel<-which(gene_exp_GO$go_id==goBP$GO.ID[i])
    if(length(sel)!=0){cells[[i]]<-gene_exp_GO[sel,]$log2.fold_change.}
  }
  
  cell.plot(x, cells, xlab="-log2(Fisher.elim)", cell.outer=1,xlab.cex=1, main=onts[iter])
  
}
dev.off()

### Get the IDs of all my genes of interest for use in cummeRbund


sigIDs<-c(as.character(sigGene_table$test_id),as.character(wt_only$test_id),as.character(ko_only$test_id))
save(sigIDs,file="sigIDs.Rdata")

GOgeneData<-read.table("GO_gene_data.txt",fill=NA)





