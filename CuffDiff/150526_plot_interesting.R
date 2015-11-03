## this script is to make plots of some interesting genes
#load all the data to query later

load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/allGenes.Rdata")
load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/regulatedGenes.Rdata")
#load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/sigGenes.Rdata")
load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/allGOgeneData.Rdata")
load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/genes.fpkm_table")
library("cummeRbund")
library("dplyr")
library("tidyr")
cuff<-readCufflinks()
std <- function(x) sd(x, na.rm=T)/sqrt(length(x)) ## function for standard error of the mean

########### grep words from gene descriptions or names to find genes of interest
geneNamesBarplot<-function(myGeneNames,all_genes=allGenes,field="name"){
  myGeneNames=genes
  all_genes=regulatedGenes
  field="description"
  # This function performs grep in either the external gene name or description field and plots those genes
  genesOfInterest<-data.frame(matrix(NA, nrow = 1, ncol = 20))
  names(genesOfInterest)<-names(all_genes)
  for(name in myGeneNames){
    if(field=="name"){related<-all_genes[grep(name,all_genes$external_gene_name),]}
    if(field=="description"){related<-all_genes[grep(name,as.character(all_genes$description)),]}
    genesOfInterest<-rbind(genesOfInterest,related)
  }
  genesOfInterest<-genesOfInterest[-1,]
  if(nrow(genesOfInterest)==0){return("sorry, no matches found!")}
  names<-genesOfInterest[,c("gene_id","external_gene_name")]
  geneIDs<-genesOfInterest$gene_id
  cb_genes<-cummeRbund::getGenes(cuff,geneIDs)
  print(expressionBarplot(cb_genes,labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE)+scale_x_discrete(labels=genesOfInterest[order(genesOfInterest$gene_id),]$external_gene_name))+theme(text = element_text(size=12,colour="black"))
  #print(expressionBarplot(cb_genes,labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE))
  #fpkm(cb_genes)
  temp<-fpkm(cb_genes)
  reps<-fpkm(cb_genes,replicates=T)
  temp<-merge(temp,names,by="gene_id",all.x=TRUE)  
  ## if sample name = wt then 4 -> ifelse
  temp$samplesize<-rep(c(3,4))
  temp$sem<-temp$stdev/sqrt(temp$samplesize)
  p=ggplot(data=temp, aes(x=external_gene_name, y=fpkm,fill=sample_name,guide=sample_name)) + 
    geom_bar(stat="identity",position=dodge, colour="black") +
    geom_errorbar(aes(ymin=fpkm-sem, ymax=fpkm+sem),width=0.25, position=dodge)+
    scale_fill_manual(values=c("darkslateblue","red"))+
    theme(text = element_text(size=12,colour="black"), 
          axis.text.x = element_text(angle=90, vjust=1, colour="black"),
          axis.text.y = element_text(colour="black"))
  print(p)
  return(temp)
}
genes=c("glucose","gaba","serotonin","noradrenaline", "dopamine","neuropeptide y","leptin")
genes=c("glutamate","ghrelin","cholecystokinin","neurotensin","vasopressin","oxytocin","glucagon","muscarinic","corticotropin","hypocretin")
genes=c("Hcrt")
out<-geneNamesBarplot(myGeneNames=genes,all_genes=regulatedGenes,field="description")
View(out)

########## the same thing for GO terms that were enriched in the experiment
GOtermBarplot<-function(identifier,go_genes=gene_exp_GO,field="ID"){
  #field="name"
  #go_genes=gene_exp_GO
  #identifier="transcription"
  GOofInterest<-data.frame(matrix(NA, nrow = 1, ncol = ncol(go_genes)))
  names(GOofInterest)<-names(go_genes)
  for(name in identifier){
    if(field=="name"){related<-go_genes[grep(name,go_genes$go_name),]}
    if(field=="ID"){related<-go_genes[grep(name,as.character(go_genes$go_id)),]}
    GOofInterest<-rbind(GOofInterest,related)
  }
  GOofInterest<-GOofInterest[-1,]
  if(nrow(GOofInterest)==0){return("sorry, no matches found!")}
  geneIDs<-GOofInterest$gene_id
  cb_genes<-cummeRbund::getGenes(cuff,geneIDs)
  print(expressionBarplot(cb_genes,labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE)+scale_x_discrete(labels=GOofInterest[order(GOofInterest$gene_id),]$external_gene_name))
  #print(expressionBarplot(cb_genes,labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE))
  #fpkm(cb_genes)
  return(GOofInterest)
}
out<-GOtermBarplot(identifier="locomotory behavior",field="name")
out<-out[,c("external_gene_name","description","log2.fold_change.")]
row.names(out)<-NULL
View(out)



######### other approach: read gene names from GOs, then search for those names. This way you can look for any GO
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
##for manually entering GO IDs of interest:
goID<-"GO:0032868"
genes<- getBM(attributes='external_gene_name', filters = 'go_id', values = goID, mart = mart)$external_gene_name
out<-geneNamesBarplot(myGeneNames=genes,all_genes=regulatedGenes,field="name")
View(out)

library("ggplot2")
dodge <- position_dodge(width=0.9)
geneIsoformsBarplot<-function(name,all_genes=allGenes){
  all_genes=allGenes
  # This function find an exact match for one gene name
  related<-all_genes[which(all_genes$external_gene_name==name),]
  if(nrow(related)==0){return("sorry, no matches found!")}
  geneIDs<-related$gene_id
  cb_genes<-cummeRbund::getGenes(cuff,geneIDs)
  print(expressionBarplot(cb_genes,labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE)+scale_x_discrete(labels=related[order(related$gene_id),]$external_gene_name))
  print(expressionBarplot(isoforms(cb_genes),labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE))
  print(expressionBarplot(TSS(cb_genes),labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE))
  
  temp<-fpkm(cb_genes)
  p=ggplot(data=temp, aes(x=gene_id, y=fpkm,fill=sample_name,guide=sample_name)) + 
    geom_bar(stat="identity",position=dodge, colour="black") +
    geom_errorbar(aes(ymin=fpkm-stdev, ymax=fpkm+stdev),width=0.25, position=dodge)+
    scale_fill_manual(values=c("darkslateblue","red"))+
    ggtitle(name)+
    theme(text = element_text(size=12,colour="black"), 
          axis.text.x = element_text(angle=0, vjust=1, colour="black"),
          axis.text.y = element_text(colour="black"))
  print(p)
  
  
  return(fpkm(cb_genes))
}
name="Mycbp2"
geneIsoformsBarplot(name)
