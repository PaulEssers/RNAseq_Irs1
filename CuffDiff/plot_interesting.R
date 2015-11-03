## this script is to make plots of some interesting genes
#load all the data to query later

load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/allGenes.Rdata")
load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/regulatedGenes.Rdata")
#load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/sigGenes.Rdata")
load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/allGOgeneData.Rdata")
library("cummeRbund")
cuff<-readCufflinks()


########### grep words from gene descriptions or names to find genes of interest
geneNamesBarplot<-function(myGeneNames,all_genes=allGenes,field="name"){
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
  geneIDs<-genesOfInterest$gene_id
  cb_genes<-cummeRbund::getGenes(cuff,geneIDs)
  print(expressionBarplot(cb_genes,labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE)+scale_x_discrete(labels=genesOfInterest[order(genesOfInterest$gene_id),]$external_gene_name))+theme(text = element_text(size=12,colour="black"))
  #print(expressionBarplot(cb_genes,labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE))
  #fpkm(cb_genes)
  return(genesOfInterest)
}
genes=c("glucose","gaba","serotonin","noradrenaline", "dopamine","neuropeptide y","leptin")
genes=c("glutamate","ghrelin","cholecystokinin","neurotensin","vasopressin","oxytocin","glucagon","muscarinic","corticotropin","hypocretin")
genes=c("dopamine")
out<-geneNamesBarplot(myGeneNames=genes,all_genes=allGenes,field="description")
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
out<-GOtermBarplot(identifier="gonadotropin",field="name")


######### other approach: read gene names from GOs, then search for those names. This way you can look for any GO
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
##for manually entering GO IDs of interest:
goID<-"GO:0032868"
genes<- getBM(attributes='external_gene_name', filters = 'go_id', values = goID, mart = mart)$external_gene_name
out<-geneNamesBarplot(myGeneNames=genes,all_genes=regulatedGenes,field="name")
View(out)


### plot transcript variants of single genes
geneTranscriptsBarplot<-function(name,all_genes=allGenes,field="name"){
  # This function performs grep in either the external gene name or description field and plots those genes
  if(field=="name"){genesOfInterest<-all_genes[which(all_genes$external_gene_name==name),]}
  if(field=="description"){genesOfInterest<-all_genes[grep(name,as.character(all_genes$description)),]}
   
  if(nrow(genesOfInterest)==0){return("sorry, no matches found!")}
  if(nrow(genesOfInterest)>1){return(c("multiple matches found, pick one of the following:",genesOfInterest$external_gene_name))}
  
  geneID<-genesOfInterest$gene_id
  cb_genes<-cummeRbund::getGene(cuff,geneID)
  print(expressionBarplot(isoforms(cb_genes),labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE))
  print(expressionBarplot(TSS(cb_genes),labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE))
  return(fpkm(isoforms(cb_genes)))
}
geneTranscriptsBarplot(name="Anks1b",all_genes=allGenes,field="name")
geneTranscriptsBarplot(name="Sergef",all_genes=allGenes,field="name")
geneTranscriptsBarplot(name="Mycbp2",all_genes=allGenes,field="name")
geneTranscriptsBarplot(name="Ptpn3",all_genes=allGenes,field="name")
