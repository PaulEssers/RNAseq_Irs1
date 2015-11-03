## this script is to make plots of some interesting genes
#load all the data to query later

load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/allGenes.Rdata")
load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/regulatedGenes.Rdata")
#load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/sigGenes.Rdata")
load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/allGOgeneData.Rdata")
load("C:/Users/PEssers/PowerFolders/Projects/IRS1-hypothalamus/Analysis/150513 CuffDiff (no bias correction)/cuffdiff_no_options/genes_fpkm_table.Rdata")
library("cummeRbund")
library("dplyr")
library("tidyr")
cuff<-readCufflinks()
std <- function(x) sd(x, na.rm=T)/sqrt(length(x)) ## function for standard error of the mean
dodge<-position_dodge(0.9)

########### grep words from gene descriptions or names to find genes of interest
geneNamesBarplot<-function(myGeneNames,all_genes=allGenes,field="name",title=""){
  #myGeneNames=genes
  #all_genes=allGenes
  #field="gene_id"
  #title=""
  # This function performs grep in either the external gene name or description field and plots those genes
  genesOfInterest<-data.frame(matrix(NA, nrow = 1, ncol = 20))
  names(genesOfInterest)<-names(all_genes)
  for(name in myGeneNames){
    if(field=="gene_id"){related<-all_genes[which(all_genes$gene_id==name),]}
    if(field=="exact"){related<-all_genes[which(all_genes$external_gene_name==name),]}
    if(field=="name"){related<-all_genes[grep(name,all_genes$external_gene_name),]}
    if(field=="description"){related<-all_genes[grep(name,as.character(all_genes$description)),]}
    genesOfInterest<-rbind(genesOfInterest,related)
  }
  genesOfInterest<-genesOfInterest[-1,]
  if(nrow(genesOfInterest)==0){return("sorry, no matches found!")}
  names<-genesOfInterest[,c("gene_id","external_gene_name")]
  if(field=="gene_id"){names$external_gene_name<-names$gene_id}
  geneIDs<-genesOfInterest$gene_id
  cb_genes<-cummeRbund::getGenes(cuff,geneIDs)
  #print(expressionBarplot(cb_genes,labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE)+scale_x_discrete(labels=genesOfInterest[order(genesOfInterest$gene_id),]$external_gene_name))+theme(text = element_text(size=12,colour="black"))
  #print(expressionBarplot(cb_genes,labRow=T, replicates=T,fullnames=FALSE,logMode=FALSE))
  #fpkm(cb_genes)
  temp<-genes.fpkm_table[which(genes.fpkm_table$gene_id %in% names$gene_id),]
  temp<-merge(temp,names,by="gene_id",all.x=TRUE)
  temps<-temp
  temps<-gather(temps,sample,fpkm,2:8)
  temps<-separate(temps,sample,into=c("genotype","sample"),sep="_")
  temps<-as.tbl(temps)
  temps<-group_by(temps, external_gene_name,genotype)
  temps<-summarize(temps, avg=mean(fpkm),sem=std(fpkm))
  temps$genotype<-factor(temps$genotype,c("wt","irs1ko"))
  if(field=="exact"){temps$external_gene_name<-factor(temps$external_gene_name,levels=myGeneNames)}
  p=ggplot(data=temps, aes(x=external_gene_name, y=avg,fill=genotype,guide=genotype)) + 
    geom_bar(stat="identity",position=dodge, colour="black") +
    geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem),width=0.25, position=dodge)+
    scale_fill_manual(values=c("darkslateblue","red"))+
    ggtitle(title)+
    ylab("fpkm")+
    xlab("")+
    theme(text = element_text(size=18,colour="black"), 
          axis.text.x = element_text(angle=90, vjust=1, colour="black"),
          axis.text.y = element_text(colour="black"))
  print(p)
  print(genesOfInterest)
  return(temp)
}


genes=c("glucose","gaba","serotonin","noradrenaline", "dopamine","neuropeptide y","leptin")
genes=c("glutamate","ghrelin","cholecystokinin","neurotensin","vasopressin","oxytocin","glucagon","muscarinic","corticotropin","hypocretin")
genes=c("Pcsk1","Pcsk2","Cpe","Pomc","Cartpt","Agrp","Npy","Pmch","Mchr1")
genes=c("Cga",)#,"Cpe","Pomc","Cartpt","Agrp","Npy","Pmch","Mchr)
genes=c("Egr1","Egr2","Egr3","Egr4")
genes=c("RP23-322M7.1","Gm11767","Gm11549")
genes=c("XLOC_004849","XLOC_015211","XLOC_024455","XLOC_027797","XLOC_038391","XLOC_042045")
genes=c("Hcrt","Pmch","Mchr1","Drd1a","Drd2","Ddc","Slc6a3")
genes=c("Crh","Oxt","Avp","Pmch","Mc4r")
genes=c("Drd1a","Drd2","Th","Cckbr","Htr2a","Oprd1","Rgs9","Rgs14","Adora2a","Adra1b","Rasd2" )
genes=c("Mao")
genes=c("monoamine")
geneNamesBarplot(myGeneNames=genes,all_genes=regulatedGenes,field="name",title="")
View(out)

View(gene_exp_GO[grep("Vip",gene_exp_GO$external_gene_name),])

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

identifier="visual learning"
out<-GOtermBarplot(identifier=identifier,field="name")
geneNamesBarplot(myGeneNames=out$external_gene_name,all_genes=allGenes,field="exact",title=identifier)
out<-out[,c("external_gene_name","description","log2.fold_change.")]
row.names(out)<-NULL
View(out)

pdf("barplots of GO term genes")


######### other approach: read gene names from GOs, then search for those names. This way you can look for any GO
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
##for manually entering GO IDs of interest:
goID<-"GO:0030431"
genes<- getBM(attributes='external_gene_name', filters = 'go_id', values = goID, mart = mart)$external_gene_name
out<-geneNamesBarplot(myGeneNames=genes,all_genes=regulatedGenes,field="name",title="Sleep")
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
