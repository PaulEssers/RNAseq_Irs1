library(BiocStyle)
library(geneplotter)
library(ggplot2)
library(plyr)
library(LSD)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(stringr)
library(topGO)
#library(genefilter)
library(biomaRt)
library(dplyr)
#library(EDASeq)
library(fdrtool)
library(cellplot)

# on the command line:
# nice featureCounts -s 2 -T 20 -t exon -g gene_id -a /beegfs/group_lp/home/PEssers/references/GRCm38/Annotation/gencode/gencode.vM4.annotation.gtf -o all.counts /beegfs/group_lp/home/PEssers/irs1ko_hypothalamus/GRCm38/tophat2/*_accepted_chr.bam

countdata <- read.table("all.rcounts", header=TRUE, row.names=1)
countdata <- countdata[,6:(ncol(countdata))]
 
# Remove path and .bam or .sam from filenames
colnames(countdata) <- gsub('_accepted.sorted.dedup.bam', "", colnames(countdata))
colnames(countdata) <- gsub("X.beegfs.group_lp.home.PEssers.irs1ko_hypothalamus.mm10.tophat2.", "", colnames(countdata))

round(colSums(countdata)/1e6,2)
countdata <- as.matrix(countdata)
head(countdata)
 
# Assign condition 
condition <- factor(c(rep("irs1ko", 3),rep("wt", 4)))
 
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
DESeq2Table <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

DESeq2Table$condition<-relevel(DESeq2Table$condition, "wt")

GeneCounts <- counts(DESeq2Table)
round(colSums(GeneCounts)/1e6,2)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
sum(idx.nz)

nz.counts <- subset(GeneCounts, idx.nz)
sam <- sample(dim(nz.counts)[1], 5)
nz.counts[sam, ]
median(nz.counts)
remove(GeneCounts)

### plot some quality control

#### estimate size factors (normalization, should be around 1 for equal size libraries)
DESeq2Table <- estimateSizeFactors(DESeq2Table)
sizeFactors(DESeq2Table)

pdf("DE_QC.pdf")
multiecdf( counts(DESeq2Table, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))
multidensity( counts(DESeq2Table, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))

## creating heatmaps etc.
rld <- rlogTransformation(DESeq2Table, blind=TRUE)
## create a distance matrix between the samples
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "OrRd"))(100)

heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
DESeq2::plotPCA(rld)+geom_text(aes(label=colnames(rld)), size=3)

# check with outlier removed (wt that clusters with mutants)
outliers <- c("wt20038","wt20042")
DESeq2Table.sub <- DESeq2Table[, !(colnames(DESeq2Table) %in% outliers)]
rld <- rlogTransformation(DESeq2Table.sub, blind=TRUE)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "OrRd"))(100)

heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

DESeq2::plotPCA(rld, intgroup=c("condition"))+geom_text(aes(label=colnames(rld)), size=3)

remove(DESeq2Table.sub)
remove(rld)
remove(distsRL)
remove(mat)
remove(hmcol)

### Differential expression analysis

# remove outlier from the analysis for real
DESeq2Table <- DESeq2Table[, !(colnames(DESeq2Table) %in% outliers)]

DESeq2Table <- estimateDispersions(DESeq2Table)
plotDispEsts(DESeq2Table)

DESeq2Table <- nbinomWaldTest(DESeq2Table)
DESeq2Res <- results(DESeq2Table, pAdjustMethod = "BH")
# number of siginificant DE-genes
table(DESeq2Res$padj < 0.1)

hist(DESeq2Res$pvalue, col = "lavender",main = "WT vs Deletion", xlab = "p-values")

remove(DESeq2Table)

## replace the FDR adjusted p values with another, stricter estimation (fdrtools)

#remove non-tested genes and genes with NA pvals (outliers)
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$padj), ]
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$pvalue), ]

# remove adjsuted pvalues, since we add the fdrtool results later on
# (based on the correct p-values)
DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]
# # use z-scores as input to FDRtool to re-estimate the p-value
FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = T)
# 
# DESeq2Res$param[1, "sd"]
DESeq2Res[,"padj"] <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
# 
hist(FDR.DESeq2Res$pval, col = "royalblue4",main = "WT vs Irs1ko, correct null model", xlab = "CORRECTED p-values")
table(DESeq2Res[,"padj"] < 0.1)
plotMA(DESeq2Res)
plot(attr(DESeq2Res,"filterNumRej"),type="b", xlab="quantiles of 'baseMean'", ylab="number of rejections")


sigGenes<-subset(DESeq2Res, padj < 0.1)
head(sigGenes)

p1 <- hist(DESeq2Res$log2FoldChange,breaks=50) 
p2 <- hist(sigGenes$log2FoldChange,breaks=50)                  
plot( p1, col=rgb(0,0,1,1/4), xlim=c(-2,2),main="Log2FC")  
plot( p2, col=rgb(1,0,0,1/4), xlim=c(-2,2), add=T)  
plot( p1, col=rgb(0,0,1,1/4), xlim=c(-2,2),ylim=c(0,100),main="Log2FC (zoom)")  
plot( p2, col=rgb(1,0,0,1/4), xlim=c(-2,2),ylim=c(0,100), add=T)  

dev.off()

write.table(sigGenes, file="sigGenes.txt",sep="\t")

### start GO term enrichment analysis. Define a background set of genes to improve performance

#find genes with similar expression levels as our sigGenes to use as background for GO
# overallBaseMean <- as.matrix(DESeq2Res[, "baseMean", drop = F])
# backG <- genefinder(overallBaseMean, rownames(sigGenes), 10, method = "manhattan")
# backG <- rownames(overallBaseMean)[as.vector(sapply(backG, function(x)x$indices))]
# total <- backG
# backG <- setdiff(backG,  rownames(sigGenes))
# length(backG)

# plot the density of the DE genes and the background to see how well they match
# multidensity( list( 
#        all= log2(DESeq2Res[,"baseMean"]) ,
#        foreground =log2(DESeq2Res[rownames(sigGenes), "baseMean"]), 
#        background =log2(DESeq2Res[backG, "baseMean"])), 
#      xlab="log2 mean normalized counts", main = "Matching for enrichment analysis")

# alg<-factor(c(rep(0,times=length(backG)),rep(1,times=nrow(sigGenes))))
# names(alg)<-c(backG,rownames(sigGenes))
# table(alg)

rownames(DESeq2Res)<-substr(rownames(DESeq2Res),1,18)
DESeq2Res$DE<-factor(as.numeric(DESeq2Res$padj<0.1))
alg<-DESeq2Res$DE
names(alg)<-rownames(DESeq2Res)
### get GO term enrichment data for each type of ontology: molecular function, biological process and cellular component. Use the "elim" algorithms, which eliminates genes from the list once they have been assigned a GO term, to reduce redundancy.

library("topGO")
onts = c( "MF", "BP", "CC" )
tab = as.list(onts)
names(tab) = onts
# get data f
for(i in 1:3){
  
  ## prepare data
  tgd <- new( "topGOdata", ontology=onts[i], allGenes = alg, nodeSize=5,
              annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl"  )
  
  ## run tests
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  
  ## look at results
  tab[[i]] <- GenTable( tgd, Fisher.elim = resultTopGO.elim, orderBy = "Fisher.elim" , topNodes = 200)
  
  
}
topGOResults <- rbind.fill(tab)
write.table(topGOResults, file = "topGOResults.txt", sep="\t")

# For each significant GO term, get the genes associated with this term combine with their expression information.


ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl") #uses mouse ensembl
pdf("GOterms.pdf")
for(iter in 1:3){
  print(onts[iter])
  goBP<-tab[[iter]]
  
  #loop over this, make a list of GO-terms and associated genes including cuffdiff data
  gene_exp_GO<-data.frame()
  nr_terms<-20
  for(i in c(1:nr_terms)){
    print(goBP$GO.ID[i])
    #gets gene name and go_id for all genes annotated with GO-term
    GO.members<- getBM(attributes=c('ensembl_gene_id', 'go_id', 'name_1006'), filters = 'go_id', values = goBP$GO.ID[i], mart = ensembl)
    #select the relevant genes and store their DE information
    sel<-rownames(DESeq2Res) %in% GO.members$ensembl_gene_id
    GO.gene_data<-data.frame(DESeq2Res[sel & DESeq2Res$padj<0.1,])
    GO.gene_data$go_id<-rep(goBP$GO.ID[i],nrow(GO.gene_data))
    GO.gene_data$Fisher.elim<-rep(goBP$Fisher.elim[i],nrow(GO.gene_data))
    gene_exp_GO<-rbind(gene_exp_GO,GO.gene_data)
    write.table(GO.gene_data,file="GO_gene_data.txt",append=(i!=1), col.names=(i==1))
  }
  print("getting ready to plot")
  # remove entries which have no genes (how did these get there??)
  x<-rep(TRUE,nrow(goBP))
  for(i in c(1:nrow(goBP))){
    sel<-which(gene_exp_GO$go_id==goBP$GO.ID[i])
    if(length(sel)==0){x[i]=FALSE}
  }
  table(x)
  goBP<-goBP[x,]
  if(nr_terms>nrow(goBP)){nr_terms<-as.numeric(nrow(goBP))}
  
  x <- -log2(as.numeric(goBP$Fisher.elim[1:nr_terms]))
  names(x)<-goBP$Term[1:nr_terms]
  cells<-list()
  for(i in c(1:nr_terms)){
    sel<-which(gene_exp_GO$go_id==goBP$GO.ID[i])
    if(length(sel)!=0){cells[[i]]<-gene_exp_GO[sel,]$log2FoldChange}
  }
  
  cell.plot(x, cells, xlab="-log2(Fisher.elim)", cell.outer=1,xlab.cex=1, main=onts[iter])
}
dev.off()

### get ncRNAs from the sigGenes
# annotated the database with "biotype"
sigGenes$ensembl_gene_id<-substr(rownames(sigGenes),1,18)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
biotypes<- getBM(attributes=c('ensembl_gene_id', 'gene_biotype', 'transcript_length'), filters = 'ensembl_gene_id', values = rownames(sigGenes), mart = ensembl)
x<-data.frame(sigGenes)
anSig<-merge(x,biotypes,by="ensembl_gene_id")