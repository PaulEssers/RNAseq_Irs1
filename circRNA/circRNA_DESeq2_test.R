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

load("circlesExpressed.Rdata")

countdata <- circlesExpressed[,3:9]
colnames(countdata) <- gsub('_circ.txt', "", colnames(countdata))
rownames(countdata)<-paste(circlesExpressed$host_gene,circlesExpressed$location,sep="_")
total_reads<-c(11153164,9010884,9420160,8877172,8576952,9527448,8998416)/4
non_circle<-total_reads - colSums(countdata)
countdata<-rbind(countdata,non_circle)
colSums(countdata)

countdata <- as.matrix(countdata)
head(countdata)

# Assign condition 
condition <- factor(c(rep("irs1ko", 3),rep("wt", 4)))

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
library(DESeq2)
DESeq2Table <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
DESeq2Table$condition<-relevel(DESeq2Table$condition, "wt")

GeneCounts <- counts(DESeq2Table)
colSums(GeneCounts)

DESeq2Table <- estimateSizeFactors(DESeq2Table)
sizeFactors(DESeq2Table)


### QC plots, don't look good, don't know if it matters..
multiecdf( counts(DESeq2Table, normalized = T),xlab="mean counts", xlim=c(0, 1000))
multidensity( counts(DESeq2Table, normalized = T),xlab="mean counts", xlim=c(0, 200))
rld <- rlogTransformation(DESeq2Table, blind=TRUE)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "OrRd"))(100)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
DESeq2::plotPCA(rld)+geom_text(aes(label=colnames(rld)), size=3)

DESeq2Table <- estimateDispersions(DESeq2Table)
plotDispEsts(DESeq2Table)

DESeq2Table <- nbinomWaldTest(DESeq2Table)
DESeq2Res <- results(DESeq2Table, pAdjustMethod = "BH")
DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]
# # use z-scores as input to FDRtool to re-estimate the p-value
FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = T)
DESeq2Res[,"padj"] <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
hist(FDR.DESeq2Res$pval, col = "royalblue4",main = "WT vs Irs1ko, correct null model", xlab = "CORRECTED p-values")
table(DESeq2Res[,"padj"] < 0.1)

View(DESeq2Res[DESeq2Res[,"padj"] < 0.1,])
#also check without the multiple hypothesis testing, to see if the others I identified come up..
View(DESeq2Res[DESeq2Res[,"pvalue"] < 0.01,])

