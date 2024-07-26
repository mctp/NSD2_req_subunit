#EYoung Outline of R based scripts for ChipSeq data 

library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(ggplot2)
library(Rsamtools)
library(rtracklayer)
library(ORFik)


save.image('Name.Rdata')
load('Name.Rdata')
setwd('/path/sample/bed')


d<-list('sample1','sample2','sample3')
for (x in d){
  assign(paste0('data',x),read.delim(paste0(x,'_macs_peaks.bed'), header=FALSE))  #match output format from macs
  }


 colnames(datasample1) <- c('space','start','end','name', 'score', 'C', 'D','E','F') 
 colnames(datasample2) <- c('space','start','end','name', 'score', 'C', 'D','E','F') 
 colnames(datasample3) <- c('space','start','end','name', 'score', 'C', 'D','E','F') 


grSample1<-toGRanges(datasample1, format='bed', header=T)
grSample1$Namey<-datasample1$name
grSample1$Opera<-datasample1$score
reduceKeepAttr(grSample1, keep.names=TRUE,min.gapwidth=500)

grSample2<-toGRanges(datasample2, format='bed', header=T)
grSample2$Namey<-datasample2$name
grSample2$Opera<-datasample2$score
reduceKeepAttr(grSample2, keep.names=TRUE,min.gapwidth=500)

grSample3<-toGRanges(datasample3, format='bed', header=T)
grSample3$Namey<-datasample3$name
grSample3$Opera<-datasample3$score
reduceKeepAttr(grSample3, keep.names=TRUE,min.gapwidth=500)

filtgrSample1<-grSample1[grSample1$score>60] 


compare1 <- findOverlapsOfPeaks(grSample1, grSample2, maxgap=-1L, minoverlap=0L, ignore.strand=TRUE, connectedPeaks=c('keepAll', 'min', 'merge'))
makeVennDiagram(compare1,c('grSample1', 'grSample2'),cex=1.5,main="Compare1 No Filter")
compare1 <- addMetadata(compare1, colNames='score', FUN=mean)
export(compare1$peaklist[['grSample1///grSample2']], 'Shared_one_two.bed')
export(compare1$peaklist$grSample1, 'Only_one.bed')
export(compare1$peaklist$grSample2, 'Only_two.bed')


files <- list(grSample1, grSample2)
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-1000, 1000), verbose=FALSE)
names(peakAnnoList) <- c('one','two')
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)
write.table(as.data.frame(peakAnnoList$grSample1@anno), sep="\t" )
#bash filter from there

