library(Sushi)
install.packages("Sushi")


##bigwigs files
covC <- read.table("/path/NSD2_EpiCrispr/bigwigs/split/chr10.SN7.AR.Ctrl.bw.bedGraph")
covKO <- read.table("/path/NSD2_EpiCrispr/bigwigs/split/chr10.SN14.AR.KO.bw.bedGraph") 
covH3<-read.table("/path//NSD2_EpiCrispr/bigwigs/split/SN3.Ctrl.H3K27Ac.bedGraph")
covH3KO<-read.table("/path/NSD2_EpiCrispr/bigwigs/split/SN11.H3K27Ac.bedGraph")

##table from motifs (assigned row per motif outside here)
bedARKOmo<-read.table(file = '/path//NSD2_EpiCrispr/IGV_Open/AR_KO_knownMatchedSushi.bed', sep = '\t', header = TRUE)
bedARCmo<-read.table(file = '/path/NSD2_EpiCrispr/IGV_Open/AR_Ctrl_knownMatchedSushi.bed', sep = '\t', header = TRUE)
bedARKOmo<-read.table(file = '/path/NSD2_EpiCrispr/IGV_Open/SplitMotifs/KO_Combo.txt', sep = '\t', header = TRUE)
bedARCmo<-read.table(file = '/path/NSD2_EpiCrispr/IGV_Open/SplitMotifs/CT_Combo.txt', sep = '\t', header = TRUE)


# chr10, IF SWITCHING CHR DOUBLE CHECK
geneSet10<-read.delim("/path/chr10.hg38.refGene2.bed", sep="\t",header=TRUE) #From ucsc goldenPath hg38
bed.sub <- geneSet10[which(geneSet10[,"chrom"] == chrom & geneSet10[,"start"] >= chromstart & geneSet10[,"stop"] <= chromend),]
write.table(bed.sub, file = "/path/sub_bed.txt", quote = F, row.names = F)
subbed<-read.delim("/path//sub_bed.txt", sep=" ", header=T)

chrom="chr10"
chromstart=71475541
chromend=71650000
paste(chrom,":",chromstart,"-",chromend)
fun_PlotMotifs(chrom, chromstart,chromend)

##Loop over all coordinates in a bed file (*Change the function to write to pdf first*)
x<- read.table("/mctp/share/users/eleanoyo/TestCoord.bed")
for (i in 1:nrow(x)) {
  fun_PlotMotifs(x[i,1],x[i,2],x[i,3])
}



############################
fun_PlotMotifs <- function(chrom, chromstart, chromend) {
  plot.new()
  ##For saving uncomment one of these
  pdf(file = paste0("/path/",chrom,"_",chromstart,".pdf"), width = 7, height = 10) #7,10 usu
  ##png(file = paste0("/path/",chrom,"_",chromstart,".png"), width = 480, height = 600)
  
  #Margins on graph 
  par(mfrow=c(6,1),mar=c(1,4,1,1)+1) ##this one
  
  ##You can remove this segment if you just want motifs. 
  ##Bigwig plot overlapping
  plotBedgraph(covC,chrom,chromstart,chromend,
               color=SushiColors(2)(2)[1],flip=FALSE,rescaleoverlay=FALSE,range=c(0,300))
  axis(side=2,las=2,tcl=.2,at=pretty(par("yaxp")[c(1,2)]),
       labels=1*pretty(par("yaxp")[c(1,2)]))
  legend("topright",inset=0.025,legend=c("AR Control","AR KO"),
         fill=SushiColors(2)(2),border=SushiColors(2)(2),text.font=1,cex=1.0)
  labelgenome(chrom,chromstart,chromend,n=3,scale="Mb")
  mtext("Peaks (bigwig)",side=2,line=1.75,cex=1,font=1)
  plotBedgraph(covKO,chrom,chromstart,chromend,
               color=SushiColors(2)(2)[2],overlay=FALSE,rescaleoverlay=FALSE,range=c(0,300),flip=FALSE)
  labelgenome(chrom,chromstart,chromend,n=3,scale="Mb")
  axis(side=2,las=2,tcl=.2,at=pretty(par("yaxp")[c(1,2)]),
       labels=1*pretty(par("yaxp")[c(1,2)]))
  mtext("Peaks (bigwig)",side=2,line=1.75,cex=0.5,font=1)
  
  #####Added
  plotBedgraph(covH3,chrom,chromstart,chromend,
               color=SushiColors(4)(4)[2],overlay=FALSE,rescaleoverlay=FALSE,range=c(0,300))
  labelgenome(chrom,chromstart,chromend,n=3,scale="Mb")
  axis(side=2,las=2,tcl=.2,at=pretty(par("yaxp")[c(1,2)]),
       labels=1*pretty(par("yaxp")[c(1,2)]))
  mtext("H3K27Ac Ctrl",side=1,line=1.75,cex=1,font=1)
  legend("topright",inset=0.025,legend=c("H3K27Ac (SN3) Control","H3K27Ac (SN10) KO"),
         fill=c(SushiColors(4)(4)[2],SushiColors(4)(4)[4]),border==c(SushiColors(4)(4)[2],SushiColors(4)(4)[4]),text.font=2,cex=1.0)
  plotBedgraph(covH3KO,chrom,chromstart,chromend
               ,color=SushiColors(4)(4)[4],overlay=FALSE,rescaleoverlay=FALSE,range=c(0,300))
  labelgenome(chrom,chromstart,chromend,n=3,scale="Mb")
  axis(side=2,las=2,tcl=.2,at=pretty(par("yaxp")[c(1,2)]),
       labels=1*pretty(par("yaxp")[c(1,2)]))
  mtext("H3K27Ac KO",side=1,line=1.75,cex=1,font=1)
  ##Motifs
 
 #side note wrote this is R >4.0.0, running in R 3.6.0 it seems to not recognize color command
  plotBed(beddata = bedARCmo,chrom = chrom,
          chromstart = chromstart,chromend = chromend,
          rownumber = bedARCmo$row, type = "circles",
          color=bedARCmo$color,row="given",
          plotbg="grey95",rowlabels=c("ARE","AR_half","FOXA1:AR", "FOX"),
          rowlabelcol=unique(bedARCmo$color),rowlabelcex=0.75),
          labelgenome(chrom,chromstart,chromend,n=3,scale="Mb")
  mtext("AR Ctrl known motifs",side=3, adj=-0.065,line=0.5,font=1,col='blue')
  
  plotBed(beddata = bedARKOmo,chrom = chrom,
          chromstart = chromstart,chromend = chromend,
          rownumber = bedARKOmo$row, type = "circles",
          color=bedARKOmo$color,row="given",
          plotbg="grey95",rowlabels=c("ARE","AR_half","FOXA1:AR", "FOX"),
          rowlabelcol=unique(bedARKOmo$color),rowlabelcex=0.75),
          labelgenome(chrom,chromstart,chromend,n=3,scale="Mb")
  mtext("AR KO known motifs",side=3, adj=-0.065,line=0.5,font=1,col='red')
  
  #Genes
  plotGenes(subbed, chrom, chromstart, chromend, types=subbed$type,bheight=0.05,
            plotgenetype="box",bentline=FALSE,col="black", labeloffset=.09,fontsize=0.8,arrowlength = 0.01,labeltext=TRUE,maxrows=10,wigglefactor=.51)
  labelgenome(chrom, chromstart,chromend,side=1,scipen=20,n=8,scale="Mb",line=.18,chromline=.5,scaleline=0.5)
  
  ##For saving pdf uncomment this
  dev.off()
  
}
fun_PlotMotifs(chrom, chromstart,chromend)
dev.off()
paste(chrom,":",chromstart,"-",chromend, sep="")



