##Title: Combining transcript read counts from multiple files into a counts matrix of gene read counts
#James George 04/21 #2022-2024 Modified by E. Young

library(tidyverse)
library(limma)
library(edgeR)
library(fgsea)
library(ggplot2)
library(EnhancedVolcano)
library(gtable)
library(gplots)
library(openxlsx)

viewer <- function(x){
  return( view(as.data.frame(x)) )
}

############### First half ###############
# Change working directory to directory with transcript read counts
wd = "//Datasets/RNAseq/SanRepeats_19762_human/"
setwd(wd)
filerange = paste0('SI_', 19762:19773)

# Specify organism for transcript to gene key
transcript_to_gene_key = read.table("//Datasets/RNAseq/hg38_transcripts_to_genes_test.txt",header = T) #swap for mouse

####Create counts matrix with transcript reads, merge with gene names, sum up all read counts for transcripts with same gene names, and write a new read counts file into the same directory ####
libs=paste0(filerange)
cts=lapply(libs, function(x) {
  y = read.table(paste0(x, '/', 'abundance.tsv'), stringsAsFactors = F, header = TRUE)
  z = right_join(transcript_to_gene_key, y, by = 'target_id')
  a = aggregate(est_counts ~ gene_id + gene_name, data = z, sum)
  write.table(a, paste0(x, '/', 'gene-abundance.tsv'))
  
  b = aggregate(tpm ~ gene_id + gene_name, data = z, sum)
  write.table(b, paste0(x, '/', 'TPM.tsv'))
})

#### Create final counts matrix ####
dir.create('read_counts')
readcounts.title = sub(".*/", "", wd)
readcounts_fullprefix = paste0(wd,'/read_counts/', readcounts.title)
libs=paste0(filerange)
cts=lapply(libs, function(x) {
  read.table(paste0(x, '/','gene-abundance.tsv'), stringsAsFactors = F, header = TRUE)
})

#make count matrix
cts.ma=do.call(cbind, lapply(cts, function(x) x$est_counts))
colnames(cts.ma) = libs
rownames(cts.ma) = cts[[1]]$gene_id
write.table(cts.ma, paste0(readcounts_fullprefix,'.cts.txt'))
write.table(cts[[1]][,1:2], paste0(readcounts_fullprefix,'.cts.genes.txt'))

# Also compile TPMs
#find libraries & load file containing gene TPMs
libs=paste0(filerange)
tpm=lapply(libs, function(x) {
  read.table(paste0(x, '/','TPM.tsv'), stringsAsFactors = F, header = TRUE)
})

#make tpm matrix
tpm.ma=do.call(cbind, lapply(tpm, function(x) x$tpm))
colnames(tpm.ma) = libs
rownames(tpm.ma) = tpm[[1]]$gene_id
write.table(tpm.ma, paste0(readcounts_fullprefix,'.TPMs.txt'))

############### Second half ###############

#### Section 1: Identify DE genes with limma-voom ####

# Change working directory
setwd("//PATHTOYOURDATA/read_counts")


# read raw count
cts2 =read.table("TMP1.csv", header = T)
cts.name = read.table("SanRepeats_19762_human.cts.genes.txt", header = T)
colnames(cts.name) = c('gene_id','gene_name')

# construct DGEList
dge=DGEList(counts = cts)

# normalize for library size
dge=calcNormFactors(dge)

#EY EDIT HERE EACH RUN
gp = c('n22RV1_Parental','n22RV1_Parental','n22RV1_57','n22RV1_57','n22RV1_WT2','n22RV1_WT2','n22RV1_47','n22RV1_47','n22RV1_73','n22RV1_73','n22RV1_WT3','n22RV1_WT3')
gp = factor(gp, levels = c('n22RV1_Parental','n22RV1_57','n22RV1_WT2','n22RV1_47','n22RV1_73','n22RV1_WT3'))

des=model.matrix(~0+gp)
colnames(des)=stringr::str_replace(colnames(des), 'gp','')
des 

#@ey change column names again
# Create average cpm matrix 
cpm = as.data.frame(cpm(dge, normalized.lib.sizes = TRUE))

#record ID column names (Paste back in if sending to collaborators)
colnames(cpm) 
write.csv(colnames(cpm) , 'CPM_Names.csv', row.names= F)
#EY EDIT HERE EACH RUN
colnames(cpm) = paste0(rep(c('n22RV1_Parental','n22RV1_57','n22RV1_WT2','n22RV1_47','n22RV1_73','n22RV1_WT3'), each = 2), rep(c('.R1','.R2'),6))
colnames(cpm)

avg_cpm <- as.data.frame(matrix(nrow= nrow(cpm)))
row.names(avg_cpm) = row.names(cpm)
counter = 1

for(x in seq(1,ncol(cpm), by=2)){
  
  # create a subset from the dataframe and compute the mean of the rows and finally cbind it to the result dataframe
  avg_cpm[ncol(avg_cpm)] <- apply(subset(cpm, select=seq(x,length.out =2)), 1, mean, na.rm = TRUE)
  
  #EY EDIT HERE 
  colnames(avg_cpm)[ncol(avg_cpm)] <- c('n22RV1_Parental','n22RV1_57','n22RV1_WT2','n22RV1_47','n22RV1_73','n22RV1_WT3')[counter]
  avg_cpm$hold <- NA
  counter = counter + 1
}

avg_cpm = avg_cpm[,-ncol(avg_cpm)]
avg_cpm = rownames_to_column(avg_cpm, var = 'gene_id')
avg_cpm = left_join(avg_cpm, cts.name, by = 'gene_id')

cpm = rownames_to_column(cpm, var = 'gene_id')
cpm = left_join(cpm, cts.name, by = 'gene_id')

write.csv(cpm, 'CPM_Values.csv', row.names= F)

# Filter 1:  Remove all genes with active CPMs less than 1 (for down-regulated genes, CPM has to be >1 in NC, for up-regulated genes, CPM has to be >1 in treated)
#cpm_flt = (with(avg_cpm, pmax(sgNC, sgTrp53)) >1) | (with(avg_cpm, pmax(sgNC, CRE_sgTrp53)) >1)
# Filter 2: keep genes with average counts greater than 10
flt=rowMeans(cts)> 10
# voom
v=voom(dge[flt, ] , design = des, plot=T )
# limma linear model
vfit=lmFit(v, des)

Color <- c("steelblue2","steelblue2","palevioletred3","palevioletred3","darkolivegreen4","darkolivegreen4","darkslateblue","darkslateblue","snow4","snow4","darkgoldenrod1","darkgoldenrod1")
pdf("PCA.pdf")
plotMDS(v$E,gene.selection="common",col=Color)
colnames(v$E)=paste0(rep(c('n22RV1_Parental','n22RV1_57','n22RV1_WT2','n22RV1_47','n22RV1_73','n22RV1_WT3'), each = 2), rep(c('.R1','.R2'),6))
dev.off()

#@EY EDIT HERE 
cnts=makeContrasts(n22RV1_57 - n22RV1_Parental,
                   n22RV1_57 - n22RV1_WT3,
                    levels = des)

cfit=contrasts.fit(vfit, cnts)
cfit=eBayes(cfit)

# extract differential analysis results
tabs=lapply(1:ncol(cfit), function(i) topTable(cfit, coef = i, number = Inf, adjust.method = "BH", sort.by = 'none'))
names(tabs) <- colnames(cfit)

# add gene names to tabs
tabs = lapply(tabs, function(x){
  x$gene_id = row.names(x)
  x = inner_join(cts.name,x, by ='gene_id')
})


# Attach cpm values
hold = strsplit(names(tabs),' - ' )
for (i in 1:length(tabs)){
  tabs[[i]][,paste0('CPM_', hold[[i]][1]) ] = avg_cpm[avg_cpm$gene_id %in% tabs[[i]]$gene_id, hold[[i]][1]]
  tabs[[i]][,paste0('CPM_', hold[[i]][2]) ] = avg_cpm[avg_cpm$gene_id %in% tabs[[i]]$gene_id, hold[[i]][2]]
}


#### Volcano plots and significant genes ####
# Establish cutoffs for significant DEGs
logFC_cutoff = 1
pval_cutoff = 0.05

#Designate significant UP and DOWN genes for all treatments
sig_genes <- lapply(tabs, function(x){
  UP= x[x$logFC> logFC_cutoff & x$adj.P.Val<pval_cutoff,]
  DOWN = x[x$logFC< -logFC_cutoff & x$adj.P.Val<pval_cutoff,]
  return(list(UP, DOWN))
})

sig_genes_final = list()
for (i in 1:length(sig_genes)){
  sig_genes_final = c(sig_genes_final, sig_genes[[i]][1], sig_genes[[i]][2])
}
names(sig_genes_final)= paste(rep(names(tabs), each=2), rep(c('UP','DN'), (length(names(tabs))) ), sep = '_')
names(sig_genes_final)=stringr::str_replace_all(names(sig_genes_final), '-', '.v.')
names(sig_genes_final)=stringr::str_replace_all(names(sig_genes_final), ' ', '')

#enhanced volcano plot
tabs_sub = tabs
s <- split(c(1:4), 1:2)
e = lapply(1:length(tabs_sub), function(i) EnhancedVolcano(toptable = tabs_sub[[i]], lab = tabs_sub[[i]]$gene_name, x = 'logFC', y = 'adj.P.Val',
                                                           xlim = c(min(tabs_sub[[i]]$logFC, na.rm = TRUE) - 0.5, max(tabs_sub[[i]]$logFC, na.rm = TRUE) + 0.5),
                                                           #xlim = c(-8, 10),
                                                           ylim = c(0, max(-log10(tabs_sub[[i]]$adj.P.Val), na.rm = TRUE) + 0.5),
                                                           xlab = bquote(~Log[2] ~ "fold change"), ylab = bquote(~-Log[10] ~ italic(FDR)), axisLabSize = 15,
                                                           title = paste(names(tabs_sub)[[i]], 'Genes'), subtitle = paste0('  Down Genes:', nrow(sig_genes_final[[s$`2`[i]]]), '          Up Genes:', nrow(sig_genes_final[[s$`1`[i]]]) ),
                                                           caption= bquote(~Log[2]~'FC cutoff:'~.(logFC_cutoff)~'; P-value cutoff:'~.(pval_cutoff)~'. P-value adjusted by FDR'),
                                                           subtitleLabSize = 10, captionLabSize = 10,
                                                           pCutoff = pval_cutoff, FCcutoff = logFC_cutoff, cutoffLineType = 'dashed',
                                                           col = c("grey30", "forestgreen", "royalblue", "red2"),
                                                           legendLabels = c("NS", expression(log[2] ~ FC), "p-value", expression(p-value ~ and ~ log[2] ~ FC)),
                                                           legendPosition = 'top'))
e
pdf(paste0('Volcano Plots for   (logFC_', logFC_cutoff,' & P-Val_', pval_cutoff, ').pdf'))
e
dev.off()



# Create a blank workbook # Excel sheet of significant DEGs
OUT <- createWorkbook()
#max 31 char EDIT HERE EACH TIME
short_hand_name = paste(rep(c('n22RV1_57n22RV1_Parental',
                   'n22RV1_57-n22RV1_WT3'), each=2), rep(c('UP','DN'), (length(names(tabs))) ), sep = '_')                     
for(x in 1:length(sig_genes_final)){
  y = as.data.frame(sig_genes_final[[x]])
  y= y[c('gene_id','gene_name','logFC','adj.P.Val', paste(names(sig_genes_final[[x]])[grep('CPM', names(sig_genes_final[[x]]))]) ) ]
  addWorksheet(OUT, sheetName = short_hand_name[x] )
  #writeData(OUT, sheet = short_hand_name[x], y)
}

# Export the file
saveWorkbook(OUT, paste0('Significant_DEGs(logFC_', logFC_cutoff,' & P-Val_', pval_cutoff, ') for_set_RNAseq.xlsx'))

hallmark = gmtPathways('/PATHTOGMTS/hallmark.gmt')

rnks=lapply(tabs, function(x) {
  x=x[order(x$logFC, decreasing =T),]
  #  x = subset(x, adj.P.Val < 0.1)
  x$gene_name = sapply(x$gene_name, function(x) strsplit(x, '[.]')[[1]][1])
  print(x[duplicated(x$gene_name), 'gene_name'])
  x=x[!duplicated(x$gene_name),]  #Cut out duplicates
  r=x$logFC
  names(r)=x$gene_name
  return(r)
})

gsea.res=lapply(rnks,function(x) fgsea(hallmark,x,10, 2000))

order.gsea.res = lapply(gsea.res, function(x){
  x= x[order(x$NES, decreasing = FALSE),]
  x$pathway = factor(x$pathway, levels = x$pathway)
  return(x)
})

# Summarize pathway NES values
p_value_NES_barplot = lapply(1: length(gsea.res), function(i) ggplot(order.gsea.res[[i]], aes(x=pathway, y=NES, label=NES)) + 
                               geom_bar(stat='identity', aes(fill=NES, alpha = factor(ifelse(padj<0.05, 1, 0.3))), width=.5)  +
                               scale_alpha_manual(values = c("1"=1, "0.3"=0.3), guide='none') +
                               scale_fill_gradient2(name="NES", 
                                                    low = "royalblue4",
                                                    mid = "white",
                                                    high = "red4") +
                               labs(title= paste(names(gsea.res)[[i]], "NES Values for Hallmark Gene Sets"))+
                               coord_flip() +
                               theme(axis.text.y = element_text(size = 6.5),
                                     plot.title = element_text(hjust = 0.5, face= "bold")) )


p_value_NES_barplot


for (i in 1: length(rnks)){
  
  current_pathway = hallmark
  
  assign(paste0('g', i), lapply(1: length(current_pathway), function(x) plotEnrichment(current_pathway[[x]], rnks[[i]]) + 
                                  labs(title= paste(names(current_pathway)[x], 'Pathway'), 
                                       x= paste(names(rnks)[[i]], 'Ranked List'),
                                       caption = c(paste('NES: ', signif(gsea.res[[i]][pathway== names(current_pathway)[x], NES], digits = 4), 
                                                         ' & p-value: ', signif( gsea.res[[i]][pathway== names(current_pathway)[x], padj], digits = 4), '   ',
                                                         '  Min:',    signif( summary(rnks[[i]])[1], digits = 4) ,
                                                         '  1st Qu:', signif( summary(rnks[[i]])[2], digits = 4) ,
                                                         '  Median:', signif( summary(rnks[[i]])[3], digits = 4) ,
                                                         '  Mean:',   signif( summary(rnks[[i]])[4], digits = 4) ,
                                                         '  3rd Qu:', signif( summary(rnks[[i]])[5], digits = 4) ,
                                                         '  Max:',    signif( summary(rnks[[i]])[6], digits = 4) )) ) +
                                  theme(plot.caption = element_text(size= 7,lineheight = 1, hjust=0),
                                        axis.title.x = element_text(size = 12),
                                        plot.title = element_text(face= 'bold', size = 16, hjust= 0.5)) ))
  
  pdf(file = paste(' Name ', names(rnks)[i], '.pdf') )
  print(p_value_NES_barplot[[i]])
  print(get(paste0('g', i)))
  dev.off()
}

##save.image("DATA.Rdata")
