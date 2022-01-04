rm(list=ls())
library(methylKit)
library(gplots)
library(ggplot2)
library(annotatr)
library(AnnotationHub)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

#Setting the directory 
setwd('/home/thomas/ErasmusProject')
list.files()

#Annotations to find
annots = c('mm10_cpgs', 'mm10_basicgenes', 'mm10_genes_intergenic', 'mm10_genes_firstexons', 'mm10_genes_intronexonboundaries')

#retrieving the annotations from the genome
annotations = build_annotations(genome = 'mm10', annotations = annots)

#reading in the DMCs
dm_file = ('myDiff25_mWT-WT-Sup_VS_mDf1-WT-Sup.txt')
TestCols = c(p_value = 'numeric', q_value = 'numeric', meth_diff = 'numeric')
dm_regions = read_regions(con = dm_file, genome = 'mm10', extraCols = TestCols, format = 'bed')

#Annotating the DMCs
dm_annotated = annotate_regions(
  regions = dm_regions,
  minoverlap = 0L,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

print(dm_annotated)


#Summarizing the results
dm_annsum = summarize_annotations(
  annotated_regions = dm_annotated,
  quiet = TRUE)
df_dm_annsum = data.frame(dm_annsum)



dm_catsum = summarize_categorical(
  annotated_regions = dm_annotated,
  by = c('annot.type'),
  quiet = TRUE)
df_dm_catsum = data.frame(dm_catsum)

df_dm_annsum



#Genomic Plot
genomicPlot_df <-df_dm_annsum[c(8,12,13,10),] 
genomicPlot_df 

sum_genomicPlot = sum(genomicPlot_df$n)
percentage_genomicPlot <- matrix(c(0), nrow = 4, ncol = 2, byrow = TRUE)
percentage_genomicPlot

i<- 1
while(i<5){
  percentage_genomicPlot[i,2]= percentage_genomicPlot[i,2] + round((genomicPlot_df[i,2]/sum_genomicPlot)*100,2)
  i = i+1
}
percentage_genomicPlot


annot.type = genomicPlot_df[,1]

n = percentage_genomicPlot[,2]

percentage_genomicPlot <- data.frame(annot.type, n)
percentage_genomicPlot


jpeg('graphs/GenomicContext_mWT-WT-Sup_VS_mDf1-WT-Sup.jpg', width = 680, height = 780)
par(mar=c(6,4,4,4))

colors = c("#cad6fe", "#90b1c9","#7590a0","#325e7a")
Genomic_plot <- barplot(height=percentage_genomicPlot$n , names.arg=c("Exons","Introns","Promoters","Intergenic"), ylim=c(0,50), col = colors, ylab="Percentage (%)", main="Annotation (Genomic context) of Differentially Methylated CpGs")
text(Genomic_plot, percentage_genomicPlot$n + 2, paste(percentage_genomicPlot$n, "%"), cex=1)
dev.off()



#CPG islands Plot

islandPlot_df <- df_dm_annsum[c(2,4,3,1),]
islandPlot_df
sum_islandPlot = sum(islandPlot_df$n)
percentage_islandPlot <- matrix(c(0), nrow = 4, ncol = 2, byrow = TRUE)
percentage_islandPlot



i<- 1
while(i<5){
  percentage_islandPlot[i,2]= percentage_islandPlot[i,2] + round((islandPlot_df[i,2]/sum_islandPlot)*100,2)
  i = i+1
}
percentage_islandPlot

annot.type = islandPlot_df[,1]

n = percentage_islandPlot[,2]

percentage_islandPlot <- data.frame(annot.type, n)
percentage_islandPlot


jpeg('graphs/CpgContext_mWT-WT-Sup_VS_mDf1-WT-Sup.jpg', width = 680, height = 780)
colors = c("#91be5f", "#dfca76","#00937d","#00637d")
Genomic_plot <- barplot(height=percentage_islandPlot$n , names.arg=c("Islands","Shores","Shelves","Open Sea") , ylim=c(0,100), col = colors, ylab="Percentage (%)", main="Annotation (CpG context) of Differentially Methylated CpGs")
text(Genomic_plot, percentage_islandPlot$n + 2, paste(percentage_islandPlot$n, "%"), cex=1)
dev.off()

save.image("mWT-WT-Sup_VS_mDf1-WT-Sup.RData")

