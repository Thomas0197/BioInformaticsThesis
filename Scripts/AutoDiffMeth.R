rm(list=ls())
library(methylKit)
library(gplots)
library(ggplot2)
library(GenomicRanges)

#setting the directory
setwd('/home/thomas/ErasmusProject')

#Reading in the files 
file.list=list(
  ("/home/thomas/ErasmusProject/Extdata/All/Alt/C66E10_Alt.tsv"),
  ("/home/thomas/ErasmusProject/Extdata/All/Alt/C66E11_Alt.tsv"),
  ("/home/thomas/ErasmusProject/Extdata/All/Alt/C64E6_Alt.tsv"),
  ("/home/thomas/ErasmusProject/Extdata/All/Alt/C64E7_Alt.tsv") 
)

myobj=methRead(file.list,
               pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                             coverage.col=5,strand.col=3,freqC.col=4),
               sample.id=list("C66E10","C66E11","C64E6","C64E7"),
               assembly="mm10",treatment=c(1,1,0,0), mincov = 10)

#Filtering the samples
myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                       hi.count=NULL,hi.perc=99.9)


myobj = normalizeCoverage(myobj)


# merge all samples to one object for base-pair locations that are covered in all samples
meth=unite(myobj, destrand=FALSE)


getCorrelation(meth,plot=FALSE)

#Correlation
pdf =('graphs/Correlation_mWT-WT-Sup_VS_mDf1-WT-Sup.txt')
capture.output(getCorrelation(meth,plot=FALSE), file = pdf)

#Dendrogram
jpeg('graphs/Dendrogram_mWT-WT-Sup_VS_mDf1-WT-Sup.jpg', width = 680, height = 780)
wardObject <-clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)
plot(wardObject, main = "Cluster Dendrogram")
dev.off()

jpeg('graphs/PCA_mWT-WT-Sup_VS_mDf1-WT-Sup.jpg', width = 680, height = 680)
PCASamples(meth,scale=TRUE,center=TRUE, adj.lim = c(0.6,01.2))
dev.off()


#Calculation of all the differential methylation present
myDiff=calculateDiffMeth(meth, slim = TRUE)

#Saving the comparison
save.image("NoneAutoDiffImage.RData")
