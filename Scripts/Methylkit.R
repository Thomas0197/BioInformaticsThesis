rm(list=ls())

library(methylKit)
library(gplots)
library(ggplot2)
library(genomation)
library(GenomicRanges)

#Set the working directory 
setwd("~/Documents/erasmus")
list.files("~/Documents/erasmus")
#load(".RData")
load('General.RData')


#A quick look into one of the Sample files
C66E10 <- read.table("Extdata/Cases/C66E10_outputForMethylKit.tsv",sep="\t",header=T, nrows=20)

dim(C66E10)
head(C66E10)
str(C66E10)
summary(C66E10)

file.exists("Extdata/Cases/C66E10_outputForMethylKit.tsv")

#Reading in one of the Sample files, assembly is mm10, fraction is TRUE because of the format given in the output file
#pipeline specified due to the file format
myobj=methRead( "/home/thomas/Documents/erasmus/Extdata/Cases/C66E10_outputForMethylKit.tsv",
                pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                              coverage.col=5,strand.col=3,freqC.col=4),
                sample.id="C66E10",assembly="mm10")

# a look at the new object
dim(myobj)
head(myobj)
str(myobj)



#Stats on the file and filtering 

getMethylationStats(myobj,plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj,plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj,plot=TRUE,both.strands=FALSE)
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)






## Looking in multiple files/ comparisons

file.list=list(
  ("/home/thomas/Documents/erasmus/Extdata/Control/C66E9_outputForMethylKit.tsv"),
  ("/home/thomas/Documents/erasmus/Extdata/Control/C66E2_outputForMethylKit.tsv"),
  ("/home/thomas/Documents/erasmus/Extdata/Cases/C66E10_outputForMethylKit.tsv"),
  ("/home/thomas/Documents/erasmus/Extdata/Cases/C66E11_outputForMethylKit.tsv") 
)

file.list

#reading in the files
myobj=methRead(file.list,
               pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                             coverage.col=5,strand.col=3,freqC.col=4),
               sample.id=list("C66E9","C66E2","C66E10","C66E11"),
               assembly="mm10",treatment=c(1,1,0,0), mincov = 10)

# getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)

#Filtering to remove those with less then 10 coverage and the top 0.1 percentile for statistical reasoning 
myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)

myobj

#normalize data
myobj = normalizeCoverage(myobj)



# merge all samples to one object for base-pair locations that are covered in all samples
meth=unite(myobj, destrand=FALSE)
head(meth)

# exporting th new object
te = getData(meth)
te

write.table(te, file = "myDiff.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)



#Correlation table and charts for samples
getCorrelation(meth,plot=TRUE)

pdf =('graphs/Correlation.txt')
capture.output(getCorrelation(meth,plot=FALSE), file = pdf)


# Dendrogram for the cluster
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
jpeg('graphs/Dendrogram.jpg',width = 680, height = 780)
wardObject <-clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)
wardObject
plot(wardObject, main = "Cluster Dendrogram")
dev.off()




#PCA for the cluster

PCASamples(meth, screeplot=TRUE)

meth
PCA = PCASamples(meth,screeplot=FALSE,adj.lim=c(0.2,0.1), scale=TRUE,center=TRUE,comp=c(1,2),transpose=TRUE,sd.filter=TRUE,
                 sd.threshold=0.1,filterByQuantile=TRUE,obj.return=TRUE)

# Adding the treatment to the principle components to be able to colour them  
name <- c("Control1","Control2","Test1","Test2")
Pcomps = PCA$x
Pcomps <- data.frame(a,name)
Pcomps

PCData<-data.frame(PCA$x,Species=a$name)
PCData


jpeg('graphs/PCA.jpg', width = 740, height = 680)

autoplot(PCA,xlim = c(-800, 800),ylim = c(-800,800), scale = 0, label = TRUE, label.size = 3, shape = FALSE, 
         data = PCi, colour = 'Species') + scale_color_manual(values = c( "#37406b", "#37406b","#ab2727","#ab2727")) 
+theme_classic()

dev.off()


#Calculation of all the differential methylation present
myDiff=calculateDiffMeth(meth, slim = TRUE)
myDiff


# needed to be able to manipulate the data + check number of DMCs
Uncompressed_myDiff=getData(myDiff)
tail(Uncompressed_myDiff)
tail(Uncompressed_myDiff[[7]])
nrow(Uncompressed_myDiff)
#toVolcPlotQval <- data.frame(Uncompressed_myDiff[[6]])
#head(toVolcPlotQval)
#toVolcPlotQval1 <- toVolcPlotQval[toVolcPlotQval[,1]>=0.01,]
#head(toVolcPlotQval1)


#Methylation for hypo and hyper based of a Q value smaller or equal then 0.01
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hyper
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
myDiff25p.hypo
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
myDiff25p



write.table(methDiff25, file = "RegionMyDiff25.csv", sep = "\t",
            row.names = FALSE, col.names = TRUE)



##  Volcano Plot ##

toVolcPlot <- data.frame(Uncompressed_myDiff[[7]], -log10(Uncompressed_myDiff[[6]]))

toHighlightHyper1 <- toVolcPlot[toVolcPlot[,1]>25 & toVolcPlot[,2]>(-log10(0.01)),]
#nrow(toHighlightHyper1)

toHighlightHypo1 <- toVolcPlot[toVolcPlot[,1]< -25 & toVolcPlot[,2]>(-log10(0.01)),]
#nrow(toHighlightHypo1)

jpeg('graphs/volcanoPlot.jpg')

plot(toVolcPlot[,1], toVolcPlot[,2], pch=16,cex=0.3, col="grey", main="Volcano plot of Differentially Methylated CPGs", xlab = "Percentage %", ylab = "-log10 (pValue)")
points(toHighlightHyper1[,1], toHighlightHyper1[,2],pch=16,cex=0.32,col="#ab2727")
points(toHighlightHypo1[,1], toHighlightHypo1[,2],pch=16,cex=0.32,col="#37406b")
abline(a=-log10(0.01),b=0,col="black")
abline(v=25,col="black")
abline(v=-25,col="black")

dev.off()



# to export the methylation data 
toShorePlot <- data.frame(Uncompressed_myDiff[[1]], Uncompressed_myDiff[[2]],Uncompressed_myDiff[[3]], -log10(Uncompressed_myDiff[[5]]), Uncompressed_myDiff[[6]], Uncompressed_myDiff[[7]])
head(toShorePlot)
methDiff25 <- toShorePlot[(toShorePlot[,6]> 25 & toVolcPlot[,2]>(-log10(0.01))) | (toShorePlot[,6]< -25 & toVolcPlot[,2]>(-log10(0.01))),]
head(methDiff25)
nrow(methDiff25)

# variables for the hypo and hypermethylated dmcs
toHyper <- methDiff25[methDiff25[,6]>25 ,]
toHypo <- methDiff25[methDiff25[,6]< -25 ,]
head(toHypo)        
nrow(toHypo)



## Pie Chart ##

jpeg('graphs/pieChart.jpg')

slices <- c(nrow(toHyper),nrow(toHypo))
lbls <- c("Hyper", "Hypo")
pct <- round(slices/sum(slices)*100, 2)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
colors = c("#ab2727", "#37406b")
pie(slices,labels = lbls, col=colors,
    main="Hyper and Hypomethylated CPGs")

dev.off()



#########  For differential methylation per chromosome ############

diffMethPerChrUncomp <-diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)
diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)
diffMethPerChr
chr <- c("chr1","chr2","chr3" ,"chr4" ,"chr5" ,"chr6" ,"chr7" ,"chr8" , "chr9",  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX",  "chrY" )
chr


BarValues <- matrix(c(diffMethPerChrUncomp[[1]]$number.of.hypermethylated, diffMethPerChrUncomp[[1]]$number.of.hypomethylated), nrow = 2, ncol = 21, byrow = TRUE)


BarValues
BarValues[c(2),c(1)]/(BarValues[c(2),c(1)] + BarValues[c(1),c(1)])*100

HyperPerChr <- matrix(c(0), nrow = 1, ncol = 21, byrow = TRUE)
HypoPerChr <- matrix(c(0), nrow = 1, ncol = 21, byrow = TRUE)

#calculating the percentages of hyper/hypo for each chromosome and storing it
i<- 1
while(i<22){
  HyperPerChr[1,i]= HyperPerChr[1,i] + round(BarValues[c(1),c(i)]/(BarValues[c(1),c(i)] + BarValues[c(2),c(i)])*100)
  i = i+1
  #print(i)
}

HyperPerChr

i<- 1
while(i<22){
  HypoPerChr[1,i]= HypoPerChr[1,i] + round(BarValues[c(2),c(i)]/(BarValues[c(2),c(i)] + BarValues[c(1),c(i)])*100)
  i = i+1
}

HypoPerChr


#creating a dataframe to use when plotting
MethDiffPerChr <- matrix(c(HyperPerChr, HypoPerChr), nrow = 2, ncol = 21, byrow = TRUE)


#plotting the graph
jpeg('graphs/methPerChromoChart.jpg')

par(las=2)
barplot(MethDiffPerChr, main = "Hyper and Hypomethlyation per chromosome", names.arg = chr, xlab = "Percentage", ylab = "Chromosome", col = colors, horiz=TRUE)

dev.off()




#### Save to End ####

save.image("General.RData")







#### Regions ###


setwd('/home/thomas/ErasmusProject')
list.files()


file.list=list(
  ("/home/thomas/ErasmusProject/Extdata/Control/C66E9Alt.tsv"),
  ("/home/thomas/ErasmusProject/Extdata/Control/C66E2Alt.tsv"),
  ("/home/thomas/ErasmusProject/Extdata/Cases/C66E10Alt.tsv"),
  ("/home/thomas/ErasmusProject/Extdata/Cases/C66E11Alt.tsv") 
)

file.list

myobj=methRead(file.list,
               pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                             coverage.col=5,strand.col=3,freqC.col=4),
               sample.id=list("C66E9","C66E2","C66E10","C66E11"),
               assembly="mm10",treatment=c(1,1,0,0), mincov = 10)


myobj=filterByCoverage(tiles,lo.count=10,lo.perc=NULL,
                       hi.count=NULL,hi.perc=99.9)

# 'tiles' the genome with windows of 1000bp length and 1000bp step-size and summarizes the methylation information on those regions 
#(DMRs)
tiles = tileMethylCounts(myobj,win.size=1000,step.size=1000,cov.bases = 0)

myobj = normalizeCoverage(myobj)

meth=unite(tiles, destrand=FALSE)
head(meth)


myDiff=calculateDiffMeth(meth, slim = TRUE)
myDiff


Uncompressed_myDiff=getData(myDiff)
head(Uncompressed_myDiff)



##### Volcano###


toVolcPlot <- data.frame(Uncompressed_myDiff[[7]], -log10(Uncompressed_myDiff[[6]]))



toHighlightHyper1 <- toVolcPlot[toVolcPlot[,1]>25 & toVolcPlot[,2]>(-log10(0.01)),]
head(toHighlightHyper1)
nrow(toHighlightHyper1)

toHighlightHypo1 <- toVolcPlot[toVolcPlot[,1]< -25 & toVolcPlot[,2]>(-log10(0.01)),]
head(toHighlightHypo1)
nrow(toHighlightHypo1)

head(toVolcPlot)



jpeg('RegionVolcanoTestPlot.jpg')

plot(toVolcPlot[,1], toVolcPlot[,2], pch=16,cex=0.3, col="grey", main="Volcano plot of Differentially Methylated Regions", xlab = "Percentage %", ylab = "-log10 (qValue)")
points(toHighlightHyper1[,1], toHighlightHyper1[,2],pch=16,cex=0.32,col="#ab2727")
points(toHighlightHypo1[,1], toHighlightHypo1[,2],pch=16,cex=0.32,col="#37406b")
abline(a=-log10(0.01),b=0,col="black")
abline(v=25,col="black")
abline(v=-25,col="black")

dev.off()


## Pie Chart ###


toShorePlot <- data.frame(Uncompressed_myDiff[[1]], Uncompressed_myDiff[[2]],Uncompressed_myDiff[[3]], -log10(Uncompressed_myDiff[[5]]), Uncompressed_myDiff[[6]], Uncompressed_myDiff[[7]])
head(toShorePlot)
methDiff25 <- toShorePlot[(toShorePlot[,6]> 25 & toVolcPlot[,2]>(-log10(0.01))) | (toShorePlot[,6]< -25 & toVolcPlot[,2]>(-log10(0.01))),]
head(methDiff25)
nrow(methDiff25)


toHyper <- methDiff25[methDiff25[,6]>25 ,]
toHypo <- methDiff25[methDiff25[,6]< -25 ,]
head(toHypo)        
nrow(toHypo)


jpeg('RegionPie2Chart.jpg')


slices <- c(nrow(toHyper),nrow(toHypo))
slices
lbls <- c("Hyper", "Hypo")
pct <- round(slices/sum(slices)*100, 2)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
colors = c("#ab2727", "#37406b")
pie(slices,labels = lbls, col=colors,
    main="Hyper and Hypomethylated CPGs")

dev.off()


### Chromosome Plot ##

diffMethPerChrUncomp <-diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)
diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)
diffMethPerChr
chr <- c("chr1","chr2","chr3" ,"chr4" ,"chr5" ,"chr6" ,"chr7" ,"chr8" , "chr9",  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX",  "chrY" )
chr


BarValues <- matrix(c(diffMethPerChrUncomp[[1]]$number.of.hypermethylated, diffMethPerChrUncomp[[1]]$number.of.hypomethylated), nrow = 2, ncol = 21, byrow = TRUE)


BarValues

BarValues[c(2),c(1)]/(BarValues[c(2),c(1)] + BarValues[c(1),c(1)])*100

HyperPerChr <- matrix(c(0), nrow = 1, ncol = 21, byrow = TRUE)
HypoPerChr <- matrix(c(0), nrow = 1, ncol = 21, byrow = TRUE)
i<- 1

while(i<22){
  HyperPerChr[1,i]= HyperPerChr[1,i] + round(BarValues[c(1),c(i)]/(BarValues[c(1),c(i)] + BarValues[c(2),c(i)])*100)
  i = i+1
  #print(i)
}
HyperPerChr
i<- 1
while(i<22){
  HypoPerChr[1,i]= HypoPerChr[1,i] + round(BarValues[c(2),c(i)]/(BarValues[c(2),c(i)] + BarValues[c(1),c(i)])*100)
  i = i+1
}
HypoPerChr

MethDiffPerChr <- matrix(c(HyperPerChr, HypoPerChr), nrow = 2, ncol = 21, byrow = TRUE)


par(las=2)
barplot(MethDiffPerChr, main = "Hyper and Hypomethlyation per chromosome", names.arg = chr, xlab = "Percentage", ylab = "Chromosome", col = colors, horiz=TRUE)


