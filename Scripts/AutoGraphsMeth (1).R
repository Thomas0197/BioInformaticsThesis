#rm(list=ls())
library(methylKit)
library(gplots)
library(ggplot2)
library(GenomicRanges)

#Setting the directory
setwd('/home/thomas/ErasmusProject')

#Loading in the data
load('AutoDiffImage.RData')

# needed to be able to manipulate the data
Uncompressed_myDiff=getData(myDiff)


# For Pie chart 

#Methylation for hypo and hyper based of a Q value smaller or equal then 0.01
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")

myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")

myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

#Extracting diff meth above 25%
write.table(myDiff25p, file = "myDiff25_mWT-WT-Sup_VS_mDf1-WT-Sup.csv", sep = "\t",
            row.names = FALSE, col.names = TRUE)


## Volcano plot

toVolcPlot <- data.frame(Uncompressed_myDiff[[7]], -log10(Uncompressed_myDiff[[6]]))

VolcData <- data.frame(Uncompressed_myDiff[[7]], Uncompressed_myDiff[[6]], -log10(Uncompressed_myDiff[[5]]) )
testVolcPlot <- VolcData[VolcData[1] & VolcData[,2] <= 0.01,]

nrow(testVolcPlot)

#Seperating the hypo and hyper values
testHyper <- testVolcPlot[testVolcPlot[,1]>25 & testVolcPlot[,3]>(-log10(0.01)),]
nrow(testHyper)

testHypo <- testVolcPlot[testVolcPlot[,1]< -25 & testVolcPlot[,3]>(-log10(0.01)),]
nrow(testHypo)


jpeg('graphs/volcanoPlot_mWT-Df1-Sup-AAD_VS_mDf1-Df1-Sup-AAD.jpg',width = 580, height = 600)

plot(VolcData[,1], VolcData[,3], pch=16,cex=0.3, col="grey", main="Volcano plot of Differentially Methylated sites", xlab = "Percentage %", ylab = "-log10 (pValue)")
points(testHyper[,1], testHyper[,3],pch=16,cex=0.32,col="#ab2727")
points(testHypo[,1], testHypo[,3],pch=16,cex=0.32,col="#37406b")
abline(a=-log10(0.01),b=0,col="black")
abline(v=25,col="black")
abline(v=-25,col="black")

dev.off()


## Pie Chart

toShorePlot <- data.frame(Uncompressed_myDiff[[1]], Uncompressed_myDiff[[2]],Uncompressed_myDiff[[3]], -log10(Uncompressed_myDiff[[5]]), Uncompressed_myDiff[[6]], Uncompressed_myDiff[[7]])

methDiff25 <- toShorePlot[(toShorePlot[,6]> 25 & toVolcPlot[,2]>(-log10(0.01))) | (toShorePlot[,6]< -25 & toVolcPlot[,2]>(-log10(0.01))),]

toHyper <- methDiff25[methDiff25[,6]>25 ,]
toHypo <- methDiff25[methDiff25[,6]< -25 ,]

jpeg('graphs/pieChart_mWT-Df1-Sup-AAD_VS_mDf1-Df1-Sup-AAD.jpg')

slices <- c(nrow(toHyper),nrow(toHypo))
lbls <- c("Hyper", "Hypo")
pct <- round(slices/sum(slices)*100, 2)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
colors = c("#ab2727", "#37406b")
pie(slices,labels = lbls, col=colors,
    main="Hyper and Hypomethylated CPGs")

dev.off()


## Chart for meth per Chromo

diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)
diffMethPerChrUncomp <-diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

#Setting the labels
chr <- c("chr1","chr2","chr3" ,"chr4" ,"chr5" ,"chr6" ,"chr7" ,"chr8" , "chr9",  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX",  "chrY" )


BarValues <- matrix(c(diffMethPerChrUncomp[[1]]$number.of.hypermethylated, diffMethPerChrUncomp[[1]]$number.of.hypomethylated), nrow = 2, ncol = 24, byrow = TRUE)
BarValues <- subset(BarValues, select = -20)


#Setting them to percentages
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

jpeg('graphs/mWTDf1Normal_SupVsDef.jpg')

par(las=2)
barplot(MethDiffPerChr, main = "Hyper and Hypomethlyation per chromosome", names.arg = chr, xlab = "Percentage", ylab = "Chromosome", col = colors, horiz=TRUE)

dev.off()



##For The annotation files ##


toShorePlot <- data.frame(Uncompressed_myDiff[[1]], Uncompressed_myDiff[[2]],Uncompressed_myDiff[[3]], -log10(Uncompressed_myDiff[[5]]), Uncompressed_myDiff[[6]], Uncompressed_myDiff[[7]])
head(toShorePlot)
methDiff25 <- toShorePlot[(toShorePlot[,6]> 25 & toVolcPlot[,2]>(-log10(0.01))) | (toShorePlot[,6]< -25 & toVolcPlot[,2]>(-log10(0.01))),]
head(methDiff25)
nrow(methDiff25)

write.table(methDiff25, file = "myDiff25_mWT-WT-Sup_VS_mDf1-WT-Sup.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

#sed -i 's/\"//g' file.txt

#sed -i '1d' file.txt


print("Finished")

