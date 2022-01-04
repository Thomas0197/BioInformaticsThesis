rm(list=ls())
library(methylKit)
library(gplots)
library(ggplot2)

#Setting the directory
setwd('/home/thomas/ErasmusProject')

#Loading in the data to split
load('mWT-WT-Sup_VS_mWT-Df1-Sup-Normal.RData')


# Splitting for 40%
myDiff40p=getMethylDiff(myDiff,difference=40,qvalue=0.01)
myDiff40p

write.table(myDiff40p, file = "DiffMeth/myDiff40_mWTDf1DefNorm_VS_mWTDf1DefAAD.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE)


#Splitting for 75%
myDiff75p=getMethylDiff(myDiff,difference=75,qvalue=0.01)
myDiff75p

write.table(myDiff75p, file = "DiffMeth/myDiff75_mWTDf1DefNorm_VS_mWTDf1DefAAD.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)



write.table(myDiff25p, file = "DiffMeth/myDiff25_mWTDf1SupNorm_VS_mWTDf1SupAAD.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE)

