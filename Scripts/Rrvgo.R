if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rrvgo")


library(rrvgo)
library(org.Mm.eg.db)


#Setting the directory
setwd("~/Documents")



## Treemaps

#reading in the file
BP_File = read.delim('BioProsInput.txt', header = TRUE, sep = "\t", dec = ".")

BP_File

BP_File$ID

# creating the similarity matrix
simMatrix <- calculateSimMatrix(BP_File$ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")

go_analysis$ID
simMatrix

# getting the scores from the file
scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
scores

#Setting the scores to use
Scores <-setNames((BP_File$Score), BP_File$ID)
Scores
scores <- setNames(-log10(BP_File$P.value), BP_File$ID)


#selects as the group representative the term with the higher score within that group
reducedTerms <- reduceSimMatrix(simMatrix,
                                Scores,
                                threshold=0.7,
                                orgdb="org.Mm.eg.db")


scatterPlot(simMatrix, reducedTerms)
reducedTerms

treemapPlot(reducedTerms)







setwd("/media/thomas/Seagate Basic/erasmus/Extdata/RrvgoAnalysis/mWT-WT-Sup vs mWT-Df1-Sup-Normal/75")


BPInputFile = read.delim('BPInput.txt', header = TRUE, sep = "\t", dec = ".")

ID <- BPInputFile$ID
ID

CScore <- BPInputFile$Old.1
CScore


simMatrix <- calculateSimMatrix(ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")


simMatrix


Scores <-setNames((CScore), ID)

Scores


reducedTerms <- reduceSimMatrix(simMatrix,
                                Scores,
                                threshold=0.7,
                                orgdb="org.Mm.eg.db")

options(ggrepel.max.overlaps = Inf)

scatterPlot(simMatrix, reducedTerms)
scatterPlot(simMatrix,reducedTerms,addLabel = TRUE,labelSize = 2.95)
ggsave("ScatterPlot.pdf", plot = last_plot(), height = 10, width = 10)

reducedTerms

treemapPlot(reducedTerms)










