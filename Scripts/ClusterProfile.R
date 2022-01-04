if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")


library(ggplot2)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(ggnewscale)
library(ggnewscale)



#Setting the directory 
setwd("~/Downloads")


x <- c("1700012A03Rik","1700018B08Rik","1700040L02Rik","1700123L14Rik","4930451G09Rik","4930595D18Rik","4933428G20Rik","Acp6","Actr10","Adam30","Adamts9","Adcy5","Agap1","Ak9","Alk","Ambn","Amdhd1","Ankfn1","Arhgap35","Asah2","Axin2","Azi2","BC048507","BC048679","Bambi","Cap2","Cd164","Cdh13","Cetn3","Cfap73","Chrna4","Chst12","Clip4","Cog6","Cpb2","Cpn2","Cpt1b","Csnk1g3","Csrp1","Cst10","Cyp11a1","Dap","Ddit4l","Ddx3x","Ddx43","Depdc1a","Dhodh","Dlst","Dlx2","Dym","Ece1","Edn2","Efna5","Eif2b2","Eif4e2","Eif4ebp1","Eif4g2","Eif4g3","Ets1","Etv6","Faiml","Fam110c","Fam19a4","Fam8a1","Fbxo11","Fgf6","Fkbp4","Fsd2","Fxn","Gabbr2","Gdi2","Gnai2","Gnas","H2-M9","H2afz","Hal","Hand1","Hao2","Heg1","Hhla1","Hif1a","Hmgb1","Hnrnpa0","Hnrnpul1","Hspb8","Htr7","Irx3","Kansl1","Kcnf1","Kcng1","Kcnq1","Kdr","Krt9","L3mbtl4","Lama2","Larp1","Lmo7","Lpp","Lrrc47","Lrrc7","Ly6g6e","Lypd4","Mapkap1","Mb21d1","Mbd2","Meaf6","Mpp6","Mrps27","Msn","Mtmr12","Mturn","Mvd","Myom1","Nabp1","Nav1","Ncf4","Ndufa8","Ndufb2","Ndufv3","Nebl","Nefl","Nsd1","Nt5c1b","Nudc","Obox1","Olfm3","Pak7","Pax3","Pbx1","Pcif1","Plagl1","Pnpla1","Pphln1","Ppif","Ppp1r15b","Ppp6r3","Prkag2","Prpf40a","Psapl1","Psg16","Psma3","Psmd7","Pwp1","Pygl","Rbfox3","Rbms1","Rell2","Rnf111","Rnf216","Rprd1b","Sash1","Sccpdh","Sec23ip","Sel1l2","Sept9","Sez6l","Shroom3","Siah3","Ski","Skint6","Slc25a3","Slc25a33","Slc27a5","Slc35c2","Slc38a4","Slc4a10","Slc6a6","Sox4","Sp3","Spag9","Spsb1","Srrm4","Ss18","Stag2","Stt3b","Sulf2","Sult1c2","Sytl1","Tagap","Tcf4","Tdpoz3","Tex29","Tgfbr3","Tgm2","Thegl","Tmem17","Tmem200c","Tmem45a2","Tmem68","Tnks2","Tns1","Tsc22d2","Ttll11","Ttn","Usp46","Usp9x","Vapa","Vwa8","Wdr41","Wdr95","Wnt9b","Xpot","Ywhag","Ywhaz","Zdhhc25","Zfp507","Zfp957","Zmat4","Zmiz1","Zmynd8","Zscan2")
eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg)
tail(eg)
Genes <- eg$ENTREZID

Genes

#converting the gene list
gene.df <- bitr(Genes, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Mm.eg.db)

head(gene.df)

#GO over-representation analysis
ggo <- groupGO(gene     = Genes,
               OrgDb    = org.Mm.eg.db,
               ont      = "BP",
               level    = 12,
               readable = TRUE)

head(ggo)
ggo
newdata <- ggo[order(-ggo$Count),] 
head(newdata)
newdata$geneID <- NULL
newdata
head(newdata,25)



ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
head(ego2)
ego2
goplot(ego2)
gene.df$ENSEMBL

y <-sort(gene.df$ENTREZID, decreasing = TRUE)
y

#gene set enrichment analysis using gene ontology.
ego3 <- gseGO(geneList     = y,
              OrgDb        = org.Mm.eg.db,
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

goplot(ego3)


edo <- enrichDGN(Genes)

edo

edo$g

## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(Genes, foldChange=)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))





