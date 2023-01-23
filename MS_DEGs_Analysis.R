setwd("~/OneDrive - University Of Houston/WU_LAB/SLE_MS/Publish_Data/")

load("~/OneDrive - University Of Houston/WU_LAB/SLE_MS/Publish_Data/MS_DEGs_Analysis.RData")
save.image("~/OneDrive - University Of Houston/WU_LAB/SLE_MS/Publish_Data/MS_DEGs_Analysis.RData")


library("ggplot2")
library("readxl")
library("dplyr")
library("reshape2")
library("stringr")
library("xlsx")
library("pheatmap")
library("heatmaply")
library("ggpubr")
library("ggrepel")
library("VennDiagram")
library(venn)
library(VennDiagram)

## AG data ------------------------------------------------
AG <- read_excel("MS_RawData.xlsx", sheet = "AG")
AG <- as.data.frame(AG)
AG[is.na(AG)] <- 0

head(AG)
dim(AG)
AG$logFC <- log2(AG$`SLE Score`/AG$`HC Score`)
str(AG)

AG_DEGs <- AG[AG$logFC> 0.5, ]
AG_DEGs <- AG_DEGs[-which(is.na(AG_DEGs$Accession)), ]
dim(AG_DEGs)
## 87

head(AG_DEGs)

write.table(AG_DEGs, file = "AG_DEGs_87.csv", quote = FALSE, row.names = FALSE, sep = '\t')

AG_LN <- AG[AG$logFC == "Inf", ]
dim(AG_LN) ## 53
write.table(AG_LN, file = "AG_LN-ONLY.csv", sep = "," ,quote =FALSE, row.names = FALSE)

## C1Q data --------------------------------------------
C1Q <- read_excel("MS_RawData.xlsx", sheet = "C1Q")
C1Q <- as.data.frame(C1Q)
C1Q[is.na(C1Q)] <- 0

C1Q$logFC <- log2(C1Q$`SLE Score`/C1Q$`HC Score`)
str(C1Q)

C1Q_DEGs <- C1Q[C1Q$logFC> 0.5, ]
C1Q_DEGs <- C1Q_DEGs[-which(is.na(C1Q_DEGs$Accession)), ]
dim(C1Q_DEGs)
## 30
write.table(C1Q_DEGs, file = "C1Q_DEGs_30.csv", sep = "\t" ,quote =FALSE, row.names = FALSE)

C1Q_LN <- C1Q[C1Q$logFC == "Inf", ]
dim(C1Q_LN) ## 27
write.table(C1Q_LN, file = "C1Q_LN-ONLY.csv", sep = "," ,quote =FALSE, row.names = FALSE)

## Combine AG & C1Q -----------------------------------------------------

head(AG_DEGs)
head(C1Q_DEGs)

inner <- inner_join(AG_DEGs, C1Q_DEGs, by = c("Accession" = "Accession"))

venn.diagram(list(AG = AG_DEGs$Accession, C1Q = C1Q_DEGs$Accession),
             resolution = 300, imagetype = "png", alpha=c(0.5,0.5),
             fill=c("red","blue"), scaled = TRUE, main= "C1Q-AG_DEGs_Compare",
             main.cex = 2,
             print.mode = "raw",
             filename = "AG_DEGs_Compare.png")
)

inner_ONLY <- inner_join(AG_LN, C1Q_LN, by = c("Accession" = "Accession"))



### GO -----------------------------------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi) ## install.packages(c("Rcpp")) if select not working
keytypes(org.Hs.eg.db)

## read DAVID convert ID list
AG_cv <- read_excel("AG_LN-ONLY.xlsx", sheet = "ID_Convert")
C1Q_cv <- read_excel("C1Q_LN-ONLY.xlsx", sheet = "ID_Convert")

GO_List <- as.character(unique(c(AG_cv$To, C1Q_cv$To)))

DEGs_ENID <- select(org.Hs.eg.db, keys=GO_List, keytype = "ENTREZID", columns= c("GENENAME","SYMBOL"))
dim(DEGs_ENID)
head(DEGs_ENID)

## MF
MF <- enrichGO(DEGs_ENID$ENTREZID, "org.Hs.eg.db", ont = "MF", keyType = "ENTREZID",
               pAdjustMethod = "BH", pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.1, readable = TRUE)
dim(MF) ## 8

MF_Sim <- simplify(MF, cutoff = 1, by = "p.adjust", select_fun = min)
dim(MF_Sim) ## 3

dotplot(MF_Sim, showCategory = 8)

## CC
CC <- enrichGO(DEGs_ENID$ENTREZID, "org.Hs.eg.db", ont = "CC", keyType = "ENTREZID",
               pAdjustMethod = "BH", pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.1, readable = TRUE)
dim(CC) ##11

CC_Sim <- simplify(CC, cutoff = 0.8, by = "p.adjust", select_fun = min)
dim(CC_Sim) ## 8

dotplot(CC_Sim, showCategory = 8)

## BP
BP <- enrichGO(DEGs_ENID$ENTREZID, "org.Hs.eg.db", ont = "BP", keyType = "ENTREZID",
               pAdjustMethod = "BH", pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.1, readable = TRUE)
dim(BP) ##55

BP_Sim <- simplify(BP, cutoff = 0.8, by = "p.adjust", select_fun = min)
dim(BP_Sim) ## 20

dotplot(BP_Sim, showCategory = 20)




## Database Cross- validation ---------------------------------------------------------

## read scRNA LN-Health DEGs ----------------------------------------------------------------------------
## 6240; (P < 0.05 & logFC >1)

scRNA <- read_excel("scRNA_LN-HC_DEGs_05-fc.xlsx", sheet = "scRNA_LN-HC_DEGs_05-fc")
colnames(scRNA)
dim(scRNA) ##6320

## read One LN Microarray GSE32591 -------------------------------------------------------
## 151; (P < 0.05 & logFC >0.5)
LN_MA_One <- read.csv("../R Analysis/Query_MS_Result/GSE32591/GSE32591-Microarray_DEGs_531.csv", header = TRUE, row.names = 1,sep = ",")
dim(LN_MA_One) ## 531
LN_MA_One$Gene.Symbol


## GSE81622 LN Illumina WGT ------------------------------------------------------
## 627; (P < 0.05 & logFC >0.5)
LN_MA_622 <- read.csv("../R Analysis/Query_MS_Result/GSE81622/GSE81622_LN_DEGs.csv", header = TRUE, row.names = 1,sep = ",")
dim(LN_MA_622) ## 627
LN_MA_622$Gene.symbol


## GSE99967 LN Affymetrix ------------------------------------------------------
## 745; (P < 0.05 & logFC >0.5)
LN_MA_967 <- read.csv("../R Analysis/Query_MS_Result/GSE99967/GSE99967-Microarray_DEGs_745.csv", header = TRUE, row.names = 1,sep = ",")
dim(LN_MA_967) ## 745
LN_MA_967$GeneSymbol


## Venn  --------------------------------------------------------------------------------
venn.diagram(
  list(LN_MS = DEGs_ENID$SYMBOL, 
       scRNA = scRNA$gene, 
       LN_MA_967 = LN_MA_967$GeneSymbol,
       LN_MA_591 = LN_MA_One$Gene.Symbol,
       LN_MA_622 = LN_MA_622$Gene.symbol
  ),
  resolution = 300, imagetype = "png", alpha = c(0.5, 0.5, 0.5, 0.5, 0.5),
  fill= c("#1B9E77", "#7570B3", "#E7298A", "#D95F02", "#E6AB02"), 
  scaled = TRUE, main= "Database Cross-Validation",
  main.cex = 2,
  print.mode = "raw",
  filename = "Database_Cross-Validation.png"
)


head(scRNA)

scRNA__LN_MS <- inner_join(scRNA, DEGs_ENID, by = c("gene" = "SYMBOL"))
dim(scRNA__LN_MS) ## 14
head(scRNA__LN_MS)

LN_MA_591__LN_MS <- inner_join(LN_MA_One, DEGs_ENID, by = c("Gene.Symbol" = "SYMBOL"))
dim(LN_MA_591__LN_MS) ## 3
head(LN_MA_591__LN_MS)

LN_MA_622__LN_MS <- inner_join(LN_MA_622, DEGs_ENID, by = c("Gene.symbol" = "SYMBOL"))
dim(LN_MA_622__LN_MS) ## 1
head(LN_MA_622__LN_MS)

LN_MA_967__LN_MS <- inner_join(LN_MA_967, DEGs_ENID, by = c("GeneSymbol" = "SYMBOL"))
dim(LN_MA_967__LN_MS) ## 1
head(LN_MA_967__LN_MS)



c1 <- append(scRNA__LN_MS$gene, LN_MA_591__LN_MS$Gene.Symbol)
c2 <- append(c1, LN_MA_622__LN_MS$Gene.symbol)
MS_Overlap_DEGs <- unique( append(c2, LN_MA_967__LN_MS$GeneSymbol))
length(MS_Overlap_DEGs) ## 16


CD <- matrix(nrow = 16, ncol = 5)

for(i in 1:length(MS_Overlap_DEGs)){
  
  gene <- MS_Overlap_DEGs[i]
  CD[i,1] <- gene
  
  CD[i,2] <- ifelse(gene %in%  scRNA__LN_MS$gene, as.character(scRNA__LN_MS[which(scRNA__LN_MS$gene == gene), 5]), NA )
  CD[i,3] <- ifelse(gene %in%  LN_MA_591__LN_MS$Gene.Symbol, as.character(LN_MA_591__LN_MS[which(LN_MA_591__LN_MS$Gene.Symbol == gene), 6]), NA ) 
  CD[i,4] <- ifelse(gene %in%  LN_MA_622__LN_MS$Gene.symbol, as.character(LN_MA_622__LN_MS[which(LN_MA_622__LN_MS$Gene.symbol == gene), 5]), NA)
  CD[i,5] <- ifelse(gene %in%  LN_MA_967__LN_MS$GeneSymbol, as.character(LN_MA_967__LN_MS[which(LN_MA_967__LN_MS$GeneSymbol == gene), 6]), NA)
  

}


CD <- as.data.frame(CD)

rownames(CD) <- CD$V1
CD <- CD[,-1]

colnames(CD) <- c("scRNA","LN_MA_967","LN_MA_591","LN_MA_622")
str(CD)

pheatmap(CD, na_col = "black", cluster_rows=FALSE, cluster_cols=FALSE)




## All MS result Venn -----------------------

AG_HC <- as.data.frame(read_excel("../SLE_MS/primary result/1.xlsx", sheet = "Sheet1"))
dim(AG_HC) ## 176

AG_LN <- as.data.frame(read_excel("../SLE_MS/primary result/2.xlsx", sheet = "Sheet1"))
dim(AG_LN) ## 142


C1Q_HC <- as.data.frame(read_excel("../SLE_MS/primary result/3.xlsx", sheet = "Sheet1"))
dim(C1Q_HC) ## 24
C1Q_LN <- as.data.frame(read_excel("../SLE_MS/primary result/4.xlsx", sheet = "Sheet1"))
dim(C1Q_LN) ## 48

venn(list(AG_HC = AG_HC$Accession, 
          AG_LN = AG_LN$Accession, 
          C1Q_HC = C1Q_HC$Accession,
          C1Q_LN = C1Q_LN$Accession
          ),
     rggplot = TRUE,
     zcolor = "style")


venn(list(AG = c(AG_HC$Accession, AG_LN$Accession),
          C1Q = c(C1Q_HC$Accession,C1Q_LN$Accession)

),
rggplot = TRUE,
zcolor = "style")

