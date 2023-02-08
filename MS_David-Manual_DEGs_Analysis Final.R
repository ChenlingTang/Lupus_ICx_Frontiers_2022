setwd("~/OneDrive - University Of Houston/WU_LAB/SLE_MS/Publish_Data/")

load("~/OneDrive - University Of Houston/WU_LAB/SLE_MS/Publish_Data/MS_David-Manual_DEGs_Analysis_Final.RDate")
save.image("~/OneDrive - University Of Houston/WU_LAB/SLE_MS/Publish_Data/MS_David-Manual_DEGs_Analysis_Final.RDate")

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
library(RColorBrewer)
library(ggpolypath)
library(caret)
library(scales)

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi) ## install.packages(c("Rcpp")) if select not working

## read DAVID convert ID list -----------------------------
AG_cv <- as.data.frame(read_excel("AG_LN-ONLY.xlsx", sheet = "AG_LN-ONLY"))
table(duplicated(AG_cv$Manul_Gene)) ## Falase
dim(AG_cv) ## 52

C1Q_cv <- as.data.frame(read_excel("C1Q_LN-ONLY.xlsx", sheet = "C1Q_LN-ONLY"))
table(duplicated(C1Q_cv$Manul_Gene)) ## Falase
dim(C1Q_cv) ## 27

## Combine AG & C1Q -----------------------------------------------------

inner <- inner_join(AG_cv, C1Q_cv, by = c("Manul_Gene" = "Manul_Gene"))

inner ## KRT14; TTR

venn(list(AG = AG_cv$Manul_Gene, C1Q = C1Q_cv$Manul_Gene),
             rggplot = TRUE, zcolor = "style")

GO_List <- as.character(unique(c(AG_cv$Manul_Gene, C1Q_cv$Manul_Gene)))

length(GO_List) ##77

keytypes(org.Hs.eg.db)

DEGs_ENID <- select(org.Hs.eg.db, keys=GO_List, keytype = "SYMBOL", columns= c("GENENAME","ENTREZID"))

head(DEGs_ENID)

write.xlsx(DEGs_ENID, "Real Publish Data/ALl_77_DEGs.xlsx", sheetName = "ID",row.names = FALSE)


### GO -----------------------------------------------------------------------------------------------
dim(DEGs_ENID) ## 77

## MF
MF <- enrichGO(DEGs_ENID$ENTREZID, "org.Hs.eg.db", ont = "MF", keyType = "ENTREZID",
               pAdjustMethod = "BH", pvalueCutoff  = 0.05,
                readable = TRUE)
dim(MF) ## 3

MF_Sim <- simplify(MF, cutoff = 1, by = "p.adjust", select_fun = min)
dim(MF_Sim) ## 3

dotplot(MF_Sim, showCategory = 8)

MF
## CC
CC <- enrichGO(DEGs_ENID$ENTREZID, "org.Hs.eg.db", ont = "CC", keyType = "ENTREZID",
               pAdjustMethod = "BH", pvalueCutoff  = 0.05,
                readable = TRUE)
dim(CC) ##9

CC_Sim <- simplify(CC, cutoff = 1, by = "p.adjust", select_fun = min)
dim(CC_Sim) ## 9

dotplot(CC_Sim, showCategory = 9)

## BP
BP <- enrichGO(DEGs_ENID$ENTREZID, "org.Hs.eg.db", ont = "BP", keyType = "ENTREZID",
               pAdjustMethod = "BH", pvalueCutoff  = 0.05,
               readable = TRUE)
dim(BP) ##52

BP_Sim <- simplify(BP, cutoff = 0.5, by = "p.adjust", select_fun = min)
dim(BP_Sim) ## 11

dotplot(BP_Sim, showCategory = 16)

## KEGG

head(DEGs_ENID)
KEGG <- enrichKEGG(DEGs_ENID$SYMBOL, keyType = 'uniprot', organism = "hsa",
                pvalueCutoff  = 0.05, minGSSize = 1)

dim(KEGG) ##3

dotplot(KEGG, showCategory = 16)

##############################################################################################
##------------------------------ HeatMap -----------------------------------------------------
##############################################################################################

head(DEGs_ENID)

## Disease class
DC <- as.data.frame(read_excel("./Real Publish Data/ALl_77_DEGs.xlsx", sheet = "Disease_Class"))
head(DC)
DEGs_ENID <- left_join(DEGs_ENID, DC,  by = c("SYMBOL" = "ID"))

## Tissue
TS <- as.data.frame(read_excel("./Real Publish Data/ALl_77_DEGs.xlsx", sheet = "Up_Tissue"))
head(TS)
DEGs_ENID <- left_join(DEGs_ENID, TS,  by = c("SYMBOL" = "ID"))

## KEGG
KEGG <- as.data.frame(read_excel("./Real Publish Data/ALl_77_DEGs.xlsx", sheet = "KEGG"))
DEGs_ENID <- left_join(DEGs_ENID, KEGG,  by = c("SYMBOL" = "ID"))


dim(DEGs_ENID) ## 77
colnames(DEGs_ENID)

## remove keratin 
dim(DEGs_ENID[grep("keratin", DEGs_ENID$GENENAME) ,]) ## 10
DEGs_1 <- DEGs_ENID[-grep("keratin", DEGs_ENID$GENENAME) ,]

## remove immunoglobulin  --16 
dim(DEGs_ENID[grep("immunoglobulin", DEGs_ENID$GENENAME) ,]) ##16
DEGs_1 <- DEGs_1[-grep("immunoglobulin", DEGs_1$GENENAME) ,]

dim(DEGs_1) ## 51 = 77-10-16

## Disease class  -----------------------------
## Key words: IMMUNE; RENAL

table(unlist(strsplit(unlist(DEGs_1$GAD_DISEASE_CLASS),",")))

DC_Pos <- unique(c(grep("IMMUNE", DEGs_1$GAD_DISEASE_CLASS),
                   grep("RENAL", DEGs_1$GAD_DISEASE_CLASS)))

DEGs_1$DC_POS <- "N"
DEGs_1[DC_Pos,]$DC_POS <- "Y"

colnames(DEGs_1)

## Tissue  -----------------------------
## Key words: Skin; muscle;bone; joints; Kidney; Spleen; Liver; Brain; 

table(unlist(strsplit(unlist(DEGs_1$UP_TISSUE),",")))

TS_Pos <- unique(c(grep("Skin", DEGs_1$UP_TISSUE, ignore.case = TRUE),
                   
                   grep("muscle", DEGs_1$UP_TISSUE, ignore.case = TRUE),
                   grep("bone", DEGs_1$UP_TISSUE, ignore.case = TRUE),
                   grep("joints", DEGs_1$UP_TISSUE, ignore.case = TRUE),
                   
                   grep("Kidney", DEGs_1$UP_TISSUE, ignore.case = TRUE),
                   grep("Spleen", DEGs_1$UP_TISSUE, ignore.case = TRUE),
                   grep("Liver", DEGs_1$UP_TISSUE, ignore.case = TRUE),
                   grep("Lymph", DEGs_1$UP_TISSUE, ignore.case = TRUE),
                   grep("Lymphocyte", DEGs_1$UP_TISSUE, ignore.case = TRUE),
                   grep("Brain", DEGs_1$UP_TISSUE, ignore.case = TRUE)
))

length(TS_Pos)

DEGs_1$TS_Pos <- "N"
DEGs_1[TS_Pos,]$TS_Pos <- "Y"

DEGs_1[DEGs_1$SYMBOL =="DOCK11",]

## KEGG  -----------------------------
## SLE (hsa05322) related pathway: hsa04060;hsa04514;hsa04610;hsa04612;hsa04630;hsa04660;hsa04662;hsa04670

table(unlist(strsplit(unlist(DEGs_1$KEGG_PATHWAY),",")))

TS_KEGG <- unique(c(grep("hsa05322", DEGs_1$KEGG_PATHWAY, ignore.case = TRUE),
                    grep("hsa04060", DEGs_1$KEGG_PATHWAY, ignore.case = TRUE),
                    grep("hsa04514", DEGs_1$KEGG_PATHWAY, ignore.case = TRUE),
                    grep("hsa04610", DEGs_1$KEGG_PATHWAY, ignore.case = TRUE),
                    grep("hsa04612", DEGs_1$KEGG_PATHWAY, ignore.case = TRUE),
                    grep("hsa04630", DEGs_1$KEGG_PATHWAY, ignore.case = TRUE),
                    grep("hsa04660", DEGs_1$KEGG_PATHWAY, ignore.case = TRUE),
                    grep("hsa04662", DEGs_1$KEGG_PATHWAY, ignore.case = TRUE),
                    grep("hsa04670", DEGs_1$KEGG_PATHWAY, ignore.case = TRUE)
                    
))

DEGs_1$KEGG_Pos <- "N"
DEGs_1[TS_KEGG,]$KEGG_Pos <- "Y"

DEGs_1[TS_KEGG,]$KEGG_PATHWAY


## concise the data -------------------------------------------------------------------
head(DEGs_1)

Pos_index <- grep("Pos", colnames(DEGs_1), ignore.case = TRUE)

colnames(DEGs_1)

DEGs_1s <- subset(DEGs_1, select = c(1,Pos_index))

colnames(DEGs_1s)[1] <- "Gene"
head(DEGs_1s)

## Immune Response related HeatMap --------------------------------------------------------------------

BP_David <- as.data.frame(read_excel("./Real Publish Data/ALl_77_DEGs.xlsx", sheet = "BP"))
colnames(BP_David)

#BP_David[grep("immune response", BP_David$GOTERM_BP_DIRECT, ignore.case = TRUE),]

ISR_set <- BP_David[grep("immune response", BP_David$GOTERM_BP_DIRECT, ignore.case = TRUE),
              which(colnames(BP_David) == "ID")]

ISR_genes <- unique(unlist(str_split(ISR_set,"/")))

DEGs_1s$BP_ImmuneResponseRelated <- "N"

for( i in 1:length(ISR_genes)){
  gene <- ISR_genes[i]
  DEGs_1s[grep(gene, DEGs_1s$Gene), 
          which(colnames(DEGs_1s) == "BP_ImmuneResponseRelated")] <- "Y"
}

table(DEGs_1s$BP_ImmuneResponseRelated) ## Y-7

## B cell activation

# BCA_set <-  BP_David[grep("B cell activation", BP_David$GOTERM_BP_DIRECT, ignore.case = TRUE),
#                      which(colnames(BP_David) == "ID")]

# BCA_genes <- unique(unlist(str_split(BCA_set,"/")))
# 
# DEGs_1s$BP_BCellActivation <- "N"
# for( i in 1:length(BCA_genes)){
#   gene <- BCA_genes[i]
#   DEGs_1s[grep(gene, DEGs_1s$Gene), 
#           which(colnames(DEGs_1s) == "BP_BCellActivation")] <- "Y"
# }
# 
# table(DEGs_1s$BP_BCellActivation) ## Y-0


## Inflammatory

INF_set <-  BP_David[grep("Inflammatory", BP_David$GOTERM_BP_DIRECT, ignore.case = TRUE),
                     which(colnames(BP_David) == "ID")]

INF_genes <- unique(unlist(str_split(INF_set,"/")))

DEGs_1s$BP_Inflammatory <- "N"
for( i in 1:length(INF_genes)){
  gene <- INF_genes[i]
  DEGs_1s[grep(gene, DEGs_1s$Gene), 
          which(colnames(DEGs_1s) == "BP_Inflammatory")] <- "Y"
}

table(DEGs_1s$BP_Inflammatory) ## Y-6


## glomerular filtration

 BP_David[grep("glomerular", BP_David$GOTERM_BP_DIRECT, ignore.case = TRUE),
          ]
# GLO_set <-  BP_David[grep("glomerular", BP_David$GOTERM_BP_DIRECT, ignore.case = TRUE),
#                      which(colnames(BP_David) == "ID")]

# GLO_genes <- unique(unlist(str_split(GLO_set,"/")))
# 
# DEGs_1s$BP_GlomerularFiltration <- "N"
# for( i in 1:length(GLO_genes)){
#   gene <- GLO_genes[i]
#   DEGs_1s[grep(gene, DEGs_1s$Gene), 
#           which(colnames(DEGs_1s) == "BP_GlomerularFiltration")] <- "Y"
# }
# 
# table(DEGs_1s$BP_GlomerularFiltration) ## All N

# ## complement activation
# 
# CA_set <-  BP_David[grep("complement activation", BP_David$GOTERM_BP_DIRECT, ignore.case = TRUE),
#                      which(colnames(BP_David) == "ID")]
# 
# CA_genes <- unique(unlist(str_split(CA_set,"/")))
# 
# DEGs_1s$BP_ComplementActiviation <- "N"
# for( i in 1:length(CA_genes)){
#   gene <- CA_genes[i]
#   DEGs_1s[grep(gene, DEGs_1s$Gene), 
#           which(colnames(DEGs_1s) == "BP_ComplementActiviation")] <- "Y"
# }
# 
# table(DEGs_1s$BP_ComplementActiviation) ## Y-1
# 

## CC enrich Gene -----------------------------------------------
CC_David <- as.data.frame(read_excel("./Real Publish Data/ALl_77_DEGs.xlsx", sheet = "CC"))
colnames(CC_David)


## immunoglobulin complex

CC_David[grep("immunoglobulin", CC_David$GOTERM_CC_DIRECT, ignore.case = TRUE),
       which(colnames(CC_David) == "ID")]

IC_set <- CC_David[grep("immunoglobulin complex", CC_David$GOTERM_CC_DIRECT, ignore.case = TRUE),
                   which(colnames(CC_David) == "ID")]

IC_genes <- unique(unlist(str_split(IC_set,"/"))) ## All Ig


## MF

MF_David <- as.data.frame(read_excel("./Real Publish Data/ALl_77_DEGs.xlsx", sheet = "MF"))
colnames(MF_David)

## antigen binding

AB_set <- MF_David[grep("antigen binding", MF_David$GOTERM_MF_DIRECT, ignore.case = TRUE),
                 which(colnames(MF_David) == "ID")]

AB_genes <- unique(unlist(str_split(AB_set,"/")))  ## All Ig


head(DEGs_1s)



## Comparing with Significant gene from Database (p<0.05) ----------------------------------------------

DEGs_ALL <- DEGs_1s
dim(DEGs_ALL)
head(DEGs_ALL)

## read scRNA Albert lupus   ----------------------------------------------------------------------------
## Der_2019_scRNA-Seq_LN_Kidney, https://www.nature.com/articles/s41590-019-0529-4
## DESeq2 found DEGs
#scRNA_Alb_DESeq2 <- as.data.frame(read.csv2("DataBase/scRNA_Albert_Kidney_DESeq2.csv",sep =","))
#head(scRNA_Alb_DESeq2)
#dim(scRNA_Alb_DESeq2) ##8894

#scRNA_Alb_DESeq2_s <- subset(scRNA_Alb_DESeq2, select = c("Gene", "log2FoldChange"))
#dim(scRNA_Alb_DESeq2_s)
#scRNA_Alb_DESeq2_s <- scRNA_Alb_DESeq2_s[!duplicated(scRNA_Alb_DESeq2_s$Gene),]
#colnames(scRNA_Alb_DESeq2_s)[2] <- "LGC_scRNA_Albert_lupus_kidney"

#DEGs_ALL <- left_join(DEGs_ALL, scRNA_Alb_DESeq2_s, by = c("Gene" = "Gene"))
#dim(DEGs_ALL)

## DESingle found DEGs
scRNA_Alb_DESingle <- as.data.frame(read.csv2("DataBase/DEGs_Albert_S_Genes_DESignle.csv",sep =";"))
head(scRNA_Alb_DESingle)
dim(scRNA_Alb_DESingle) ##4432

scRNA_Alb_DESingle$log2FoldChange <- log2(scRNA_Alb_DESingle$norm_total_mean_2 / scRNA_Alb_DESingle$norm_total_mean_1)

scRNA_Alb_DESingle_s <- subset(scRNA_Alb_DESingle, select = c("Gene", "log2FoldChange"))
dim(scRNA_Alb_DESingle_s)

colnames(scRNA_Alb_DESingle_s)[2] <- "Der_2019_scRNA-Seq_LN_Kidney"

DEGs_ALL <- left_join(DEGs_ALL, scRNA_Alb_DESingle_s, by = c("Gene" = "Gene"))
dim(DEGs_ALL)
head(DEGs_ALL)
## read scRNA Jackson lupus Serum  ----------------------------------------------------------------------------
## Nehar-Belaid_2020_scRNA-Seq_SLE_PBMC
## DESeq2
# scRNA_Jackson <- as.data.frame(read.csv2("DataBase/DEGs_Jackson_S_Genes2.csv",sep =","))
# head(scRNA_Jackson)
# dim(scRNA_Jackson) ##3157
# scRNA_Jackson$Gene <- rownames(scRNA_Jackson)
# 
# scRNA_JacksonS<- subset(scRNA_Jackson, select = c("Gene", "log2FoldChange"))
# head(scRNA_JacksonS)
# colnames(scRNA_JacksonS)[2] <- "LGC_scRNA_Jackson_SLE_Serum"
# 
# DEGs_ALL <- left_join(DEGs_ALL, scRNA_JacksonS, by = c("Gene" = "Gene"))
# dim(DEGs_ALL) 
# head(DEGs_ALL)

## DESignle
scRNA_Jackson_DESingle <- as.data.frame(read.csv2("DataBase/DEGs_Jackson_S_Genes_DESignle.csv",sep =";"))
colnames(scRNA_Jackson_DESingle)[1] <- "Gene" 
dim(scRNA_Jackson_DESingle) ##8406
head(scRNA_Jackson_DESingle)

scRNA_Jackson_DESingle$log2FoldChange <- log2(scRNA_Jackson_DESingle$total_mean_2 / scRNA_Jackson_DESingle$total_mean_1)
scRNA_Jackson_DESingleS<- subset(scRNA_Jackson_DESingle, select = c("Gene", "log2FoldChange"))
head(scRNA_Jackson_DESingleS)
colnames(scRNA_Jackson_DESingleS)[2] <- "Nehar-Belaid_2020_scRNA-Seq_SLE_PBMC"


DEGs_ALL <- left_join(DEGs_ALL, scRNA_Jackson_DESingleS, by = c("Gene" = "Gene"))
dim(DEGs_ALL) 
head(DEGs_ALL)

## read RNA Zhengzhou lupus   ----------------------------------------------------------------------------
## Yao_2020_RNA-Seq_LN_Kidney,https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157293
RNA_ZZ <- as.data.frame(read.csv2("DataBase/GSE157293_RNA_S_Gene.csv",sep =","))
head(RNA_ZZ)
dim(RNA_ZZ) ##2247

RNA_ZZ_s <- subset(RNA_ZZ, select = c("Gene", "log2FoldChange"))
head(RNA_ZZ_s)
colnames(RNA_ZZ_s)[2] <- "Yao_2020_RNA-Seq_LN_Kidney"


DEGs_ALL <- left_join(DEGs_ALL, RNA_ZZ_s, by = c("Gene" = "Gene"))
head(DEGs_ALL)

## read RNA Imperial SLE   ----------------------------------------------------------------------------
## Buang_2021_RNA-Seq_SLE_CD8+TCell,https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97264
RNA_IC <- read.csv2("DataBase/GSE97264_HC_LN_S_Gene.csv",sep ="\t")
head(RNA_IC)
dim(RNA_IC) ##2312

RNA_IC_s <- subset(RNA_IC, select = c("SYMBOL", "GSE97264_LN_HC_log2FC"))
head(RNA_IC_s)
colnames(RNA_IC_s)[2] <- "Buang_2021_RNA-Seq_SLE_CD8+TCell"


DEGs_ALL <- left_join(DEGs_ALL, RNA_IC_s, by = c("Gene" = "SYMBOL"))
head(DEGs_ALL)

## Microarray, Xiangya, SLE, GSE81622 --------------------------------------------------------
## Zhu_2016_Expression-Array_SLE_PBMC  ,https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
MA_xiangya_SLE <- read.csv2("DataBase/GSE81622_S_Gene.csv",sep =",")
head(MA_xiangya_SLE)
dim(MA_xiangya_SLE) ##1638

MA_xiangya_SLE_s <- subset(MA_xiangya_SLE, select = c("Gene.symbol", "logFC"))
head(MA_xiangya_SLE_s)
colnames(MA_xiangya_SLE_s)[2] <- "Zhu_2016_Expression-Array_SLE_PBMC"


DEGs_ALL <- left_join(DEGs_ALL, MA_xiangya_SLE_s, by = c("Gene" = "Gene.symbol"))
head(DEGs_ALL)

## Microarray, Toronto, LN, GSE99967 --------------------------------------------------------
##  Wither_2018_Expression-Array_LN_Blood ,https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
# MA_Tornoto_LN <- read.csv2("DataBase/GSE99967-Microarray_S_Gene.csv",sep =",")
# head(MA_Tornoto_LN)
# dim(MA_Tornoto_LN) ##1538
# 
# MA_Tornoto_LN_s <- subset(MA_Tornoto_LN, select = c("GeneSymbol", "logFC"))
# head(MA_Tornoto_LN_s)
# colnames(MA_Tornoto_LN_s)[2] <- "Wither_2018_Expression-Array_LN_Blood"
# 
# 
# DEGs_ALL <- left_join(DEGs_ALL, MA_Tornoto_LN_s, by = c("Gene" = "GeneSymbol"))
# head(DEGs_ALL)

## Microarray, MichiganU, LN, GSE32591 --------------------------------------------------------
##  Berthier_2012_Expression-Array_LN_Kidney  ,https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32591
MA_MU_LN <- read.csv2("DataBase/GSE32591-Microarray_S_Gene.csv",sep =",")
head(MA_MU_LN)
dim(MA_MU_LN) ##5823

MA_MU_LN_s <- subset(MA_MU_LN, select = c("Gene.Symbol", "logFC"))
head(MA_MU_LN_s)
colnames(MA_MU_LN_s)[2] <- "Berthier_2012_Expression-Array_LN_Kidney"


DEGs_ALL <- left_join(DEGs_ALL, MA_MU_LN_s, by = c("Gene" = "Gene.Symbol"))
head(DEGs_ALL)
dim(DEGs_ALL)


## ANtigen Array, UTSW, SLE --------------------------------------------------------

## 437 IgG AutoAbs
# AG_UTSW_SLE <- as.data.frame(read_xlsx("DataBase/DEGs_UTSW_AAB_SLE.xlsx", sheet = "sheet 1", skip =1))
# head(AG_UTSW_SLE)
# AG_UTSW_SLE$logFC <- log2(AG_UTSW_SLE$`Mean SLE` / AG_UTSW_SLE$`Mean NC`)
# dim(AG_UTSW_SLE) ## 437
# head(AG_UTSW_SLE)
#
# AG_UTSW_SLE_s <- subset(AG_UTSW_SLE, select = c("Protein target", "logFC"))
# head(AG_UTSW_SLE_s)
# colnames(AG_UTSW_SLE_s)[2] <- "LGC_AgArray_AutoIgG_UTSW_SLE_Serum"
# 
# 
# DEGs_ALL <- left_join(DEGs_ALL, AG_UTSW_SLE_s, by = c("Gene" = "Protein target"))
# head(DEGs_ALL)
# dim(DEGs_ALL)
# 
# ## 1213 IgM AutoAbs
# AG_UTSW_SLE_IgM <- as.data.frame(read_xlsx("DataBase/DEGs_UTSW_AAB_SLE.xlsx", sheet = "sheet 2", skip =1))
# head(AG_UTSW_SLE_IgM)
# AG_UTSW_SLE_IgM$logFC <- log2(AG_UTSW_SLE_IgM$`Mean SLE` / AG_UTSW_SLE_IgM$`Mean NC`)
# dim(AG_UTSW_SLE_IgM) ## 1213
# head(AG_UTSW_SLE_IgM)
# 
# AG_UTSW_SLE_IgM_s <- subset(AG_UTSW_SLE_IgM, select = c("Protein target", "logFC"))
# head(AG_UTSW_SLE_IgM_s)
# colnames(AG_UTSW_SLE_IgM_s)[2] <- "LGC_AgArray_AutoIgM_UTSW_SLE_Serum"
# 
# 
# DEGs_ALL <- left_join(DEGs_ALL, AG_UTSW_SLE_IgM_s, by = c("Gene" = "Protein target"))
# head(DEGs_ALL)
# dim(DEGs_ALL)


### Make Heatmap ---------------------------
DEGs_ALL_HM <- DEGs_ALL

head(DEGs_ALL_HM)


DEGs_ALL_HM[DEGs_ALL_HM == "Y"] <- 5
DEGs_ALL_HM[DEGs_ALL_HM == "N"] <- NA

DEGs_ALL_HM[DEGs_ALL_HM == "Inf"] <- 5

DEGs_ALL_HM[,-1] <- as.data.frame(apply(DEGs_ALL_HM[,-1] , 2, as.numeric))

rownames(DEGs_ALL_HM) <- DEGs_ALL_HM$Gene


## add SLE score in MS 
head(DEGs_ALL)
DEGs_ALL_MS <- subset(DEGs_ALL, select= c("Gene"))

DEGs_ALL_MS$label <- ifelse(DEGs_ALL_MS$Gene %in% AG_cv$Manul_Gene, "Protein A/G", "C1Q")

DEGs_ALL_MS[DEGs_ALL_MS$Gene == "TTR", ]$label <- "Both"

table(DEGs_ALL_MS$label)

DEGs_ALL_MS <- DEGs_ALL_MS[c(1:20,22:30,21,31:51),]

## AG
DEGs_ALL_MS_AG <- DEGs_ALL_MS[DEGs_ALL_MS$label == "Protein A/G",]

DEGs_ALL_MS_AG <- left_join(DEGs_ALL_MS_AG, AG_cv, by = c("Gene" = "Manul_Gene"))

DEGs_ALL_MS_AG <- subset(DEGs_ALL_MS_AG, select= c("Gene","label","SLE Score"))

DEGs_ALL_MS_AG <- DEGs_ALL_MS_AG[order(DEGs_ALL_MS_AG$`SLE Score`, decreasing = TRUE),]

## C1Q
DEGs_ALL_MS_C1Q <- DEGs_ALL_MS[DEGs_ALL_MS$label == "C1Q",]

DEGs_ALL_MS_C1Q <- left_join(DEGs_ALL_MS_C1Q, C1Q_cv, by = c("Gene" = "Manul_Gene"))

DEGs_ALL_MS_C1Q <- subset(DEGs_ALL_MS_C1Q, select= c("Gene","label","SLE Score"))

DEGs_ALL_MS_C1Q <- DEGs_ALL_MS_C1Q[order(DEGs_ALL_MS_C1Q$`SLE Score`, decreasing = TRUE),]

## Both 
DEGs_ALL_MS_Both <- DEGs_ALL_MS[DEGs_ALL_MS$label == "Both",]

TTR_AG <- AG_cv[AG_cv$Manul_Gene == "TTR",]$`SLE Score`
TTR_C1Q <- C1Q_cv[C1Q_cv$Manul_Gene == "TTR",]$`SLE Score`

DEGs_ALL_MS_Both$"SLE Score" <- mean(TTR_AG,TTR_C1Q)
  
DEGs_ALL_MS_score <- as.data.frame(
  rbind(
    rbind(DEGs_ALL_MS_AG, DEGs_ALL_MS_Both),
        DEGs_ALL_MS_C1Q))

## remove "housekeep gene" or structure protein
DEGs_ALL_MS_score <- DEGs_ALL_MS_score[-which(DEGs_ALL_MS_score$Gene %in% c("SCFV", "CFH",
                                                                       "ALB","H2AC12","H2BC5","C3")),]

DEGs_ALL_MS_score$`SLE ScoreN` <- rescale(DEGs_ALL_MS_score$`SLE Score`, to= c(0,5) )
head(DEGs_ALL_MS_score)
table(DEGs_ALL_MS_score$label)


## combine 
DEGs_ALL_Combine <- left_join(DEGs_ALL_MS_score,DEGs_ALL_HM, by = c("Gene" = "Gene") )

rownames(DEGs_ALL_Combine) <- DEGs_ALL_Combine$Gene

DEGs_ALL_Final <- DEGs_ALL_Combine[,-c(1,2,3)]

colnames(DEGs_ALL_Final)

DEGs_ALL_Final <- DEGs_ALL_Final[,c(2:6,1,7:12)]
colnames(DEGs_ALL_Final)

colnames(DEGs_ALL_Final)[1] <- "Disease-Class_DAVID"
colnames(DEGs_ALL_Final)[2] <- "Up-Tissue_DAVID"
colnames(DEGs_ALL_Final)[3] <- "SLE-Related-Pathway_KEGG"

colnames(DEGs_ALL_Final)[4] <- "ImmuneResponseRelated_GO-BP"
colnames(DEGs_ALL_Final)[5] <- "Inflammatory_GO-BP"

colnames(DEGs_ALL_Final)[6] <- "LN MS Score"

colnames(DEGs_ALL_Final)
str(DEGs_ALL_Final)

### Plot the heatmap -----------------------------------------------------------------------
paletteLength <- 44
myColor <- colorRampPalette(c("#1100EE", "white", "orangered"))(paletteLength)
myBreaks <- c(seq(-2, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(5/paletteLength, max(5), length.out=floor(paletteLength/2)))
pheatmap(DEGs_ALL_Final, color=myColor, breaks=myBreaks, 
         cluster_cols = FALSE, cluster_rows = FALSE, na_col = "white",
         gaps_col = c(5,6), gaps_row = c(26,27)
        )
### Find the overlap gene across the six database ------------------------------------------

venn(list(scRNASeq_LN_Kidney = scRNA_Alb_DESingle_s$Gene,
                  scRNASeq_SLE_PBMC= scRNA_Jackson_DESingleS$Gene,
                  RNASeq_LN_Kidney = RNA_ZZ_s$Gene,
                  RNASeq_SLE_CD8TCell = RNA_IC_s$SYMBOL,
                  ExprArray_SLE_PBMC = MA_xiangya_SLE_s$Gene.symbol,
                  ExprArray_LN_Kidney = MA_MU_LN_s$Gene.Symbol
        ),
     ggplot = TRUE,
     zcolor = "style"
     )


SixDatabase <- inner_join(scRNA_Alb_DESingle_s, scRNA_Jackson_DESingleS, by = c("Gene" = "Gene") )
SixDatabase <- inner_join(SixDatabase, RNA_ZZ_s, by = c("Gene" = "Gene") )
SixDatabase <- inner_join(SixDatabase, RNA_IC_s, by = c("Gene" = "SYMBOL") )
SixDatabase <- inner_join(SixDatabase, MA_xiangya_SLE_s, by = c("Gene" = "Gene.symbol") )
SixDatabase <- inner_join(SixDatabase, MA_MU_LN_s, by = c("Gene" = "Gene.Symbol") )

SixDatabase_uniq <- SixDatabase[!duplicated(SixDatabase$Gene),]

SixDatabase_uniq[SixDatabase_uniq == "Inf"] <- 9
SixDatabase_uniq[,-1] <- as.data.frame(apply(SixDatabase_uniq[,-1] , 2, as.numeric))

rownames(SixDatabase_uniq) <- SixDatabase_uniq$Gene
SixDatabase_uniq <- SixDatabase_uniq[,-1]

str(SixDatabase_uniq)

pheatmap(SixDatabase_uniq, cluster_cols = FALSE)


## sub cellular localization
head(AG_cv)
head(C1Q_cv)

AccList <- unique(c(AG_cv$Accession, C1Q_cv$Accession))
length(AccList) ## 77

library(UniprotR)
SubCellLoc <- GetSubcellular_location(AccList)
colnames(SubCellLoc)
SubCellLoc_sub <- SubCellLoc$Subcellular.location..CC.


head(DEGs_ENID)
COMP_data <- read.table(file = "human_compartment_integrated_full.tsv", sep = "\t")
#COMP_data <- read.table(file = "human_compartment_experiments_full.tsv", sep = "\t")

head(COMP_data)

COMP_list <- as.data.frame(DEGs_ENID$SYMBOL)
colnames(COMP_list)[1] <- "DEGs"
head(COMP_list)

COMP_CellLoc <- right_join(COMP_data, COMP_list, by = c("V2" = "DEGs"))
#colnames(COMP_CellLoc) <- c("ID","Name","GO","CellLoc","Score")
length(unique(COMP_CellLoc$V2)) ## 77
table(COMP_CellLoc$V2)

head(COMP_CellLoc)



