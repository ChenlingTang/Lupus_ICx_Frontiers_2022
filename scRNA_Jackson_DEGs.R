#library(dplyr)
library(Seurat)
#library(patchwork)
#library(dplyr)
#library(cowplot)
#library(ggplot2)
#library(ggpubr)
#library(reshape2)
library(DESeq2)
library(stringr)


## for macOX
#setwd("~/OneDrive - University Of Houston/WU_LAB/Single-Cell_LN/scRNA_SLE_Jackson/")

## for Ubuntu
setwd("/home/cougarnet.uh.edu/ctang9/Desktop/scRNA_SLE_Jackson/")

load("./scRNA_Jackson_DGEs.RData")
save.image("./scRNA_Jackson_DGEs.RData")


## Read Meta file (Cell-Infor, cellType, Group)
meta <- as.data.frame(read.csv2("scRNA_Jackson_meta.cvs"))
head(meta)
meta_sub <- subset(meta, select =c("Ident","group") )
rownames(meta_sub) <- meta_sub$Ident
head(meta_sub)

## Read GeneMatrix (count)
CountMatrix <- read.csv2("scRNA_Jackson_GeneMatrix_Count.cvs")
colnames(CountMatrix) <- str_replace(colnames(CountMatrix),"\\.","-") 
rownames(CountMatrix) <- CountMatrix$X
CountMatrix <- CountMatrix[,-1]

colnames(CountMatrix)[1:10]
dim(CountMatrix) ## 23257 genes across 88606 cells
CountMatrix[1:10,1:10]


## Use DESeq2 to find DEGs ---------------------------------------
CountMatrix_1 <- CountMatrix + 1 
all(colnames(CountMatrix_1) == rownames(meta_sub))

DEGs_Jackson <- DESeqDataSetFromMatrix(countData = CountMatrix_1,
                                       colData = meta_sub,
                                       design = ~group)


DEGs_Jackson$group <- relevel(DEGs_Jackson$group, ref = "aHD")
DEGs_Jackson <- DESeq(DEGs_Jackson)
DEGs_Jackson <- results(DEGs_Jackson)

DEGs_Jackson_sign <- DEGs_Jackson[-which(is.na(DEGs_Jackson$padj)), ]
## p.adj < 0.05
DEGs_Jackson_sign <- DEGs_Jackson_sign[(DEGs_Jackson_sign$padj < 0.05),] 


## write cvs2
write.csv2(DEGs_Jackson_sign, file = "DEGs_Jackson_S_Genes.csv", quote = FALSE)

save.image("./scRNA_Jackson_DGEs.RData")


######################################################################################
### Using Seurat to find DEGs ---------------------------------------------------

Jackson <- readRDS(file = "./scRNA_SLE_Jackson.rds")

Jackson

DimPlot(Jackson, reduction = "umap")

Idents(Jackson)
