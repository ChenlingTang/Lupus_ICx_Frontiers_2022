## Packages 
{
library(dplyr)
library(Seurat)
library(patchwork)
library(dplyr)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(DESeq2)
library(SingleR)
library(celldex)
library(stringr)
library(DEsingle)
}

## for macOX
#setwd("~/OneDrive - University Of Houston/WU_LAB/Single-Cell_LN/scRNA_SLE_Jackson/")

## for Ubuntu
setwd("/home/cougarnet.uh.edu/ctang9/Desktop/scRNA_SLE_Jackson/")

load("./scRNA_SLE_Jackson_SingleR.RData")
save.image("./scRNA_SLE_Jackson_SingleR.RData")


## scRNA SLE adult from Jackson ----------------------------------------------------------
## Paper link: https://www.nature.com/articles/s41590-020-0743-0#Sec15
## Analysis follow Seurat tuition
## Tuition Link: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# Load the PBMC dataset -----------------------------------------------------
pbmc.data <- Read10X(data.dir = "./")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "SLE_Jackson", min.cells = 3, min.features = 200)

pbmc
## 23257 genes across 88606 cells


## Pre-Processing ----------------------------------------------------------
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## In the original Paper
## all genes that were not detected in ≥3 cells were discarded
## 400 < nFeature < 2500
## cells in which >20% of the transcripts mapped to the mitochondrial genes were filtered out
pbmc <- subset(pbmc, subset = nFeature_RNA > 400 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc
## 23257 gene across 82748 cells

## Normalizing the data
## log-tranformed and 10000 scale factor
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


## feature selection
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes) ## running on server only
# pbmc <- ScaleData(pbmc) ## running on MacOS

## Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

## Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100) ## only on server
#pbmc <- JackStraw(pbmc, num.replicate = 30) ## On macOS
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:20)
ElbowPlot(pbmc) ## 12 PC


## Cluster the Cells
pbmc <- FindNeighbors(pbmc, dims = 1:18) ## dims from previous step
pbmc <- FindClusters(pbmc, resolution = 0.3) ## 15 clusters
levels(pbmc)

table(pbmc$seurat_clusters)

## Run non-linear dimensional reduction
## UMAP
pbmc <- RunUMAP(pbmc, dims = 1:18) ## same PCs as cluster
DimPlot(pbmc, reduction = "umap", label = TRUE)

## TSNE
pbmc <- RunTSNE(pbmc, dims = 1:18)
DimPlot(pbmc, reduction = "tsne",label = TRUE)

save.image("./scRNA_SLE_Jackson_SingleR.RData")

# ## Finding differentially expressed features
# # find markers for every cluster compared to all remaining cells, report only the positive ones
# pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top10
# dim(top10)
# DoHeatmap(pbmc, features = top10$gene) + NoLegend()



canonial_markers <-  c("MS4A1", ## B
                       
                       "GNLY","NKG7", ## NK

                       "CD8A", ## CD8+ T
                       "IL7R", "CCR7", ## Naive CD4+ T
                       #"S100A4", #"IL7R",## Memory CD4+ T
                       "PTPRC", ## general immune cell
                       "FCER1A","CST3", ## DCs
                       
                       "CD14","LYZ", ## CD14+ Mono
                       "FCGR3A","MS4A7", ## FCGR3A+ Mono
                       
                       
                       "PPBP" ## Platelet
                       
                      )

# VlnPlot(pbmc, features = canonial_markers )
# VlnPlot(pbmc, features = canonial_markers, slot = "counts", log = TRUE)

VlnPlot(pbmc, features =c("FCER1A","CST3")) ## cluster 10 is DC

cells_DCs <- WhichCells(pbmc, idents = '10')
length(cells_DCs)

FeaturePlot(pbmc, features = canonial_markers)


# new.cluster.ids <- c("Memory CD4+ T", #0
#                       "Naive CD4+ T", #1 
#                      "CD8+ T",  #2
#                      "CD14+ Mono", #3
#                      "Naive CD4+ T", #4
#                      "NK", #5
#                      "CD14+ Mono", #6
#                      "B Cell", #7
#                      "FCGR3A+ Mono", #8
#                      "Platelet", #9
#                      "plasmacytoid DCs", #10
#                      "Erythroid Cell", #11
#                       "Plasma Cell", #12
#                      "DCs" #13
#                      )
# 
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# 
# DotPlot(pbmc, features = canonial_markers) + RotatedAxis()
# 
# DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) 
# DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) 



###############################################################################################
## SingleR Annotation -------------------------------------------------------------------
# Loading reference data with Ensembl annotations.

## UMAP
pbmc <- RunUMAP(pbmc, dims = 1:18) ## same PCs as cluster
DimPlot(pbmc, reduction = "umap", label = TRUE)

## TSNE
pbmc <- RunTSNE(pbmc, dims = 1:18)
DimPlot(pbmc, reduction = "tsne",label = TRUE)

## PrimaryCell ----------------------------------------------------------
ref.data <- HumanPrimaryCellAtlasData()

table(ref.data$label.main)

length(intersect(rownames(pbmc) ,rownames(ref.data)))

# Performing predictions.
labels <- factor(ref.data$label.main)

predictions <- SingleR(test = as.SingleCellExperiment(pbmc), 
                       ref=ref.data, labels= labels, fine.tune = TRUE)

table(predictions$labels)
# B_cell               BM              CMP               DC     Erythroblast              GMP       HSC_-G-CSF              MEP 
# 3526                4               69                5                3               37               55                3 
# Monocyte      Neutrophils          NK_cell        Platelets Pre-B_cell_CD34- Pro-B_cell_CD34+          T_cells 
# 18956               17             7849              224              729                9            51262 
pbmc$SingleR_Primary.labels<- predictions$labels

Idents(object = pbmc) <- pbmc$SingleR_Primary.labels

cells_Platelets <- WhichCells(pbmc, idents = 'Platelets')

DimPlot(pbmc, reduction = "umap" , label = TRUE, pt.size = 0.5)# + NoLegend()
DimPlot(pbmc, reduction = "tsne" , label = TRUE, pt.size = 0.5)# + NoLegend()


## Immune Cell ------------------------------------------
ref_imm <- DatabaseImmuneCellExpressionData()

table(ref_imm$label.fine)

length(intersect(rownames(pbmc) ,rownames(ref_imm)))

# Performing predictions
labels_imm <- factor(ref_imm$label.fine)

predictions_Immune <- SingleR(test = as.SingleCellExperiment(pbmc), 
                       ref=ref_imm, labels= labels_imm, fine.tune = TRUE)

table(predictions_Immune$labels)
# B cells, naive                 Monocytes, CD14+                 Monocytes, CD16+                         NK cells 
# 3543                            17111                             2818                            15789 
# T cells, CD4+, memory TREG             T cells, CD4+, naive        T cells, CD4+, naive TREG               T cells, CD4+, TFH 
# 2752                             8286                             1496                             9240 
# T cells, CD4+, Th1            T cells, CD4+, Th1_17              T cells, CD4+, Th17               T cells, CD4+, Th2 
# 5383                             1694                             3934                             2207 
# T cells, CD8+, naive T cells, CD8+, naive, stimulated 
# 8492                                3 

pbmc$SingleR_Immune.labels<- predictions_Immune$labels

Idents(object = pbmc) <- pbmc$SingleR_Immune.labels

## Set platelets and DCs cell type
pbmc <- SetIdent(pbmc, cells = cells_Platelets, value = 'Platelets')
pbmc <- SetIdent(pbmc, cells = cells_DCs, value = 'DC')

## remove redundant cell type
pbmc <- RenameIdents(object = pbmc, `T cells, CD8+, naive, stimulated` = "T cells, CD8+")
pbmc <- RenameIdents(object = pbmc, `T cells, CD8+, naive` = "T cells, CD8+")

pbmc <- RenameIdents(object = pbmc, `T cells, CD4+, Th1_17` = "T cells, CD4+")
pbmc <- RenameIdents(object = pbmc, `T cells, CD4+, naive TREG` = "T cells, CD4+")
pbmc <- RenameIdents(object = pbmc, `T cells, CD4+, memory TREG` = "T cells, CD4+")
pbmc <- RenameIdents(object = pbmc, `T cells, CD4+, TFH` = "T cells, CD4+")
pbmc <- RenameIdents(object = pbmc, `T cells, CD4+, naive` = "T cells, CD4+")
pbmc <- RenameIdents(object = pbmc, `T cells, CD4+, Th17` = "T cells, CD4+")
pbmc <- RenameIdents(object = pbmc, `T cells, CD4+, Th1` = "T cells, CD4+")
pbmc <- RenameIdents(object = pbmc, `T cells, CD4+, Th2` = "T cells, CD4+")

pbmc <- RenameIdents(object = pbmc, `B cells, naive` = "B cells")

table(Idents(pbmc))
# T cells, CD4+    T cells, CD8+               DC        Platelets         NK cells Monocytes, CD14+ Monocytes, CD16+   B cells 
# 34983             8495              514              224            15732            16494             2764             3542 

Idents(pbmc) <- factor(Idents(pbmc), levels = c("T cells, CD4+","T cells, CD8+", 
                                                "B cells", 
                                                "NK cells",
                                                "Monocytes, CD14+","Monocytes, CD16+",
                                                "DC",
                                                "Platelets"
                ))

DimPlot(pbmc, reduction = "umap" , label = TRUE, pt.size = 0.5)# + NoLegend()
DimPlot(pbmc, reduction = "tsne" , label = TRUE, pt.size = 0.5)# + NoLegend()


save.image("./scRNA_SLE_Jackson_SingleR.RData")
save.image("./scRNA_SLE_Jackson_SingleR_buckup.RData")

## check with canonial markers
canonial_markers <-  c("MS4A1", ## B
                       
                       "GNLY","NKG7", ## NK
                       
                       "CD8A", ## CD8+ T
                       "IL7R", "CCR7", ## Naive CD4+ T
                       #"S100A4", #"IL7R",## Memory CD4+ T
                       "PTPRC", ## general immune cell
                       "FCER1A","CST3", ## DCs
                       
                       "CD14","LYZ", ## CD14+ Mono
                       "FCGR3A","MS4A7", ## FCGR3A+ Mono
                       
                       
                       "PPBP" ## Platelet
                       
)

DotPlot(pbmc, features = canonial_markers) + RotatedAxis()



########################################################################################
##-------------------- Check Specific Genes in this Dataset --------------------------##
########################################################################################

Qgenes <- c("CD14",
            "CD34",
            "BST1",
            "CSTA",
            #"GUK1",
            "MEF2C",
            "FN1",
            "P3H1",
            "PHACTR4",
           "RGS12",
            "UBC"
)

pbmc$Group <- pbmc_ActiveIdent$group
table(pbmc$Group)

VlnPlot(pbmc, features = Qgenes, split.by = "Group", ncol = 5, pt.size = 0.1) 

DotPlot(pbmc, features = Qgenes, cols = c("blue","red"), split.by = "Group") 

DimPlot(pbmc, reduction = "umap" , label = TRUE, pt.size = 0.5, split.by = "Group")

FeaturePlot(pbmc, reduction = "tsne", features = Qgenes, ncol = 5 ,pt.size = 0.5) & NoLegend() & NoAxes()

DoHeatmap(subset(pbmc, downsample = 1000), features = Qgenes)

### scRNA 10 Immune complex boxplot -----------------------------------------------------------------

head(GEM)
dim(GEM) ## Genes * Cells
GEM[1:5,1:5]

GEM_10IC <- GEM[which(rownames(GEM) %in% Qgenes),]

GEM_10IC <- as.data.frame(GEM_10IC)
GEM_10IC[1:5, 1:5]
GEM_10IC$Gene <- rownames(GEM_10IC)

GEM_10IC_M <- melt(GEM_10IC)
head(GEM_10IC_M)
dim(GEM_10IC_M) ## 827480

## drop al the 0 cell
GEM_10IC_MC <- GEM_10IC_M[-which(GEM_10IC_M$value == 0),]
dim(GEM_10IC_MC) ## 122290
colnames(GEM_10IC_MC) <- c("Gene", "Cells" , "Expression Level")
head(GEM_10IC_MC)

head(pbmc_ActiveIdent)

table(pbmc_ActiveIdent$`pbmc@active.ident`)
pbmc_ActiveIdent_GEM <- left_join(GEM_10IC_MC, pbmc_ActiveIdent, by = c("Cells" = "Ident"))

head(pbmc_ActiveIdent_GEM)
str(pbmc_ActiveIdent_GEM)

colnames(pbmc_ActiveIdent_GEM)[4] <- "CellType"

my_comparision <- list(c("aHD", "aSLE"))

pbmc_ActiveIdent_GEM$group <- factor(pbmc_ActiveIdent_GEM$group, levels = c("aSLE" ,"aHD"))

ggplot(pbmc_ActiveIdent_GEM, aes( x = group, y = `Expression Level`)) + 
  #geom_jitter(aes(color = CellType),height = 0, size = 0.5, alpha = 0.5) + 
  geom_boxplot(aes(fill = group)) +
  facet_wrap(~ Gene, ncol =5) + theme_bw() + xlab("") +
  stat_compare_means(comparisons = my_comparision, paired = FALSE, method = "wilcox.test")

## read DEsingle result 

DEGs <- as.data.frame(read.csv2("Final/DEGs_Jackson_S_Genes_DESignle.csv"))

head(DEGs)

DEGs

DEGs_10IC <- DEGs[which(DEGs$X %in% Qgenes),]

DEGs_10IC$FC_SLEHC <- DEGs_10IC$norm_total_mean_2 /DEGs_10IC$norm_total_mean_1

head(DEGs_10IC)

write.csv2(DEGs_10IC, "DEGs_10_ImmuneComplex_0812.csv")

## Check the SLE Vs. HD expression level in cell type of each 10 Biomarkers

head(GEM_10IC_M)

head(pbmc_ActiveIdent)
table(pbmc_ActiveIdent$`pbmc@active.ident`)

pbmc_ActiveIdent_GEM_Raw <- left_join(GEM_10IC_M, pbmc_ActiveIdent, by = c("variable" = "Ident"))
colnames(pbmc_ActiveIdent_GEM_Raw) <- c("Gene" ,"Cells","ExpressLevel","CellType","SampleID","Sample_id","Group")

head(pbmc_ActiveIdent_GEM_Raw)

Gene_list <- unique(pbmc_ActiveIdent_GEM_Raw$Gene)
CellType_list <- unique(pbmc_ActiveIdent_GEM_Raw$CellType)

IC_CellType <- matrix(ncol = 6, nrow = length(Gene_list)*length(CellType_list))

z = 0
for(i in 1:length(Gene_list)){
  
  gene <- Gene_list[i]
  
  for(j in 1:length(CellType_list)){
    
    z <- z +1
    celltype <- CellType_list[j]
    
    sub <- pbmc_ActiveIdent_GEM_Raw[pbmc_ActiveIdent_GEM_Raw$Gene == gene & pbmc_ActiveIdent_GEM_Raw$CellType == celltype, ]
    
    SLE <- sub[sub$Group == "aSLE", 3]
    HC <- sub[sub$Group == "aHD", 3]
    
    SLE_mean <- mean(SLE)
    HC_mean <- mean(HC)
    Pvalue <- wilcox.test(SLE, HC)$p.value
    
    IC_CellType[z,1] <- as.character(gene)
    IC_CellType[z,2] <- as.character(celltype)
    IC_CellType[z,3] <- as.numeric(SLE_mean)
    IC_CellType[z,4] <- as.numeric(HC_mean)
    IC_CellType[z,5] <- as.numeric(Pvalue)
    IC_CellType[z,6] <- ifelse( SLE_mean > HC_mean, "Up", "Down")
      
  }
}

IC_CellType <- as.data.frame(IC_CellType)
colnames(IC_CellType) <- c("Gene", "CellType", "SLE_mean", "HC_mean", "P_value")

str(IC_CellType)
IC_CellType$P_value <- as.numeric(IC_CellType$P_value)

IC_CellType$stars <- cut(IC_CellType$P_value, breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                  label = c("****", "***", "**", "*", ""))

table(IC_CellType$stars)

IC_CellType[which(IC_CellType$stars == "****"),]
IC_CellType[which(IC_CellType$stars == "***"),]
IC_CellType[which(IC_CellType$stars == "**"),]
IC_CellType[which(IC_CellType$stars == "*"),]

write.csv2(IC_CellType, "Final/Ten_IC_CellType_ExpressionComparing.csv")

IC_CellType_sing <- IC_CellType[grep("\\*", IC_CellType$stars),]
CellType_Sign_Genes <- unique(IC_CellType_sing$Gene)

pbmc$Group <- factor(pbmc$Group, levels = c("aSLE" ,"aHD"))

VlnPlot(pbmc, features = CellType_Sign_Genes[1:3], split.by = "Group", ncol = 3, pt.size = 0.2) +
  theme(legend.position = "right")

VlnPlot(pbmc, features = CellType_Sign_Genes[4:6], split.by = "Group", ncol = 3, pt.size = 0.2)+
  theme(legend.position = "right")


## Quering the 2D gel Modular targets

VlnPlot(pbmc, features = c("ACTA2", "ADAR", "ZC3HAV1", "MOV10", "LAP3"), split.by = "Group", ncol = 5, pt.size = 0.2) + 
  theme(legend.position = "right")
  



### Extract the T Cell CD4+; T Cell CD8+; B Cell; Mono for find DEGs ---------------------------------------------

#library(BiocParallel)

pbmc$SingleR_Immune.labels
head(pbmc@active.ident)

table(Idents(pbmc))

T_CD4 <- subset(x = pbmc, idents = "T cells, CD4+")
table(T_CD4$Group) ## SLE:HC 15086:19897
T_CD4 <- as.SingleCellExperiment(T_CD4)
T_CD4_DEGs <- DEsingle(counts = T_CD4, group = factor(T_CD4$Group))
T_CD4_DEGs_classified <- DEtype(results = T_CD4_DEGs, threshold = 0.05)
T_CD4_DEGs_classified_sig <- T_CD4_DEGs_classified[T_CD4_DEGs_classified$pvalue.adj.FDR < 0.05, ]
dim(T_CD4_DEGs_classified_sig) ## 6011
write.csv2(T_CD4_DEGs_classified_sig, file = "DEGs_Jackson_S_Genes_DESignle_T_CD4.csv", quote = FALSE)



T_CD8 <- subset(x = pbmc, idents = "T cells, CD8+")
table(T_CD8$Group) ## SLE:HC 4471:4024
T_CD8 <- as.SingleCellExperiment(T_CD8)
T_CD8_DEGs <- DEsingle(counts = T_CD8, group = factor(T_CD8$Group))
T_CD8_DEGs_classified <- DEtype(results = T_CD8_DEGs, threshold = 0.05)
T_CD8_DEGs_classified_sig <- T_CD8_DEGs_classified[T_CD8_DEGs_classified$pvalue.adj.FDR < 0.05, ]
dim(T_CD8_DEGs_classified_sig) ## 2967
write.csv2(T_CD8_DEGs_classified_sig, file = "DEGs_Jackson_S_Genes_DESignle_T_CD8.csv", quote = FALSE)


B_Cell <- subset(x = pbmc, idents = "B cells")
table(B_Cell$SingleR_Immune.labels)
table(B_Cell$Group) ## SLE:HC 1124:2418
B_Cell <- as.SingleCellExperiment(B_Cell)
B_Cell_DEGs <- DEsingle(counts = B_Cell, group = factor(B_Cell$Group))
B_Cell_DEGs_classified <- DEtype(results = B_Cell_DEGs, threshold = 0.05)
B_Cell_DEGs_classified_sig <- B_Cell_DEGs_classified[B_Cell_DEGs_classified$pvalue.adj.FDR < 0.05, ]
dim(B_Cell_DEGs_classified_sig) ## 1223
write.csv2(B_Cell_DEGs_classified_sig, file = "DEGs_Jackson_S_Genes_DESignle_BCell.csv", quote = FALSE)



pbmc$SingleR_Immune.labels.Mono <- Idents(pbmc)
table(pbmc$SingleR_Immune.labels.Mono)
pbmc$SingleR_Immune.labels.Mono <- sub( "Monocytes, CD14+", "Monocytes", pbmc$SingleR_Immune.labels.Mono)
pbmc$SingleR_Immune.labels.Mono <- sub( "Monocytes, CD16+", "Monocytes", pbmc$SingleR_Immune.labels.Mono)
table(pbmc$SingleR_Immune.labels.Mono)

Idents(object = pbmc) <- pbmc$SingleR_Immune.labels.Mono

table(Idents(pbmc))
Mono_Cell <- subset(x = pbmc, idents = "Monocytes+")
table(Mono_Cell$SingleR_Immune.labels.Mono)
table(Mono_Cell$Group) ## SLE:HC 12108:7150
Mono_Cell <- as.SingleCellExperiment(Mono_Cell)
Mono_Cell_DEGs <- DEsingle(counts = Mono_Cell, group = factor(Mono_Cell$Group))
Mono_Cell_DEGs_classified <- DEtype(results = Mono_Cell_DEGs, threshold = 0.05)
Mono_Cell_DEGs_classified_sig <- Mono_Cell_DEGs_classified[Mono_Cell_DEGs_classified$pvalue.adj.FDR < 0.05, ]
dim(Mono_Cell_DEGs_classified_sig) ## 2481
head(Mono_Cell_DEGs_classified_sig)
write.csv2(Mono_Cell_DEGs_classified_sig, file = "DEGs_Jackson_S_Genes_DESignle_MonoCell.csv", quote = FALSE)






########################################################################################################################
########################## ---------- For BPMA Array -------------------------#####################################

BPMA <- c("LGALS9","TNFRSF1B","IGFBP2","PF4","BST1","CD40","TNFRSF13C","ALCAM","TFPI","CRP","VCAM1","SPP1","CD14","CSTA","VSIG4","CD177")
length(BPMA) ## 16

#pbmc$Group <- pbmc_ActiveIdent$group
#table(pbmc$Group)

VlnPlot(pbmc, features = BPMA, split.by = "Group", pt.size = 1) 
## CRP & CD177 didn't find
VlnPlot(pbmc, features = BPMA, pt.size = 1)

DotPlot(pbmc, features = BPMA, cols = c("blue","red"), split.by = "Group") 
DotPlot(pbmc, features = BPMA) 

DimPlot(pbmc, reduction = "umap" , label = TRUE, pt.size = 0.5, split.by = "Group")

FeaturePlot(pbmc, reduction = "tsne", features = BPMA, ncol = 4 ,pt.size = 0.5) & NoLegend() & NoAxes()

DoHeatmap(subset(pbmc, downsample = 5000), features = BPMA)

DoHeatmap(pbmc, features = BPMA)


### scRNA expression of 16 BPMA list -----------------------------------------------------------------

head(GEM)
dim(GEM) ## Genes * Cells
GEM[1:5,1:5]

GEM_BPMA <- GEM[which(rownames(GEM) %in% BPMA),]

GEM_BPMA <- as.data.frame(GEM_BPMA)
GEM_BPMA$Gene <- rownames(GEM_BPMA)

GEM_BPMA_M <- melt(GEM_BPMA)
head(GEM_BPMA_M)
dim(GEM_BPMA_M) ## 82748


## Check the SLE Vs. HD expression level in cell type of each 10 Biomarkers

head(GEM_BPMA_M)
#GEM_BPMA_M$variable <- GEM_BPMA_M$Gene

head(pbmc_ActiveIdent)
table(pbmc_ActiveIdent$`pbmc@active.ident`)

pbmc_ActiveIdent_GEM_Raw <- left_join(GEM_BPMA_M, pbmc_ActiveIdent, by = c("variable" = "Ident"))
head(pbmc_ActiveIdent_GEM_Raw)
colnames(pbmc_ActiveIdent_GEM_Raw) <- c("Gene" ,"Cells","ExpressLevel","CellType","SampleID","Sample_id","Group")

head(pbmc_ActiveIdent_GEM_Raw)
dim(pbmc_ActiveIdent_GEM_Raw)

write.csv(pbmc_ActiveIdent_GEM_Raw,"BPMA_Jackson_scRNA_serum_expression.csv")

Gene_list <- unique(pbmc_ActiveIdent_GEM_Raw$Gene)
CellType_list <- unique(pbmc_ActiveIdent_GEM_Raw$CellType)

IC_CellType <- matrix(ncol = 6, nrow = length(Gene_list)*length(CellType_list))

z = 0
for(i in 1:length(Gene_list)){
  
  gene <- Gene_list[i]
  
  for(j in 1:length(CellType_list)){
    
    z <- z +1
    celltype <- CellType_list[j]
    
    sub <- pbmc_ActiveIdent_GEM_Raw[pbmc_ActiveIdent_GEM_Raw$Gene == gene & pbmc_ActiveIdent_GEM_Raw$CellType == celltype, ]
    
    SLE <- sub[sub$Group == "aSLE", 3]
    HC <- sub[sub$Group == "aHD", 3]
    
    SLE_mean <- mean(SLE)
    HC_mean <- mean(HC)
    Pvalue <- wilcox.test(SLE, HC)$p.value
    
    IC_CellType[z,1] <- as.character(gene)
    IC_CellType[z,2] <- as.character(celltype)
    IC_CellType[z,3] <- as.numeric(SLE_mean)
    IC_CellType[z,4] <- as.numeric(HC_mean)
    IC_CellType[z,5] <- as.numeric(Pvalue)
    IC_CellType[z,6] <- ifelse( SLE_mean > HC_mean, "Up", "Down")
    
  }
}
head(IC_CellType)
IC_CellType <- as.data.frame(IC_CellType)
colnames(IC_CellType) <- c("Gene", "CellType", "SLE_mean", "HC_mean", "P_value","UP-Down")

str(IC_CellType)
IC_CellType$P_value <- as.numeric(IC_CellType$P_value)

IC_CellType$stars <- cut(IC_CellType$P_value, breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                         label = c("****", "***", "**", "*", ""))

table(IC_CellType$stars)

IC_CellType[which(IC_CellType$stars == "****"),]
IC_CellType[which(IC_CellType$stars == "***"),]
IC_CellType[which(IC_CellType$stars == "**"),]
IC_CellType[which(IC_CellType$stars == "*"),]

write.csv2(IC_CellType, "Final/BPMA_Jackson_CellType_ExpressionComparing.csv")


