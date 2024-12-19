#Install immunedeconv from CRAN
install.packages("remotes")
remotes::install_github("omnideconv/immunedeconv")

library(SeuratDisk)
library(SeuratData)
library(sctransform)
library(ggplot2)
library(magrittr)
library(edgeR)
library(stringr)
library(scales)
library(dplyr)
library(readr)
library(Matrix)
library(gridExtra)
library(Seurat)
library(ggpubr)
library(RColorBrewer)
library(readxl)
library(patchwork)
library(harmony)
library(purrr)
library(tidyverse)
library(openxlsx)
library(immunedeconv)
library(sctransform)
library(RcppArmadillo)
library(harmony)
library(cowplot) 
library(HGNChelper)
library(viridis)
library(scCustomize)
library(qs)
library(Matrix)

#Make three dataset seurat object
#1--NatureCancer=====
GBM.data <- read_csv("~/data/Richards_NatureCancer_GBM_scRNAseq_counts.csv") %>% as.data.frame()
row.names(GBM.data) <- GBM.data[,1] %>% as.vector()
GBM.data$'...1' <- NULL

seurat_ob <- CreateSeuratObject(counts = GBM.data, project = "gbm", min.cells = 3, min.features = 200)
seurat_ob[["percent.mt"]] <- PercentageFeatureSet(seurat_ob, pattern = "^MT-")
seurat_ob <- subset(seurat_ob, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 30)

VlnPlot(seurat_ob, features = c("nFeature_RNA", 
                                "nCount_RNA", 
                                "percent.mt"), ncol = 3)

#2--GSE131928=====
cells1_proc <- read_csv("~data/GSE131928/cells1.proc.tsv", col_names = FALSE) %>% as.data.frame()
cells2_proc <- read_csv("~data/GSE131928/cells2.proc.tsv", col_names = FALSE) %>% as.data.frame()
genes1 <- read_delim("~data/GSE131928/genes1.tsv", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE) %>% as.data.frame()
genes2 <- read_delim("~data/GSE131928/genes2.tsv", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE) %>% as.data.frame()
count1 <- Matrix::readMM("~data/GSE131928/IDHwtGBM.processed.10X.counts.mtx") %>% as.matrix() %>% as.data.frame()
count2 <- Matrix::readMM("~data/GSE131928/IDHwtGBM.processed.10X.counts.2.mtx") %>% as.matrix() %>% as.data.frame()
merge1 <- cbind(genes1,count1)
merge1 <- merge1[,-1]
merge1_clean <- merge1 %>% group_by(X2) %>% summarise(across(everything(), sum))
merge2 <- cbind(genes2,count2)
merge2 <- merge2[,-1]
merge2_clean <- merge2 %>% group_by(X2) %>% summarise(across(everything(), sum))
GBM.data <- merge(merge1_clean, merge2_clean, by.x = "X2", by.y = "X2", all = TRUE)
GBM.data[is.na(GBM.data)] <- 0
#remove 105
row.names(GBM.data) <- GBM.data$X2
GBM.data$X2 <- NULL
colnames(GBM.data) <- append(cells1_proc$X1, cells2_proc$X1)
GBM.data_8sample <- GBM.data[,-grep(colnames(GBM.data),pattern = "105_")]



seurat_ob <- CreateSeuratObject(counts = GBM.data_8sample, project = "gbm", min.cells = 3, min.features = 200)
seurat_ob[["percent.mt"]] <- PercentageFeatureSet(seurat_ob, pattern = "^MT-")
seurat_ob <- subset(seurat_ob, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 15)
VlnPlot(seurat_ob, features = c("nFeature_RNA", 
                                "nCount_RNA", 
                                "percent.mt"), ncol = 3)


#3--GSE182109 without "LGG-04" ,"LGG-03", "ndGBM-11","ndGBM-10"=====
count <- Matrix::readMM("~data/GSE182109/matrix.mtx") %>% as.matrix() %>% as.data.frame()
barcodes <- read_csv("~data/GSE182109/barcodes.tsv", 
                     col_names = FALSE) %>% as.data.frame()
genes <- read_delim("~data/GSE182109/genes.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    col_names = FALSE, trim_ws = TRUE) %>% as.data.frame()

GBM.data <- cbind(genes,count)
row.names(GBM.data) <-  GBM.data$X1
GBM.data <- GBM.data[,-c(1,2)]
colnames(GBM.data) <- barcodes$X1

seurat_ob <- CreateSeuratObject(counts = GBM.data, project = "gbm", min.cells = 3, min.features = 200)
seurat_ob <- PercentageFeatureSet(seurat_ob, pattern = "^MT-", col.name = "percent.mt")
seurat_ob <- subset(seurat_ob, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 20)
seurat_ob <- subset(x = seurat_ob, subset = Patient %in% c("ndGBM-01","ndGBM-02","ndGBM-03","ndGBM-04","ndGBM-05",
                                                           "ndGBM-06","ndGBM-07","ndGBM-08","ndGBM-09",
                                                           "rGBM-01","rGBM-02","rGBM-03", "rGBM-04","rGBM-05"))

VlnPlot(seurat_ob, features = c("nFeature_RNA", 
                                "nCount_RNA", 
                                "percent.mt"), ncol = 3)

Final_Metadata <- read_csv("~data/Meta_GBM.txt") %>% as.data.frame()
Final_Metadata <- Final_Metadata[-1,] %>% as.data.frame()
Patient <- factor(Final_Metadata$Patient)
names(Patient) <- Final_Metadata$NAME
seurat_ob <- AddMetaData(
  object = seurat_ob,
  metadata = list(Patient),
  col.name = c("Patient")
)

#Runnning sctransform on a Seurat object------------------------------------------------------------------------------------------------------
library(future)
plan("sequential", workers = 24) 
options(future.globals.maxSize= 20000 * 1024^2)

seurat_ob_SCT <- SCTransform(seurat_ob,
                             vars.to.regress = "percent.mt",
                             ncells = 10000,
                             return.only.var.genes = FALSE,
                             variable.features.n = 5000,
                             seed.use = 42)
seurat_ob_SCT <- RunPCA(seurat_ob_SCT,npcs = 50) 
seurat_ob_SCT@meta.data$orig.ident <- as.character(seurat_ob_SCT@meta.data$orig.ident)
seurat_ob_SCT@meta.data$orig.ident <- as.factor(seurat_ob_SCT@meta.data$orig.ident)

seurat_ob_SCT_hmy <- RunHarmony(seurat_ob_SCT,
                                assay.use = "SCT",
                                Reductions.use = "pca",
                                dims.use = 1:50,
                                group.by.vars = "orig.ident")
seurat_ob_SCT_hmy <- RunUMAP(seurat_ob_SCT_hmy,
                             assay = "SCT",
                             reduction = "harmony",
                             dims = 1:50,
                             n.neighbors = 30, min.dist = 0.15, seed.use = 42) #GSE182109 n.neighbors = 20


DimPlot(seurat_ob_SCT_hmy, reduction = "umap", label = TRUE)