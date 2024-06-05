library(devtools)
install_github("immunogenomics/harmony")
#restart
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("granulator")
#restart

#-----------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(magrittr)
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
library(tidyverse)
library(sctransform)

count <- Matrix::readMM("data2/GSE182109/matrix.mtx") %>% as.matrix() %>% as.data.frame()
barcodes <- read_csv("data2/GSE182109/barcodes.tsv", 
                     col_names = FALSE) %>% as.data.frame()
genes <- read_delim("data2/GSE182109/genes.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    col_names = FALSE, trim_ws = TRUE) %>% as.data.frame()
GBM.data <- cbind(genes,count)
row.names(GBM.data) <-  GBM.data$X1
GBM.data <- GBM.data[,-c(1,2)]
colnames(GBM.data) <- barcodes$X1

Final_Metadata <- read_csv("data/Meta_GBM.txt") %>% as.data.frame()
Final_Metadata <- Final_Metadata[-1,] %>% as.data.frame()

seurat_ob <- CreateSeuratObject(counts = GBM.data, project = "gbm", min.cells = 3, min.features = 200)
seurat_ob <- PercentageFeatureSet(seurat_ob, pattern = "^MT-", col.name = "percent.mt")

target <- colnames(seurat_ob)
Final_Metadata <- Final_Metadata[match(target, Final_Metadata$NAME),] %>% as.data.frame()
Patient <- factor(Final_Metadata$Patient)
names(Patient) <- colnames(seurat_ob)
Type <- factor(Final_Metadata$Type)
names(Type) <- colnames(seurat_ob)
Assignment <- factor(Final_Metadata$Assignment)
names(Assignment) <- colnames(seurat_ob)
seurat_ob <- AddMetaData(
  object = seurat_ob,
  metadata = list(Patient,Type,Assignment),
  col.name = c("Patient","Type","Assignment")
)
seurat_ob <- subset(seurat_ob, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 20)
remove(barcodes,count,genes,Assignment,Patient,Type)

seurat_ob <- subset(x = seurat_ob, subset = Patient %in% c("ndGBM-01","ndGBM-02","ndGBM-03","ndGBM-04","ndGBM-05","ndGBM-06","ndGBM-07","ndGBM-08","ndGBM-09",
                                                           "rGBM-01","rGBM-02","rGBM-03","rGBM-04","rGBM-05"))
VlnPlot(seurat_ob, features = c("nFeature_RNA"), ncol = 1)
#--------------------------------------------------------------------------------------------------------------------------------------
seurat_ob_SCT <- SCTransform(seurat_ob,
                             vars.to.regress = "percent.mt",
                             ncells = 10000,
                             variable.features.n = 5000,
                             return.only.var.genes = FALSE,
                             seed.use = 42)
seurat_ob_SCT <- RunPCA(seurat_ob_SCT,npcs = 50) #npcs = 50
seurat_ob_SCT@meta.data$orig.ident <- as.character(seurat_ob_SCT@meta.data$orig.ident)
seurat_ob_SCT@meta.data$orig.ident <- as.factor(seurat_ob_SCT@meta.data$orig.ident)
seurat_ob_SCT_hmy <- RunHarmony(seurat_ob_SCT,
                                assay.use = "SCT",
                                reduction = "pca",
                                dims.use = 1:50,
                                group.by.vars = "orig.ident")
seurat_ob_SCT_hmy <- RunUMAP(seurat_ob_SCT_hmy,
                             assay = "SCT",
                             reduction = "harmony",
                             dims = 1:50,
                             n.neighbors = 20, min.dist = 0.15, seed.use = 42)
seurat_ob_SCT_hmy <- FindNeighbors(seurat_ob_SCT_hmy, reduction = "harmony", dims = 1:50)
seurat_ob_SCT_hmy <- FindClusters(seurat_ob_SCT_hmy, resolution = 1.0)
DimPlot(seurat_ob_SCT_hmy, label = TRUE, repel = TRUE, group.by = "seurat_clusters", raster=FALSE)

GSE182109_metadata <- read_csv("data/GSE182109_ndGBM_rGBM_Nond1011_metadata.txt") %>% as.data.frame()
target <- colnames(seurat_ob_SCT_hmy)
GSE182109_metadata <- GSE182109_metadata[match(target, GSE182109_metadata$cells),] %>% as.data.frame()
CellType1 <- factor(GSE182109_metadata$CellType1)
names(CellType1) <- colnames(seurat_ob_SCT_hmy)
seurat_ob_SCT_hmy <- AddMetaData(
  object = seurat_ob_SCT_hmy,
  metadata = list(CellType1),
  col.name = c("CellType1")
)
remove(target,CellType1)
DimPlot(seurat_ob_SCT_hmy, label = TRUE, repel = TRUE, group.by = "CellType1", raster=FALSE)

#scType----------------------------------------------------------------------------------------------------------------
install.packages("HGNChelper")
library(openxlsx)
library(HGNChelper)
library(dplyr)

scData <- seurat_Tumor
#Run gene_sets_prepare.R, sctype_score_.R
db_ <- "data/ScTypeDB_GBM_integrate.xlsx"
tissue <- "GBM" #Immune system, GBM, GAM, Brain
gs_list <- gene_sets_prepare(db_, tissue)
es.max <- sctype_score(scRNAseqData = scData[["SCT"]]@scale.data,
                       scaled = TRUE, 
                       gs = gs_list$gs_positive,
                       gs2 = gs_list$gs_negative)
cL_resutls <- do.call("rbind", lapply(unique(scData@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(scData@meta.data[scData@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scData@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"

scData@meta.data$Brain <- ""
for(j in unique(sctype_scores$cluster)){
  cl_type <- sctype_scores[sctype_scores$cluster==j,] 
  scData@meta.data$Brain[scData@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
}
DimPlot(scData, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Brain', raster=FALSE)

scData@meta.data$GBM <- ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  scData@meta.data$GBM[scData@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
DimPlot(scData, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'GBM', raster=FALSE)

scData@meta.data$Immune <- ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  scData@meta.data$Immune[scData@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
DimPlot(scData, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Immune', raster=FALSE)

seurat_immune <- scData
remove(scData,cL_resutls,cl_type,es.max,gs_list,sctype_scores,j,tissue)
#--------------------------------------------------------------------
seurat_ob_SCT_hmy <- scData

target_cluster <- seurat_ob_SCT_hmy@meta.data$SCT_snn_res.1 %>% as.character()

T_cells <- c("6","8","23","24","35","42")
B_cells <- c("29")
Endothelial_cells <- c("30")
Mural_cells <- c("33")
Oligodendrocytes <- c("20")
Tumor_cells <- c("2","4","7","10","11","13","16","17","22","25","26","28")
Proliferating_tumor_cells <- c("9","14","40")
GAMs <- c("0","1","3","5","12","15","18","19","21","31","32","34","36","37","41")
Others <- c("27","38","39")

target_cluster[which(target_cluster %in% T_cells)] <- "T cells"
target_cluster[which(target_cluster %in% B_cells)] <- "B cells"
target_cluster[which(target_cluster %in% Endothelial_cells)] <- "Endothelial cells"
target_cluster[which(target_cluster %in% Mural_cells)] <- "Mural cells"
target_cluster[which(target_cluster %in% Oligodendrocytes)] <- "Oligodendrocytes"
target_cluster[which(target_cluster %in% Tumor_cells)] <- "Tumor cells"
target_cluster[which(target_cluster %in% GAMs)] <- "GAMs"
target_cluster[which(target_cluster %in% Proliferating_tumor_cells)] <- "Proliferating tumor cells"
target_cluster[which(target_cluster %in% Others)] <- "Others"

SCT_res.1_board <- factor(target_cluster)
names(SCT_res.1_board) <- colnames(seurat_ob_SCT_hmy)
seurat_ob_SCT_hmy <- AddMetaData(
  object = seurat_ob_SCT_hmy,
  metadata = list(SCT_res.1_board),
  col.name = c("SCT_res.1_board")
)
DimPlot(seurat_ob_SCT_hmy, label = TRUE, repel = TRUE, group.by = "SCT_res.1_board", raster=FALSE)

#BT cells-----------------------------------------------------------------------------
seurat_immune <- subset(x = seurat_ob_SCT_hmy, subset = SCT_res.1_board %in% c("B cells","T cells"))
seurat_immune <- FindNeighbors(seurat_immune, reduction = "harmony", dims = 1:50)
seurat_immune <- FindClusters(seurat_immune, resolution = 1.0)
DimPlot(seurat_immune, reduction = "umap", label = TRUE)

seurat_immune <- scData

BT_immune <- seurat_immune@meta.data$Immune %>% as.character()
BT_immune[grep("T cells",BT_immune)] <- "T cells"
BT_immune[grep("NKT-like",BT_immune)] <- "NKT-like cells"

BT_immune <- factor(BT_immune)
names(BT_immune) <- colnames(seurat_immune)
seurat_immune <- AddMetaData(
  object = seurat_immune,
  metadata = list(BT_immune),
  col.name = c("BT_immune")
)
DimPlot(seurat_immune, label = TRUE, repel = TRUE, group.by = "BT_immune", raster=FALSE)
#GAM-----------------------------------------------------------------------------
seurat_GAM <- subset(x = seurat_ob_SCT_hmy, subset = SCT_res.1_board %in% c("GAMs"))
seurat_GAM <- FindNeighbors(seurat_GAM, reduction = "harmony", dims = 1:50)
seurat_GAM <- FindClusters(seurat_GAM, resolution = 1.0)#0.4 general
DimPlot(seurat_GAM, reduction = "umap", label = TRUE)

seurat_GAM <- scData

target_cluster <- rep("GAMs",length(colnames(seurat_GAM)))
target_cluster[which((seurat_GAM@meta.data$SCT_snn_res.1 %>% as.character()) %in% c("12","25"))] <- "Tumor cells"
target_cluster[which((seurat_GAM@meta.data$SCT_snn_res.1 %>% as.character()) %in% c("18"))] <- "Dendritic cells"
GAM_temp <- factor(target_cluster)
names(GAM_temp) <- colnames(seurat_GAM)
seurat_GAM <- AddMetaData(
  object = seurat_GAM,
  metadata = list(GAM_temp),
  col.name = c("GAM_temp")
)
DimPlot(seurat_GAM, label = TRUE, repel = TRUE, group.by = "GAM_temp", raster=FALSE)

seurat_GAM_2 <- subset(x = seurat_GAM, subset = SCT_snn_res.1 %in% c("6"))
seurat_GAM_2 <- FindNeighbors(seurat_GAM_2, reduction = "harmony", dims = 1:50)
seurat_GAM_2 <- FindClusters(seurat_GAM_2, resolution = 1.0)
DimPlot(seurat_GAM_2, reduction = "umap", label = TRUE)
FeaturePlot(seurat_GAM_2, features = c("TMEM119","CX3CR1"))

target_cluster <- seurat_GAM_2@meta.data$SCT_snn_res.1 %>% as.character()
target_cluster[which(target_cluster %in% c("0","5","7","9","2"))] <- "Macrophage-like GAMs"
target_cluster[which(target_cluster %in% c("1","3","4","6","8"))] <- "Microglia-like GAMs"
GAM_group <- factor(target_cluster)
names(GAM_group) <- colnames(seurat_GAM_2)
seurat_GAM_2 <- AddMetaData(
  object = seurat_GAM_2,
  metadata = list(GAM_group),
  col.name = c("GAM_group")
)
DimPlot(seurat_GAM_2, label = TRUE, repel = TRUE, group.by = "GAM_group", raster=FALSE)

GAM_meta <- data.frame(cell_ID = colnames(seurat_GAM),
                       GAM_group = seurat_GAM@meta.data$GAM_temp)
temp <- data.frame(cell_ID = colnames(seurat_GAM_2),
                   GAM_group_2 = seurat_GAM_2@meta.data$GAM_group)
GAM_meta <- merge(x=GAM_meta, y=temp, by.x = "cell_ID", by.y = "cell_ID", all = TRUE) %>% as.data.frame()
GAM_meta$GAM_group <- as.character(GAM_meta$GAM_group)
GAM_meta$GAM_group_2 <- as.character(GAM_meta$GAM_group_2)
GAM_meta$GAM_group[which(is.na(GAM_meta$GAM_group_2)==FALSE)] <- GAM_meta$GAM_group_2[which(is.na(GAM_meta$GAM_group_2)==FALSE)]
GAM_meta$GAM_group <- as.factor(GAM_meta$GAM_group)
GAM_meta$GAM_group_2 <- as.factor(GAM_meta$GAM_group_2)
GAM_group_2 <- factor(GAM_meta$GAM_group )
names(GAM_group_2) <- colnames(seurat_GAM)
seurat_GAM <- AddMetaData(
  object = seurat_GAM,
  metadata = list(GAM_group_2),
  col.name = c("GAM_group_2")
)
DimPlot(seurat_GAM, label = TRUE, repel = TRUE, group.by = "GAM_group_2", raster=FALSE)

seurat_GAM_3 <- subset(x = seurat_GAM, subset = GAM_group_2 %in% c("GAMs"))
seurat_GAM_3 <- FindNeighbors(seurat_GAM_3, reduction = "harmony", dims = 1:50)
seurat_GAM_3 <- FindClusters(seurat_GAM_3, resolution = 1.0)
DimPlot(seurat_GAM_3, reduction = "umap", label = TRUE)
FeaturePlot(seurat_GAM_3, features = c("TMEM119","CX3CR1","TGFBI"))

target_cluster <- rep("Microglia-like GAMs",length(colnames(seurat_GAM_3)))
target_cluster[which((seurat_GAM_3@meta.data$SCT_snn_res.1 %>% as.character()) %in% c("12","14","5","7","9","2","15","20","10","16","22","21","14","11","13","8","17"))] <- "Macrophage-like GAMs"

GAM_temp_2 <- factor(target_cluster)
names(GAM_temp_2) <- colnames(seurat_GAM_3)
seurat_GAM_3 <- AddMetaData(
  object = seurat_GAM_3,
  metadata = list(GAM_temp_2),
  col.name = c("GAM_temp_2")
)
DimPlot(seurat_GAM_3, label = TRUE, repel = TRUE, group.by = "GAM_temp_2", raster=FALSE)

temp <- data.frame(cell_ID = colnames(seurat_GAM_3),
                   GAM_group_3 = seurat_GAM_3@meta.data$GAM_temp_2)
GAM_meta <- merge(x=GAM_meta, y=temp, by.x = "cell_ID", by.y = "cell_ID", all = TRUE) %>% as.data.frame()
GAM_meta$GAM_group <- as.character(GAM_meta$GAM_group)
GAM_meta$GAM_group_3 <- as.character(GAM_meta$GAM_group_3)
GAM_meta$GAM_group[which(is.na(GAM_meta$GAM_group_3)==FALSE)] <- GAM_meta$GAM_group_3[which(is.na(GAM_meta$GAM_group_3)==FALSE)]
GAM_meta$GAM_group <- as.factor(GAM_meta$GAM_group)
GAM_meta$GAM_group_3 <- as.factor(GAM_meta$GAM_group_3)
GAM_group_2 <- factor(GAM_meta$GAM_group )
names(GAM_group_2) <- colnames(seurat_GAM)
GAM_meta$GAM_group[grep("GAMs",GAM_meta$GAM_group)] <- GAM_meta$GAM_group_3[grep("GAMs",GAM_meta$GAM_group)]
GAM_meta$GAM_group <- as.factor(GAM_meta$GAM_group)
GAM_meta$GAM_group_2 <- as.factor(GAM_meta$GAM_group_2)

GAM_group <- factor(GAM_meta$GAM_group)
names(GAM_group) <- colnames(seurat_GAM)
seurat_GAM <- AddMetaData(
  object = seurat_GAM,
  metadata = list(GAM_group),
  col.name = c("GAM_group")
)
DimPlot(seurat_GAM, label = TRUE, repel = TRUE, group.by = "GAM_group", raster=FALSE)

#Tumor cells-------------------------------------------------------------------------------------------------------------------------------
seurat_Tumor <- subset(x = seurat_ob_SCT_hmy, subset = SCT_res.1_board %in% c("Tumor cells","Proliferating tumor/GAM cells","Mural cells"))
seurat_Tumor <- FindNeighbors(seurat_Tumor, reduction = "harmony", dims = 1:50)
seurat_Tumor <- FindClusters(seurat_Tumor, resolution = 1.0)
DimPlot(seurat_Tumor, reduction = "umap", label = TRUE)

Const_meta <- data.frame(cells = colnames(seurat_ob_SCT_hmy),
                         Patient =  seurat_ob_SCT_hmy@meta.data$Patient,
                         Assignment = seurat_ob_SCT_hmy@meta.data$Assignment,
                         Type = seurat_ob_SCT_hmy@meta.data$Type,
                         Broad_celltype = seurat_ob_SCT_hmy@meta.data$SCT_res.1_board)
BT_immune <- data.frame(cell_ID = colnames(seurat_immune),
                        BT_group = seurat_immune@meta.data$BT_immune)
GAM_meta <- data.frame(cell_ID = colnames(seurat_GAM),
                       GAM_group = seurat_GAM@meta.data$GAM_group)

Const_meta <- merge(x=Const_meta, y=BT_immune, by.x = "cells", by.y = "cell_ID", all = TRUE)
Const_meta <- merge(x=Const_meta, y=GAM_meta, by.x = "cells", by.y = "cell_ID", all = TRUE)
target <- colnames(seurat_ob_SCT_hmy)
Const_meta <- Const_meta[match(target, Const_meta$cells),] %>% as.data.frame()
Const_meta$CellType1 <- Const_meta$Broad_celltype

Const_meta$CellType1 <- as.character(Const_meta$CellType1)
Const_meta$BT_group <- as.character(Const_meta$BT_group)
Const_meta$GAM_group <- as.character(Const_meta$GAM_group)
Const_meta$CellType1[which(is.na(Const_meta$BT_group)==FALSE)] <- Const_meta$BT_group[which(is.na(Const_meta$BT_group)==FALSE)]
Const_meta$CellType1[which(is.na(Const_meta$GAM_group)==FALSE)] <- Const_meta$GAM_group[which(is.na(Const_meta$GAM_group)==FALSE)]
Const_meta$CellType1 <- as.factor(Const_meta$CellType1)
Const_meta$BT_group <- as.factor(Const_meta$BT_group)
Const_meta$GAM_group <- as.factor(Const_meta$GAM_group)

CellType1 <- factor(Const_meta$CellType1)
names(CellType1) <- colnames(seurat_ob_SCT_hmy)
seurat_ob_SCT_hmy <- AddMetaData(
  object = seurat_ob_SCT_hmy,
  metadata = list(CellType1),
  col.name = c("CellType1")
)
DimPlot(seurat_ob_SCT_hmy, label = TRUE, repel = TRUE, group.by = "CellType1", raster=FALSE)

write.csv(Const_meta, file = "data/GSE131928_SCT_Celltype11_metadata.txt")

#Findmarker---------------------------------------------------------------------------------------------------
CellType<- c("Proliferating tumor cells","Microglia-like GAMs","Macrophage-like GAMs","NKT-like cells","Tumor cells","T cells","B cells","Oligodendrocytes",
             "Dendritic cells","Endothelial cells","Mural cells")
seurat_ob_SCT_hmy <- subset(x = seurat_ob_SCT_hmy, subset = CellType1 %in% CellType)
Idents(seurat_ob_SCT_hmy) <- seurat_ob_SCT_hmy@meta.data$CellType1
seurat_ob_SCT_hmy <- PrepSCTFindMarkers(seurat_ob_SCT_hmy)
cluster.markers <- FindAllMarkers(seurat_ob_SCT_hmy,
                                  assay = "SCT",
                                  only.pos = TRUE,
                                  test.use = "wilcox",
                                  recorrect_umi = FALSE,
                                  min.pct = 0.6, #0.2
                                  logfc.threshold = 0.1, #default 0.25
                                  min.diff.pct = 0.1, #0.1
                                  return.thresh = 0.05)
cluster.markers <- cluster.markers[-which(cluster.markers$cluster=="Others"),] %>% as.data.frame()
ex_matrix <- seurat_ob_SCT_hmy[["SCT"]]@data %>% as.data.frame()
target_celltype <- cluster.markers$cluster %>% unique()
for (i in 1:length(target_celltype)) {
  target_gene <- ex_matrix[cluster.markers$gene[which(cluster.markers$cluster==target_celltype[i])],] %>% t() %>% as.data.frame()
  target_gene$celltype <- seurat_ob_SCT_hmy@meta.data$CellType1
  target_gene$celltype[which(target_gene$celltype!=target_celltype[i])] <- "All_others"
  gene_mean_ex <- target_gene %>% group_by(celltype) %>% summarise_all(mean) %>% t() %>%  as.data.frame()
  if(i==1){
    gene_mean_ex_total <- gene_mean_ex[-1,]
  }else{
    gene_mean_ex_total <- rbind(gene_mean_ex_total,gene_mean_ex[-1,])
  }
}
colnames(gene_mean_ex_total) <- c("target_expression_mean","others_expression_mean")
cluster.markers <- cbind(cluster.markers,gene_mean_ex_total) %>% as.data.frame()

write.csv(cluster.markers, file = "~/Findallmarker_GSE182109_11celltype_211_DatasetWhole_MinDist.csv")
remove(cluster.markers)

df <- cluster.markers %>% as.data.frame()
gene_n <- 20
df$pct.diff <- df$pct.1 - df$pct.2
target_celltype <- cluster.markers$cluster %>% unique()
pval0_num_2 <- c()
pval0_num <- c()
celltype_name <- c()
for (i in 1:length(target_celltype)) {
  total_gene_num <- length(df$gene[which(df$cluster==target_celltype[i])]) %>% as.integer()
  pval0_gene_num <- length(which(df$p_val_adj[which(df$cluster==target_celltype[i])]==0))
  pval0_num <- append(pval0_num,pval0_gene_num)
  pval0_num_2 <- append(pval0_num_2,total_gene_num)
  celltype_name <- append(celltype_name,target_celltype[i])
}
sum_table <- data.frame(celltype = celltype_name,
                        gene_num = pval0_num_2,
                        pval0 = pval0_num)

df <- df[which(df$p_val_adj==0),] %>% as.data.frame()
df <- df %>% group_by(cluster) %>%
  top_n(gene_n*4, target_expression_mean) %>%  #(p_val_adj很多都為0)
  top_n(gene_n*2, pct.diff) %>%
  top_n(gene_n, avg_log2FC)

str_c(df$gene[which(df$cluster=="Macrophage-like GAMs")] %>% unlist() %>% as.vector(),collapse = ",")
str_c(df$gene[which(df$cluster=="Microglia-like GAMs")] %>% unlist() %>% as.vector(),collapse = ",")
str_c(df$gene[which(df$cluster=="Dendritic cells")] %>% unlist() %>% as.vector(),collapse = ",")
str_c(df$gene[which(df$cluster=="Mural cells")] %>% unlist() %>% as.vector(),collapse = ",")
str_c(df$gene[which(df$cluster=="Endothelial cells")] %>% unlist() %>% as.vector(),collapse = ",")
str_c(df$gene[which(df$cluster=="Proliferating tumor cells")] %>% unlist() %>% as.vector(),collapse = ",")
str_c(df$gene[which(df$cluster=="Tumor cells")] %>% unlist() %>% as.vector(),collapse = ",")
str_c(df$gene[which(df$cluster=="Oligodendrocytes")] %>% unlist() %>% as.vector(),collapse = ",")
str_c(df$gene[which(df$cluster=="NKT-like cells")] %>% unlist() %>% as.vector(),collapse = ",")
str_c(df$gene[which(df$cluster=="B cells")] %>% unlist() %>% as.vector(),collapse = ",")
str_c(df$gene[which(df$cluster=="T cells")] %>% unlist() %>% as.vector(),collapse = ",")

#Sampling of GSE182109--------------------------------------------------------------------------------------------------------------------------------------------
Celltype <- c("Proliferating tumor cells","Microglia-like GAMs","Macrophage-like GAMs","NKT-like cells","Tumor cells","T cells","B cells","Oligodendrocytes",
              "Dendritic cells","Endothelial cells","Mural cells")
seurat_ob_SCT_hmy <- subset(x = seurat_ob_SCT_hmy, subset = CellType1 %in% Celltype)
Numb <- lapply(Celltype, function(x)length(which(GSE182109_metadata$CellType1==x))) %>% unlist() %>% as.vector()

single_dataset <- seurat_ob_SCT_hmy@reductions[["umap"]]@cell.embeddings %>% as.data.frame()
single_dataset$UMAP_1 <- single_dataset$UMAP_1 + 20
single_dataset$UMAP_2 <- single_dataset$UMAP_2 + 20
single_dataset$Celltype <- seurat_ob_SCT_hmy@meta.data$CellType1

Final_sample <- c()
samplingNumer <- 1000  #Cell number for each targeted cell type

for (celltype_num in c(1:11)) {
  single_dataset_cell <- single_dataset[which(single_dataset$Celltype==Celltype[celltype_num]),] #單一dataset的資料及特定細胞種類
  sample_pool <- row.names(single_dataset_cell) #targeted cell ID
  X_TotalGM <- exp(mean(log(single_dataset_cell$UMAP_1)))
  Y_TotalGM <- exp(mean(log(single_dataset_cell$UMAP_2)))
  single_dataset_cell$Distance <- ((X_TotalGM-single_dataset_cell$UMAP_1)^2 + (Y_TotalGM-single_dataset_cell$UMAP_2)^2)^0.5
  res.df <- single_dataset_cell %>%
    top_n(-samplingNumer, Distance)
  Final_sample <- append(Final_sample,row.names(res.df))
  print(paste0(Celltype[celltype_num],"GSE182109, 1000顆細胞加入Final_sample"))
  remove(res.df,X_TotalGM,Y_TotalGM,sample_pool)
}

a <- data.frame(cell_ID = Final_sample)
write.csv(a, file = "~/Cell_sampling_11celltype_GSE182109_DatasetWhole_MinDist.csv")

#--------------------------------------------------------------------------------------------
Celltype <- c("Proliferating tumor cells","Microglia-like GAMs","Macrophage-like GAMs","NKT-like cells","Tumor cells","T cells","B cells","Oligodendrocytes",
              "Dendritic cells","Endothelial cells","Mural cells")
FAM <- "611"  #211,411,611
Sampling <- "DatasetWhole_MinDist"  #ALL,DatasetWhole_MinDist
gene_n <- 20  #20,50,100
ex_matrix <- sob[["SCT"]]@data %>% as.data.frame() #ALL=seurat_ob_SCT_hmy, MinDist=sob

cluster.markers <- read_csv(paste0("data/FeatureGeneSet/Findallmarker_GSE182109_11celltype_",FAM,"_",Sampling,".csv")) %>% as.data.frame()
a <- cluster.markers %>% group_by(cluster) %>% summarise(n = length(which(p_val_adj==0)))
for (i in 1:length(Celltype)) {
  if(i==1 & a$n[which(a$cluster==Celltype[i])]<gene_n){
    b <- cluster.markers[which(cluster.markers$cluster==Celltype[i]),] %>% top_n(-gene_n*4, p_val_adj)
  }else if(i==1 & a$n[which(a$cluster==Celltype[i])]>gene_n){
    b <- cluster.markers[which(cluster.markers$cluster==Celltype[i]),]
  }else if(i>1 & a$n[which(a$cluster==Celltype[i])]<gene_n){
    b <- rbind(b,cluster.markers[which(cluster.markers$cluster==Celltype[i]),] %>% top_n(-gene_n*4, p_val_adj))
  }else if(i>1 & a$n[which(a$cluster==Celltype[i])]>gene_n){
    b <- rbind(b,cluster.markers[which(cluster.markers$cluster==Celltype[i]),])
  }
}
remove(a,i)

target_celltype <- b$cluster %>% unique()
for (i in 1:length(target_celltype)) {
  target_gene <- ex_matrix[b$gene[which(b$cluster==target_celltype[i])],] %>% t() %>% as.data.frame()
  target_gene$celltype <- sob@meta.data$CellType1  #改data
  target_gene$celltype[which(target_gene$celltype!=target_celltype[i])] <- "All_others"
  gene_mean_ex <- target_gene %>% group_by(celltype) %>% summarise_all(mean) %>% t() %>%  as.data.frame()
  if(i==1){
    gene_mean_ex_total <- gene_mean_ex[-1,]
  }else{
    gene_mean_ex_total <- rbind(gene_mean_ex_total,gene_mean_ex[-1,])
  }
}
colnames(gene_mean_ex_total) <- c("target_expression_mean","others_expression_mean")
b <- cbind(b,gene_mean_ex_total) %>% as.data.frame()
remove(ex_matrix,gene_mean_ex,gene_mean_ex_total,i,target_celltype,target_gene)

b$pct.diff <- b$pct.1 - b$pct.2
df <- b %>% group_by(cluster) %>%
  top_n(gene_n*4, target_expression_mean) %>%  #(p_val_adj many are 0)
  top_n(gene_n*2, pct.diff) %>%
  top_n(gene_n, avg_log2FC)
write.csv(df, file = paste0("~/Findallmarker_GSE182109_11celltype_",FAM,"_",Sampling,"_n",gene_n,".csv"))
remove(b,cluster.markers,df,Celltype,FAM,gene_n,Sampling)

#--------------------------------------------------------------------------------
TMM <- read_csv("data/GSE182109_norm_tmm_df.csv") %>% as.data.frame()
row.names(TMM) <- TMM$...1
TMM$...1 <- NULL
ex_matrix <- TMM[,which(colnames(TMM) %in% colnames(seurat_ob_SCT_hmy))] %>% as.data.frame() #TMM for ALL
target <- colnames(seurat_ob_SCT_hmy)
ex_matrix <- ex_matrix[,match(target, colnames(ex_matrix))] %>% as.data.frame()

ex_matrix <- TMM[,which(colnames(TMM) %in% colnames(sob))] %>% as.data.frame() #TMM for DatasetWhole_MinDist
target <- colnames(sob)
ex_matrix <- ex_matrix[,match(target, colnames(ex_matrix))] %>% as.data.frame()

TPM <- read_csv("data/GSE182109_TPM_nor.csv") %>% as.data.frame()
row.names(TPM) <- TPM$...1
TPM$...1 <- NULL
ex_matrix <- TPM[,which(colnames(TPM) %in% colnames(seurat_ob_SCT_hmy))] %>% as.data.frame() #TPM for ALL
target <- colnames(seurat_ob_SCT_hmy)
ex_matrix <- ex_matrix[,match(target, colnames(ex_matrix))] %>% as.data.frame()

ex_matrix <- TPM[,which(colnames(TPM) %in% colnames(sob))] %>% as.data.frame() #TPM for DatasetWhole_MinDist
target <- colnames(sob)
ex_matrix <- ex_matrix[,match(target, colnames(ex_matrix))] %>% as.data.frame()


Meth <- "SCT" #RawCount, LogNormalize, SCT, TMM, TPM
Sampling <- "ALL"  #ALL,DatasetWhole_MinDist
#ALL=seurat_ob_SCT_hmy, DatasetWhole_MinDist=sob
#ex_matrix <- seurat_ob_SCT_hmy[["SCT"]]@data %>% as.data.frame()
ex_matrix <- seurat_ob_SCT_hmy[["SCT"]]@scale.data %>% as.data.frame()
#ex_matrix <- seurat_ob_SCT_hmy[["SCT"]]@counts %>% as.data.frame()

FAM <- "611"  #211,411,611
for (gene_n in c(20,50,100)) {
  FM <- read_csv(paste0("data/FeatureGeneSet/Findallmarker_GSE182109_11celltype_",FAM,"_",Sampling,"_n",gene_n,".csv")) %>% as.data.frame()
  target_gene <- FM$gene
  df <- ex_matrix[which(row.names(ex_matrix) %in% target_gene),] %>% t() %>% as.data.frame()
  df$celltype <- seurat_ob_SCT_hmy@meta.data$CellType1    #改檔案
  df2 <- df %>% group_by(celltype) %>% summarise_all(mean) %>% t() %>% as.data.frame()
  colnames(df2) <- df2[1,]
  df2 <- df2[-1,]
  colnames(df2)[grep(pattern = "Proliferating",colnames(df2))] <- "Proliferating GAMs/tumor cells"
  write.csv(df2, file = paste0("data/Reference_20230529/Reference_GSE182109_11celltype_",FAM,"_",Sampling,"_n",gene_n,"_",Meth,".csv"))
  remove(gene_n,FM,target_gene,df,df2)
}
#-------------------------------------------------------------------------------------------
#細胞數量
seurat_ob_SCT_hmy@meta.data$cell_name <- colnames(seurat_ob_SCT_hmy)
cell_num <- c(100,1000,2500,6000)

#Tumor cell
dataset <- "GSE131928"
cell_num <- 100  #6000,2500,1000,100
perc <- seq(from = 60, to = 100, by = 2)

RaTio <- c(0.2,0.4,0.6,0.8)
TumorCell_ID <- seurat_ob_SCT_hmy@meta.data$cell_name[seurat_ob_SCT_hmy@meta.data$CellType=="Tumor cells"]
GAMCell_ID <- seurat_ob_SCT_hmy@meta.data$cell_name[seurat_ob_SCT_hmy@meta.data$CellType %in% c("Microglia-like GAMs", "Macrophage-like GAMs")]
OtherCell_ID <- seurat_ob_SCT_hmy@meta.data$cell_name[seurat_ob_SCT_hmy@meta.data$CellType %in% c("Proliferating tumor cells","NKT-like cells","T cells","B cells",
                                                                                                  "Oligodendrocytes","Dendritic cells","Endothelial cells","Mural cells")]

a <- data.frame(cell_ID = seq(from = 1, to = cell_num, by = 1))

if(dataset=="GSE182109"){
  Nn <- 12
}else{
  Nn <- 6
}
for (i in perc) {
  TG_num <- cell_num*i*0.01
  Other_num <- cell_num-TG_num
  print(paste0("START--",i,"% Tumor/GAM=",TG_num," Others=",Other_num))
  for (j in RaTio) {
    print(paste0("START Raio=",j))
    TumorCell_num <- round(TG_num*j, digits=0) %>% as.numeric()
    GAMCell_num <- TG_num-TumorCell_num
    print(GAMCell_num)
    for (k in 1:Nn) {
      print(k)
      Tumor_ID <- sample(TumorCell_ID,TumorCell_num,replace=TRUE)
      GAM_ID <- sample(GAMCell_ID,GAMCell_num,replace=TRUE)
      Others_ID <- sample(OtherCell_ID,Other_num,replace=TRUE)
      Total_cell_ID <- append(Tumor_ID,GAM_ID)
      Total_cell_ID <- append(Total_cell_ID,Others_ID)
      a <- a %>% add_column(Total_cell_ID)
      remove(Tumor_ID,GAM_ID,Others_ID,Total_cell_ID)
    }
    print(paste0("Completed Raio=",j))
    remove(TumorCell_num,GAMCell_num)
  }
  print(paste0("Completed--",i))
  remove(TG_num,Other_num)
}
write.csv(a, file = paste0("data/BulkRNA_noPro/",dataset,"_BulkRNA_SampleID_N",cell_num,".csv"))
#--------------------
scData <- seurat_ob_SCT_hmy

ex_matrix <- scData[["SCT"]]@counts %>% as.data.frame()

#enter here
dataset <- "GSE131928"
cell_num <- 100 

a <- read_csv(paste0("data/BulkRNA_noPro/",dataset,"_BulkRNA_SampleID_N",cell_num,".csv")) %>% as.data.frame()
a$...1 <- NULL
row.names(a) <- a$cell_ID
a$cell_ID <- NULL
if(dataset=="GSE182109"){
  NN <- 1008
}else{
  NN <- 504
}
colnames(a) <- lapply(c(1:NN), function(x)paste0("GBM",x)) %>% unlist() %>% as.vector()
temp_table <- data.frame(gene_name = row.names(ex_matrix))

for (i in c(1:length(colnames(a)))) {
  temp_table <- temp_table %>% add_column(ex_matrix[,colnames(ex_matrix) %in% a[,i]] %>% rowSums())
}
row.names(temp_table) <- temp_table$gene_name
temp_table$gene_name <- NULL
colnames(temp_table) <- lapply(c(1:NN), function(x)paste0("GBM",x)) %>% unlist() %>% as.vector()
write.csv(temp_table, file = paste0("data/BulkRNA_noPro/",dataset,"_BulkRNA_RawCount_N",cell_num,".csv"))
print(paste0("data/BulkRNA_noPro/",dataset,"_BulkRNA_RawCount_N",cell_num,".csv"))
remove(temp_table,i)
#--------------------
celln <- 100/cell_num

metadata <- data.frame(Cell_ID = scData@meta.data$cell_name,
                       Celltype = scData@meta.data$CellType1) #注意這邊要改
temp_table <- data.frame(Celltype = levels(metadata$Celltype))

for (i in c(1:length(colnames(a)))) {
  b <- metadata[metadata$Cell_ID %in% a[,i],]
  temp_table <- temp_table %>% add_column(summary(b$Celltype)*celln)
}
row.names(temp_table) <- temp_table$Celltype
temp_table$Celltype <- NULL

colnames(temp_table) <- lapply(c(1:NN), function(x)paste0("GBM",x)) %>% unlist() %>% as.vector()
write.csv(temp_table, file = paste0("data/BulkRNA_20230601/",dataset,"_BulkRNA_Celltype_N",cell_num,".csv"))
print(paste0("data/BulkRNA_20230601/",dataset,"_BulkRNA_Celltype_N",cell_num,".csv"))
remove(b,metadata,celln,i)


a <- read_csv("data/BulkRNA_20230601/GSE131928_BulkRNA_Celltype_N100.csv") %>% as.data.frame()
row.names(a) <- a$...1
a$...1 <- NULL
colnames(a) <- lapply(c(1:length(colnames(a))), function(x)paste0("cGBM",x)) %>% unlist() %>% as.vector()
write.csv(a, file = "data/BulkRNA_20230601/GSE131928_BulkRNA_Celltype_N100.csv")

#BulkRNA transformation------------------------------------------------------------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("geneplotter")

library(geneplotter)
library(DESeq2)

dataset <- "MergeDataset" #MergeDataset, GSE182109
NN <- 6000  #100,1000,2500,6000

CountMatrix <- read_csv(paste0("data/BulkRNA_noPro/",dataset,"_BulkRNA_RawCount_N",NN,".csv")) %>% as.data.frame()
#row.names(CountMatrix) <- CountMatrix$geneID
row.names(CountMatrix) <- CountMatrix$...1
#CountMatrix$geneID <- NULL
CountMatrix$...1 <- NULL
coldata <- data.frame(cell_name = colnames(CountMatrix),
                      dataset_name = str_extract(colnames(CountMatrix), pattern = "[:alpha:]+"))

dds <- DESeqDataSetFromMatrix(countData = CountMatrix,
                              colData = coldata,
                              design= ~ 1)
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]
dds <- DESeq(dds, parallel = T)
vsd <- vst(dds, blind=TRUE)
vsd_matrix <- assay(vsd) %>% as.data.frame()
write.csv(vsd_matrix, file = paste0("data/BulkRNA_noPro/",dataset,"_BulkRNA_vst_N",NN,".csv"))

normcount <- counts(dds,normalized=T) %>% as.data.frame()
write.csv(normcount, file = paste0("data/BulkRNA_noPro/",dataset,"_BulkRNA_NormCount_N",NN,".csv"))

ntd <- normTransform(dds)
ntd_matrix <- assay(ntd) %>% as.data.frame()
write.csv(ntd_matrix, file = paste0("data/BulkRNA_noPro/",dataset,"_BulkRNA_ntd_N",NN,".csv"))

remove(dataset,NN,CountMatrix,coldata,dds,vsd,vsd_matrix,ntd,ntd_matrix,normcount)

#rlogd <- rlog(dds, blind=TRUE)
#rlogd_matrix <- assay(rlogd) %>% as.data.frame()
#write.csv(rlogd_matrix, file = "data/BulkRNA_20230601/MergeDataset_BulkRNA_rlog_N1000.csv")

#deconvolution-----------------------------------------------------------------------------------------------------------------------------------------
install.packages("remotes")
remotes::install_github("omnideconv/immunedeconv")

library(immunedeconv)
library(magrittr)
library(tidyverse)
library(stringr)

#deconvolution method
decon_Meth <- "EPIC"   #EPIC

#for expression
Dataset <- "MergeDataset"  #MergeDataset, GSE182109
norMeth_bulk <- "NormCount"  #TPM, TMM, vst, NormCount, RawCount, ntd, logTPM, RPKM
cell_num <- 100   #100,1000,2500,6000

exp_matrix <- read_csv(paste0("data/BulkRNA_noPro/",Dataset,"_BulkRNA_",norMeth_bulk,"_N",cell_num,".csv")) %>% as.data.frame()
row.names(exp_matrix) <- exp_matrix$...1
exp_matrix$...1 <- NULL

#for references
FAM <- c(211,411,611)  #211,411,611
num <- c(20,50,100)   #20, 50, 100
norMeth <- c("TPM","TMM","LogNormalize","SCT","RawCount")

Sampling <- "ALL"  #ALL, DatasetWhole_MinDist

for (i in FAM) {
  for (j in num) {
    for (k in norMeth) {
      ref_matrix <- read_csv(paste0("data/Reference_noPro/Reference_GSE182109_11celltype_",i,"_",Sampling,"_n",j,"_",k,".csv"),
                             show_col_types = FALSE) %>% as.data.frame()
      row.names(ref_matrix) <- ref_matrix$...1
      ref_matrix$...1 <- NULL
      
      a <- deconvolute_epic_custom(exp_matrix,
                                   signature_matrix = ref_matrix,
                                   signature_genes = row.names(ref_matrix))
      a <- a %>% as.data.frame()
      b <- round(a*100, digits = 2)
      write.csv(b, file = paste0("data/Deconvolution_noPro/DeconRes_ref_",i,Sampling,j,k,"_BulkRNA_",Dataset,cell_num,norMeth_bulk,"_",decon_Meth,".csv"))
      print(paste0("data/Deconvolution_noPro/DeconRes_ref_",i,Sampling,j,k,"_BulkRNA_",Dataset,cell_num,norMeth_bulk,"_",decon_Meth,".csv"))
      remove(a,b,ref_matrix)
    }
  }
}
#deconvolution method-mclapply---------------------------------------------------------------------------------------------------------
library(immunedeconv)
library(magrittr)
library(tidyverse)
library(stringr)
library(parallel)

decon_Meth <- "EPIC"   #EPIC

#for expression
Dataset <- "MergeDataset"  #MergeDataset, GSE182109
norMeth_bulk <- "RPKM"  #TPM, TMM, vst, NormCount, RawCount, ntd, logTPM, RPKM

#cell_num <- c(100,1000,2500,6000)   #100,1000,2500,6000
cell_num <- c(2500)   #100,1000,2500,6000
#for references
FAM <- c(211,411,611)  #211,411,611
num <- c(20,50,100)   #20, 50, 100
norMeth <- c("TPM","TMM","LogNormalize","SCT","RawCount")
Sampling <- c("DatasetWhole_MinDist","ALL")

mclapply(cell_num, function(CELL_NUM){
  exp_matrix <- read_csv(paste0("~/",Dataset,"_BulkRNA_",norMeth_bulk,"_N",CELL_NUM,".csv"),
                         show_col_types = FALSE) %>% as.data.frame()
  row.names(exp_matrix) <- exp_matrix$...1
  exp_matrix$...1 <- NULL
  lapply(FAM, function(fam){
    lapply(num, function(NUM){
      mclapply(norMeth, function(NORMETH){
        lapply(Sampling, function(SAMPLING){
          ref_matrix <- read_csv(paste0("data/Reference_noPro/Reference_GSE182109_11celltype_",fam,"_",SAMPLING,"_n",NUM,"_",NORMETH,".csv"),
                                 show_col_types = FALSE) %>% as.data.frame()
          row.names(ref_matrix) <- ref_matrix$...1
          ref_matrix$...1 <- NULL
          
          a <- deconvolute_epic_custom(exp_matrix,
                                       signature_matrix = ref_matrix,
                                       signature_genes = row.names(ref_matrix))
          a <- a %>% as.data.frame()
          b <- round(a*100, digits = 2)
          write.csv(b, file = paste0("data/Deconvolution_noPro/DeconRes_ref_",fam,SAMPLING,NUM,NORMETH,"_BulkRNA_",Dataset,CELL_NUM,norMeth_bulk,"_",decon_Meth,".csv"))
          print(paste0("data/Deconvolution_noPro/DeconRes_ref_",fam,SAMPLING,NUM,NORMETH,"_BulkRNA_",Dataset,CELL_NUM,norMeth_bulk,"_",decon_Meth,".csv"))
          remove(a,b,ref_matrix)
        })
      },mc.cores=5)
    })
  })
},mc.cores=3)


#deconvolution method-------------------------------------------------------------------------------------------------------------
decon_Meth <- "ConsesnusTME"   #EPIC, ConsesnusTME

#for expression
Dataset <- "MergeDataset"  #MergeDataset, GSE182109
norMeth_bulk <- "logTPM"  #TPM, TMM, vst, NormCount, RawCount, ntd, logTPM, rpkm
cell_num <- 6000   #100,1000,2500,6000

exp_matrix <- read_csv(paste0("data/BulkRNA_noPro/",Dataset,"_BulkRNA_",norMeth_bulk,"_N",cell_num,".csv"),
                       show_col_types = FALSE) %>% as.data.frame()
row.names(exp_matrix) <- exp_matrix$...1
exp_matrix$...1 <- NULL
exp_matrix <- data.matrix(exp_matrix)

#for references
Sampling <- c("DatasetWhole_MinDist")  #ALL, DatasetWhole_MinDist

FAM <- c(211,411,611)  #211,411,611
#FAM <- c(411,611)  #211,411,611
num <- c(20,50,100)   #20, 50, 100
#num <- c(100)   #20, 50, 100
norMeth <- c("TPM","TMM","LogNormalize","SCT","RawCount")
#norMeth <- c("TPM")

for (i in FAM) {
  for (j in num) {
    for (k in norMeth) {
      ref_matrix <- read_csv(paste0("data/FeatureGeneSet/Findallmarker_GSE182109_11celltype_",FAM,"_",Sampling,"_n",num,".csv"),
                             show_col_types = FALSE) %>% as.data.frame()
      #      ref_matrix$cluster[grep(pattern = "Proliferating",ref_matrix$cluster)] <- "Proliferating GAMs/tumor cells"
      ref_matrix <- ref_matrix[-grep(pattern = "Proliferating",ref_matrix$cluster),]
      refgene_list <-lapply(unique(ref_matrix$cluster),function(x){
        a <- ref_matrix$gene[ref_matrix$cluster==x]
        return(a)
      })
      names(refgene_list) <- unique(ref_matrix$cluster)
      
      aa <- deconvolute_consensus_tme_custom(exp_matrix,
                                             signature_genes = refgene_list,
                                             stat_method = "ssgsea")
      aa <- aa %>% as.data.frame()
      bb <- colSums(aa)
      temp_a <- data.frame(celltype = row.names(aa))
      
      for(m in 1:ncol(aa)){
        temp_a <- temp_a %>% add_column(aa[,m]/bb[m])
      }
      row.names(temp_a) <- temp_a$celltype
      temp_a$celltype <- NULL
      colnames(temp_a) <- colnames(aa)
      temp_a <- round(temp_a*100, digits = 2)
      write.csv(temp_a, file = paste0("data/Deconvolution_noPro/DeconRes_ref_",i,Sampling,j,k,"_BulkRNA_",Dataset,cell_num,norMeth_bulk,"_",decon_Meth,".csv"))
      print(paste0("data/Deconvolution_noPro/DeconRes_ref_",i,Sampling,j,k,"_BulkRNA_",Dataset,cell_num,norMeth_bulk,"_",decon_Meth,".csv"))
      remove(ref_matrix,refgene_list,aa,bb,temp_a)
    }
  }
}

#deconvolution method: Cibersort----------------------------------------------------------------------------------------------------------
source("~/data/CIBERSORT_ian.R")
decon_Meth <- "Cibersort"   

#for expression
Dataset <- "MergeDataset"  #MergeDataset, GSE182109
norMeth_bulk <- "RPKM"  #TPM, TMM, vst, NormCount, RawCount, ntd, logTPM, rpkm
cell_num <- 100   #100,1000,2500,6000

exp_matrix <- read_csv("data/exome_G70_Rawcounts.csv",
                       show_col_types = FALSE) %>% as.data.frame()
#exp_matrix <- read_csv(paste0("data/BulkRNA_20230601/",Dataset,"_BulkRNA_",norMeth_bulk,"_N",cell_num,".csv"),
#                       show_col_types = FALSE) %>% as.data.frame()
#row.names(exp_matrix) <- exp_matrix$...1
row.names(exp_matrix) <- exp_matrix$gene_name
exp_matrix$...1 <- NULL
exp_matrix$gene_name <- NULL
exp_matrix[is.na(exp_matrix)] <- 0
exp_matrix <- exp_matrix[!is.infinite(rowSums(exp_matrix)),]

#for references
Sampling <- c("ALL")  #ALL, DatasetWhole_MinDist

#FAM <- c(211,411,611)  #211,411,611
FAM <- c(611)  #211,411,611
#num <- c(20,50,100)   #20, 50, 100
num <- c(100)   #20, 50, 100
#norMeth <- c("TPM","TMM","LogNormalize","SCT","RawCount")
norMeth <- c("SCT")

i <- FAM
j <- num
k <- norMeth

for (i in FAM) {
  for (j in num) {
    for (k in norMeth) {
      ref_matrix <- read_csv(paste0("data/Reference_noPro/Reference_GSE182109_11celltype_",i,"_",Sampling,"_n",j,"_",k,".csv"),
                             show_col_types = FALSE) %>% as.data.frame()
      row.names(ref_matrix) <- ref_matrix$...1
      ref_matrix$...1 <- NULL
      
      aa <- CIBERSORT_ian(mixture_file = exp_matrix,
                          sig_matrix = ref_matrix,
                          perm = 100,
                          QN = FALSE,
                          absolute = FALSE,
                          abs_method = "sig.score")
      
      aa <- aa[1:11,] %>% as.data.frame()
      temp_a <- round(aa*100, digits = 2)
      write.csv(temp_a, file = paste0("data/Deconvolution_result_new/DeconRes_ref_",i,Sampling,j,k,"_BulkRNA_",Dataset,cell_num,norMeth_bulk,"_",decon_Meth,".csv"))
      print(paste0("data/Deconvolution_result_new/DeconRes_ref_",i,Sampling,j,k,"_BulkRNA_",Dataset,cell_num,norMeth_bulk,"_",decon_Meth,".csv"))
      remove(ref_matrix,aa,temp_a)
    }
  }
}

#----------------------------------------------
install.packages('Metrics')
library(Metrics)

#for bulk true celltype-----------------------------------------------------------------------------------------------------------
library(ggpubr)

Dataset_true <- "MergeDataset"  #MergeDataset, GSE182109
Dataset <- Dataset_true
norMeth_bulk <- c("TPM", "TMM", "vst", "NormCount", "RawCount", "ntd", "logTPM", "RPKM")  #TPM, TMM, vst, NormCount, RawCount, ntd, logTPM, rpkm
#norMeth_bulk <- c("RawCount")  #TPM, TMM, vst, NormCount, RawCount, ntd, logTPM, RPKM)
Sampling <- c("ALL","DatasetWhole_MinDist")  #ALL, DatasetWhole_MinDist
#Sampling <- c("DatasetWhole_MinDist")  #ALL, DatasetWhole_MinDist
norMeth <- c("TPM","TMM","LogNormalize","SCT","RawCount")   #TMM, TPM, LogNormalize, SCT, RawCount
FAM <- c(211,411,611)  #211,411,611
#FAM <- c(611)  #211,411,611
num <- c(20,50,100)   #20, 50, 100
#num <- c(50,100)   #20, 50, 100

decon_Meth <- "EPIC"   #EPIC, ConsesnusTME, CiberAbs, Cibersort
cell_num_true <- 2500   #100,1000,2500,6000
cell_num <- cell_num_true   #100,1000,2500,6000

True_res <- read_csv(paste0("data/BulkRNA_noPro/",Dataset_true,"_BulkRNA_Celltype_N",cell_num_true,".csv"),
                     show_col_types = FALSE) %>% as.data.frame()
True_res <- True_res[-which(True_res$Celltype=="Others"),] %>% as.data.frame()
row.names(True_res) <- True_res$Celltype
True_res$...1 <- NULL
True_res$Celltype <- NULL
True_res[is.na(True_res)==TRUE] <- 0
True_res[nrow(True_res) + 1,] <- 0
row.names(True_res)[nrow(True_res)] <- "otherCells"   #EPIC: otherCells

for (i in norMeth_bulk) {
  for (j in FAM) {
    for (k in num) {
      for (l in norMeth) {
        for (m in Sampling) {
          Decon_res <- read_csv(paste0("data/Deconvolution_noPro/DeconRes_ref_",j,m,k,l,"_BulkRNA_",Dataset,cell_num,i,"_",decon_Meth,".csv"),
                                show_col_types = FALSE) %>% as.data.frame()
          row.names(Decon_res) <- Decon_res$...1
          Decon_res$...1 <- NULL
          Decon_res$Length <- NULL
          
          target <- row.names(True_res)
          Decon_res <- Decon_res[match(target, row.names(Decon_res)),] %>% as.data.frame()
          row.names(Decon_res) <- row.names(True_res)
          Decon_res[is.na(Decon_res)==TRUE] <- 0
          
          temp_df_a <- data.frame(Celltype = c(row.names(True_res),"SumUp"))
          temp_df_pval <- data.frame(Celltype = c(row.names(True_res),"SumUp"))
          temp_df_corr <- data.frame(Celltype = c(row.names(True_res),"SumUp"))
          
          for (ds_name in c(".", "NC", "cGBM", "^GBM")) {
            True_res_dataset <- True_res[,grep(pattern = ds_name,colnames(True_res))]   
            Decon_res_dataset <- Decon_res[,grep(pattern = ds_name,colnames(Decon_res))]   
            a <- lapply(1:nrow(True_res_dataset),function(x){
              actual    <- True_res_dataset[x,] %>% unlist() %>% as.vector()
              predicted <- Decon_res_dataset[x,] %>% unlist() %>% as.vector()
              if (length(actual) != length(predicted)) {
                print('The legnth of two array is not equal')
              }
              RMSE <- rmse(actual, predicted)
              return(RMSE)
            }) %>% unlist() %>% as.vector()
            a <- append(a,sqrt(sum((a^2)*ncol(True_res_dataset))/(ncol(True_res_dataset)*nrow(True_res_dataset))))
            temp_df_a <- temp_df_a %>% add_column(a)
            
            b <- lapply(1:nrow(True_res_dataset),function(x){
              actual    <- True_res_dataset[x,] %>% unlist() %>% as.vector()
              predicted <- Decon_res_dataset[x,] %>% unlist() %>% as.vector()
              if (length(actual) != length(predicted)) {
                print('The legnth of two array is not equal')
              }
              res <- cor.test(actual, predicted, 
                              method = "pearson")
              pval <- res$p.value
              corR <- res$estimate
              sta <- append(pval,corR)
              return(sta)
            })
            actual <- as.vector(unlist(True_res_dataset))
            predicted <- as.vector(unlist(Decon_res_dataset))
            res <- cor.test(actual, predicted, 
                            method = "pearson")
            temp_b <- c()
            b_pval <- lapply(1:length(b),function(x){
              temp_b <- append(temp_b,b[[x]][1])
              return(temp_b)
            }) %>% unlist() %>% as.vector()
            b_pval <- append(b_pval,res$p.value)
            temp_df_pval <- temp_df_pval %>% add_column(b_pval)
            temp_b <- c()
            b_corr <- lapply(1:length(b),function(x){
              temp_b <- append(temp_b,b[[x]][2])
              return(temp_b)
            }) %>% unlist() %>% as.vector()
            b_corr <- append(b_corr,res$estimate)
            temp_df_corr <- temp_df_corr %>% add_column(b_corr)
            
            remove(True_res_dataset,Decon_res_dataset,a,b,b_pval,b_corr,res,actual,predicted)
          }
          remove(target,Decon_res)
          
          temp_df_b <- cbind(temp_df_corr,temp_df_pval[,-1]) %>% as.data.frame()
          final_temp <- merge(x=temp_df_a,y=temp_df_b, by.x = "Celltype", by.y = "Celltype", all.x = TRUE) %>% as.data.frame()
          colnames(final_temp) <- c("Celltype",
                                    "MergeDataset_RMSE","NC_RMSE","GSE131928_RMSE","GSE182109_RMSE",
                                    "MergeDataset_PearsonR","NC_PearsonR","GSE131928_PearsonR","GSE182109_PearsonR",
                                    "MergeDataset_CorP","NC_CorP","GSE131928_CorP","GSE182109_CorP")
          target <- temp_df_a$Celltype
          final_temp <- final_temp[match(target, final_temp$Celltype),] %>% as.data.frame()
          row.names(final_temp ) <- final_temp$Celltype
          final_temp$Celltype <- NULL
          
          write.csv(final_temp, file = paste0("data/RMSE_noPro/DeconStat_ref_",j,m,k,l,"_BulkRNA_",Dataset,cell_num,i,"_",decon_Meth,".csv"))
          print(paste0("data/RMSE_noPro/DeconStat_ref_",j,m,k,l,"_BulkRNA_",Dataset,cell_num,i,"_",decon_Meth,".csv"))
          remove(temp_df_a,temp_df_b,temp_df_corr, temp_df_pval,ds_name,target,final_temp)
        }
      }
    }
  }
}

list.files(path = "data/RMSE_noPro/", pattern = "MergeDataset100TPM_EPIC")

#conclusion---------------------------------------------------------------------------------------------------
library(magrittr)
library(stringr)
library(tidyverse)
library(openxlsx)

Dataset_true <- "MergeDataset"
Dataset <- Dataset_true
decon_Meth <- c("EPIC","ConsesnusTME", "Cibersort")  #EPIC, ConsesnusTME, CiberAbs, Cibersort
Sampling <- c("ALL","DatasetWhole_MinDist")  #ALL, DatasetWhole_MinDist
norMeth <- c("TPM","TMM","LogNormalize","SCT","RawCount")   #TMM, TPM, LogNormalize, SCT, RawCount
FAM <- c(211,411,611)  #211,411,611
num <- c(20,50,100)   #20, 50, 100
dataset_name <- c("MergeDataset", "NC", "GSE131928", "GSE182109")

norMeth_bulk <- "RPKM"   #TPM, TMM, vst, NormCount, RawCount, ntd, logTPM, RPKM
cell_num_true <- 1000   #100,1000,2500,6000
cell_num <- cell_num_true   #100,1000,2500,6000

for (n in dataset_name) {
  sum_table <- data.frame(Celltype = c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells",
                                       "Oligodendrocytes","T cells","Tumor cells","B cells","Mural cells","otherCells","SumUp"))
  sum_table_r <- data.frame(Celltype = c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells",
                                         "Oligodendrocytes","T cells","Tumor cells","B cells","Mural cells","otherCells","SumUp"))
  sample_name <- c()
  
  for(i in decon_Meth){
    for(j in FAM){
      for (k in num) {
        for (l in norMeth) {
          for (m in Sampling) {
            temp_file <- read_csv(paste0("data/RMSE_noPro/DeconStat_ref_",j,m,k,l,"_BulkRNA_",Dataset,cell_num,norMeth_bulk,"_",i,".csv"),
                                  col_select = c(paste0(n,"_RMSE"),paste0(n,"_PearsonR")),
                                  show_col_types = FALSE)
            sum_table <- sum_table %>% add_column(temp_file[,1] %>% as.vector())
            sum_table_r <- sum_table_r %>% add_column(temp_file[,2] %>% as.vector())
            sample_name <- append(sample_name,paste(i,j,m,k,l,sep = "_"))
          }
        }
      }
    }
  }
  row.names(sum_table) <- sum_table$Celltype
  sum_table$Celltype <- NULL
  colnames(sum_table) <- sample_name
  a <- sum_table %>% t() %>% as.data.frame()
  write.xlsx(a, file = paste0("data/RMSE_noPro/SumTable_BulkRNA_",n,"_",cell_num_true,norMeth_bulk,".xlsx"),
             rowNames=TRUE)
  print(paste0("data/RMSE_noPro/SumTable_BulkRNA_",n,"_",cell_num_true,norMeth_bulk,".xlsx"))
  
  row.names(sum_table_r) <- sum_table_r$Celltype
  sum_table_r$Celltype <- NULL
  colnames(sum_table_r) <- sample_name
  a <- sum_table_r %>% t() %>% as.data.frame()
  write.xlsx(a, file = paste0("data/RMSE_noPro/SumTable_BulkRNA_",n,"_",cell_num_true,norMeth_bulk,"_PearsonR.xlsx"),
             rowNames=TRUE)
  print(paste0("data/RMSE_noPro/SumTable_BulkRNA_",n,"_",cell_num_true,norMeth_bulk,"_PearsonR.xlsx"))
  
  remove(sum_table,sum_table_r,sample_name,a,temp_file)
}

#--------------------------------------------------------------------------------------------------
library(ggplot2)
library(ggpubr)

dataset_name <- c("MergeDataset", "NC", "GSE131928", "GSE182109")

norMeth_bulk <- "TPM"   #TPM, TMM, vst, NormCount, RawCount, ntd, logTPM, RPKM

for (i in dataset_name) {
  temp100_RMSE <- read_excel(paste0("data/RMSE_noPro/SumTable_BulkRNA_",i,"_100",norMeth_bulk,".xlsx")) %>% as.data.frame()
  temp100_RMSE$...1 <- paste0(temp100_RMSE$...1,"_",norMeth_bulk,"100")
  temp100_Pearson <- read_excel(paste0("data/RMSE_noPro/SumTable_BulkRNA_",i,"_100",norMeth_bulk,"_PearsonR.xlsx")) %>% as.data.frame()
  temp100_Pearson$...1 <- paste0(temp100_Pearson$...1,"_",norMeth_bulk,"100")
  
  temp1000_RMSE <- read_excel(paste0("data/RMSE_noPro/SumTable_BulkRNA_",i,"_1000",norMeth_bulk,".xlsx")) %>% as.data.frame()
  temp1000_RMSE$...1 <- paste0(temp1000_RMSE$...1,"_",norMeth_bulk,"1000")
  temp1000_Pearson <- read_excel(paste0("data/RMSE_noPro/SumTable_BulkRNA_",i,"_1000",norMeth_bulk,"_PearsonR.xlsx")) %>% as.data.frame()
  temp1000_Pearson$...1 <- paste0(temp1000_Pearson$...1,"_",norMeth_bulk,"1000")
  
  temp2500_RMSE <- read_excel(paste0("data/RMSE_noPro/SumTable_BulkRNA_",i,"_2500",norMeth_bulk,".xlsx")) %>% as.data.frame()
  temp2500_RMSE$...1 <- paste0(temp2500_RMSE$...1,"_",norMeth_bulk,"2500")
  temp2500_Pearson <- read_excel(paste0("data/RMSE_noPro/SumTable_BulkRNA_",i,"_2500",norMeth_bulk,"_PearsonR.xlsx")) %>% as.data.frame()
  temp2500_Pearson$...1 <- paste0(temp2500_Pearson$...1,"_",norMeth_bulk,"2500")
  
  temp6000_RMSE <- read_excel(paste0("data/RMSE_noPro/SumTable_BulkRNA_",i,"_6000",norMeth_bulk,".xlsx")) %>% as.data.frame()
  temp6000_RMSE$...1 <- paste0(temp6000_RMSE$...1,"_",norMeth_bulk,"6000")
  temp6000_Pearson <- read_excel(paste0("data/RMSE_noPro/SumTable_BulkRNA_",i,"_6000",norMeth_bulk,"_PearsonR.xlsx")) %>% as.data.frame()
  temp6000_Pearson$...1 <- paste0(temp6000_Pearson$...1,"_",norMeth_bulk,"6000")
  
  temp_RMSE <- rbind(temp100_RMSE,temp1000_RMSE,temp2500_RMSE,temp6000_RMSE) %>% as.data.frame()
  remove(temp100_RMSE,temp1000_RMSE,temp2500_RMSE,temp6000_RMSE)
  row.names(temp_RMSE) <- temp_RMSE$...1
  temp_RMSE$...1 <- NULL
  temp_Pearson <- rbind(temp100_Pearson,temp1000_Pearson,temp2500_Pearson,temp6000_Pearson) %>% as.data.frame()
  remove(temp100_Pearson,temp1000_Pearson,temp2500_Pearson,temp6000_Pearson)
  row.names(temp_Pearson) <- temp_Pearson$...1
  temp_Pearson$...1 <- NULL
  
  temp_RMSE$dataset <- i
  temp_Pearson$dataset <- i
  
  if(i=="MergeDataset"){
    temp_RMSE_MergeDataset <- temp_RMSE
    temp_Pearson_MergeDataset <- temp_Pearson
    print("MergeDataset")
  }else if(i=="NC"){
    temp_RMSE_NC <- temp_RMSE
    temp_Pearson_NC <- temp_Pearson
    print("NC")
  }else if(i=="GSE131928"){
    temp_RMSE_GSE131928 <- temp_RMSE
    temp_Pearson_GSE131928 <- temp_Pearson
    print("GSE131928")
  }else if(i=="GSE182109"){
    temp_RMSE_GSE182109 <- temp_RMSE
    temp_Pearson_GSE182109 <- temp_Pearson
    print("GSE182109")
  }
}
dataset_RMSE <- list(temp_RMSE_MergeDataset,temp_RMSE_NC,temp_RMSE_GSE131928, temp_RMSE_GSE182109)
dataset_Pearson <- list(temp_Pearson_MergeDataset,temp_Pearson_NC,temp_Pearson_GSE131928, temp_Pearson_GSE182109)
remove(temp_RMSE_MergeDataset,temp_RMSE_NC,temp_RMSE_GSE131928, temp_RMSE_GSE182109,temp_Pearson_MergeDataset,temp_Pearson_NC,
       temp_Pearson_GSE131928, temp_Pearson_GSE182109,temp_Pearson,temp_RMSE)


Celltype <- colnames(dataset_RMSE[[1]])[1:13]

temp_RMSE <- data.frame(SampleMethod = row.names(dataset_RMSE[[1]]))
temp_Pearson <- data.frame(SampleMethod = row.names(dataset_Pearson[[1]]))
c_name <- c()

for (i in Celltype) {
  for (j in 1:length(dataset_RMSE)) {
    temp_RMSE <- temp_RMSE %>% add_column(dataset_RMSE[[j]][,i])
    temp_Pearson <- temp_Pearson %>% add_column(dataset_Pearson[[j]][,i])
    c_name <- append(c_name, paste(i,dataset_name[j],sep = "_"))
  }
}
colnames(temp_RMSE) <- c("SampleMethod",c_name)
colnames(temp_Pearson) <- c("SampleMethod",c_name)

row.names(temp_RMSE) <- temp_RMSE$SampleMethod
temp_RMSE$SampleMethod <- NULL
row.names(temp_Pearson) <- temp_Pearson$SampleMethod
temp_Pearson$SampleMethod <- NULL

write.xlsx(temp_RMSE, file = paste0("data/RMSE_noPro/SumTable_BulkRNA_RMSE_AllDataset_Celltype_",norMeth_bulk,".xlsx"),
           rowNames=TRUE)
print(paste0("data/RMSE_noPro/SumTable_BulkRNA_RMSE_AllDataset_Celltype_",norMeth_bulk,".xlsx"))
write.xlsx(temp_Pearson, file = paste0("data/RMSE_noPro/SumTable_BulkRNA_Pearson_AllDataset_Celltype_",norMeth_bulk,".xlsx"),
           rowNames=TRUE)
print(paste0("data/RMSE_noPro/SumTable_BulkRNA_Pearson_AllDataset_Celltype_",norMeth_bulk,".xlsx"))
#-------------------------------------
df_all <- temp_df[-grep(pattern = "DatasetWhole_MinDist",row.names(temp_df)),] %>% as.data.frame()
df_dwmd <- temp_df[grep(pattern = "DatasetWhole_MinDist",row.names(temp_df)),] %>% as.data.frame()

a <- str_split(row.names(df_all), pattern = "_")
temp_a <- data.frame(sample = row.names(df_all),
                     decon_meth = lapply(1:length(a),function(x)a[[x]][1]) %>% unlist() %>% as.vector(),
                     fam = lapply(1:length(a),function(x)a[[x]][2]) %>% unlist() %>% as.vector(),
                     sampling = lapply(1:length(a),function(x)a[[x]][3]) %>% unlist() %>% as.vector(),
                     gene_num = lapply(1:length(a),function(x)a[[x]][4]) %>% unlist() %>% as.vector(),
                     ref_normeth = lapply(1:length(a),function(x)a[[x]][5]) %>% unlist() %>% as.vector(),
                     bulk_data = lapply(1:length(a),function(x)a[[x]][6]) %>% unlist() %>% as.vector())
df_all <- cbind(df_all,temp_a) %>% as.data.frame()
df_all$sample <- NULL
remove(a,temp_a)
a <- str_split(row.names(df_dwmd), pattern = "_")
temp_a <- data.frame(sample = row.names(df_dwmd),
                     decon_meth = lapply(1:length(a),function(x)a[[x]][1]) %>% unlist() %>% as.vector(),
                     fam = lapply(1:length(a),function(x)a[[x]][2]) %>% unlist() %>% as.vector(),
                     sampling = "DatasetWhole_MinDist",
                     gene_num = lapply(1:length(a),function(x)a[[x]][5]) %>% unlist() %>% as.vector(),
                     ref_normeth = lapply(1:length(a),function(x)a[[x]][6]) %>% unlist() %>% as.vector(),
                     bulk_data = lapply(1:length(a),function(x)a[[x]][7]) %>% unlist() %>% as.vector())
df_dwmd <- cbind(df_dwmd,temp_a) %>% as.data.frame()
df_dwmd$sample <- NULL
remove(a,temp_a)
temp_df <- rbind(df_all,df_dwmd) %>% as.data.frame()

factor_in <- temp_df$`Mural cells`

pic_ref_normeth <- ggplot(temp_df, aes(x=ref_normeth, y=factor_in, fill=decon_meth)) +
  geom_boxplot(position=position_dodge(0.85)) +
  theme(legend.position="none")

pic_ref_gene_num <- ggplot(temp_df, aes(x=gene_num, y=factor_in, fill=decon_meth)) +
  geom_boxplot(position=position_dodge(0.85)) +
  theme(legend.position="none")

pic_ref_sampling <- ggplot(temp_df, aes(x=sampling, y=factor_in, fill=decon_meth)) +
  geom_boxplot(position=position_dodge(0.85)) +
  theme(legend.position="none")

pic_ref_fam <- ggplot(temp_df, aes(x=fam, y=factor_in, fill=decon_meth)) +
  geom_boxplot(position=position_dodge(0.85)) +
  theme(legend.position="none")

pic_ref_bulk_data <- ggplot(temp_df, aes(x=bulk_data, y=factor_in, fill=decon_meth)) +
  geom_boxplot(position=position_dodge(0.85)) +
  theme(legend.position="none")

ggarrange(pic_ref_normeth, pic_ref_gene_num, pic_ref_sampling, pic_ref_fam, pic_ref_bulk_data,
          ncol = 5, nrow = 1)
#Best result of each bulk transformation------------------------------------------------------------------------------------------------------
library(readxl)

BulkNorm_Meth <- "RawCount"   #"TPM","TMM","vst","NormCount","RawCount","ntd","logTPM","RPKM"

oridata <- read_excel(paste0("data/RMSE_noPro/SumTable_BulkRNA_RMSE_AllDataset_Celltype_",BulkNorm_Meth,".xlsx")) %>% as.data.frame()
print(paste0("open file--data/RMSE_noPro/SumTable_BulkRNA_RMSE_AllDataset_Celltype_",BulkNorm_Meth,".xlsx"))
row.names(oridata) <- oridata$...1
oridata$...1 <- NULL
oridata <- oridata[-grep(pattern = "CibersortEPIC_Best_mix_", row.names(oridata)),] %>% as.data.frame()
tempdf_all <- oridata[-grep(pattern = "DatasetWhole_MinDist",row.names(oridata)),] %>% as.data.frame()
tempdf_dwmd <- oridata[grep(pattern = "DatasetWhole_MinDist",row.names(oridata)),] %>% as.data.frame()
a <- str_split(row.names(tempdf_all), pattern = "_")
b <- str_split(row.names(tempdf_dwmd), pattern = "_")
for (j in 1:nrow(tempdf_all)) {
  tempdf_all[j,c(49:54)] <- a[[j]]
  tempdf_dwmd[j,c(49:54)] <- b[[j]][c(1,2,4,5,6,7)]
}
tempdf <- rbind(tempdf_all,tempdf_dwmd) %>% as.data.frame()
print("Done RMSE")
remove(oridata,tempdf_all,tempdf_dwmd,a,b)

colnames(tempdf)[c(49:54)] <- c("DeconMeth","Min.pct","Sampling","GeneNumber","RefNorm","BulkRNA")
tempdf$BulkNorm <- BulkNorm_Meth
tempdf$BulkCellNumber <- c(rep(100,135),rep(1000,135),rep(2500,135),rep(6000,135),
                           rep(100,135),rep(1000,135),rep(2500,135),rep(6000,135))
tempdf$Value <- c(rep("RMSE",1080))
mergetable <- tempdf[-which(tempdf$BulkCellNumber==100),] %>% as.data.frame()
mergetable$Value <- NULL
mergetable$BulkCellNumber <- NULL
mergetable$BulkRNA <- NULL
mergetable <- mergetable %>% group_by(DeconMeth,Min.pct,Sampling,GeneNumber,RefNorm,BulkNorm) %>% summarise_all(mean) %>% as.data.frame()
View(mergetable[,c(1,2,3,4,5,6,7)])

#Combine results-----------------modified----------------------
decon_Meth <- c("Cibersort","EPIC","ConsesnusTME")
FAM <- c(211,411,611)
Sampling <- c("ALL","DatasetWhole_MinDist")
num <- num <- c(20,50,100)
norMeth <- c("TPM","TMM","LogNormalize","SCT","RawCount")
Celltype <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells",
              "Oligodendrocytes","T cells","Tumor cells","B cells","Mural cells","otherCells")

#bulk RNAseq normalization
norMeth_bulk <- c("RPKM")  #TPM, TMM, vst, NormCount, RawCount, ntd, logTPM, RPKM

#enter c(decon_Meth,FAM,Sampling,num,norMeth)
factorlist <- list(Dendritic_cells = c(1,2,1,3,4),
                   Endothelial_cells = c(1,3,1,2,5),
                   Macrophage_GAMs = c(1,1,2,2,5),
                   Microglia_GAMs = c(1,3,1,3,5),
                   NKT_cells = c(2,1,1,1,2),
                   Oligo = c(1,1,2,3,2),
                   T_cells = c(2,2,1,3,1),
                   Tumor_cells = c(1,3,1,3,2),
                   B_cells = c(1,3,2,3,2),
                   Mural_cells = c(1,1,1,3,5))

for (cell_num in c(100,1000,2500,6000)) {
  filelist <- c()
  for (j in 1:(length(Celltype)-1)) {
    a <- decon_Meth[factorlist[[j]][1]]
    b <- FAM[factorlist[[j]][2]]
    c <- Sampling[factorlist[[j]][3]]
    d <- num[factorlist[[j]][4]]
    e <- norMeth[factorlist[[j]][5]]
    decon_res <- read_csv(paste0("data/Deconvolution_noPro/DeconRes_ref_",b,c,d,e,"_BulkRNA_MergeDataset",cell_num,norMeth_bulk,"_",a,".csv"),
                          show_col_types = FALSE) %>% as.data.frame()
    if(j==1){
      combine_res <- decon_res[which(decon_res$...1==Celltype[j]),]
    }else{
      combine_res <- rbind(combine_res,decon_res[which(decon_res$...1==Celltype[j]),])
    }
    filelist <- append(filelist,paste0(a,b,c,d,e))
    remove(a,b,c,d,e)
  }
  print(filelist)
  combine_res[11,] <- c(Celltype[11],rep(0,length(colnames(combine_res))-1))
  row.names(combine_res) <- combine_res$...1
  combine_res$...1 <- NULL
  re_res <- combine_res
  for (i in 1:length(colnames(combine_res))) {
    combine_res[,i] <- combine_res[,i] %>% as.double()
    re_res[,i] <- NA
    re_res[,i] <- round(combine_res[,i]/sum(combine_res[,i]),digits = 4)*100
  }
  write.csv(re_res, file = paste0("data/Deconvolution_noPro/DeconRes_ref_Best_mix_",cell_num,"_",norMeth_bulk,"_CibersortEPIC_noPro.csv"))
  print(paste0("data/Deconvolution_noPro/DeconRes_ref_Best_mix_",cell_num,"_",norMeth_bulk,"_CibersortEPIC_noPro.csv"))
  remove(filelist,combine_res,re_res,i,j)
}

#add the "mix" result into sumtable---------------------------------------------------------------------------------------------
install.packages('Metrics')
library(Metrics)
library(ggpubr)

Dataset_true <- "MergeDataset"
norMeth_bulk <- c("RPKM")  #TPM, TMM, vst, NormCount, RawCount, ntd, logTPM, RPKM

for (cell_num in c(100,1000,2500,6000)) {
  True_res <- read_csv(paste0("~/",Dataset_true,"_BulkRNA_Celltype_N",cell_num,".csv"),
                       show_col_types = FALSE) %>% as.data.frame()
  True_res <- True_res[-which(True_res$Celltype=="Others"),] %>% as.data.frame()
  row.names(True_res) <- True_res$Celltype
  True_res$...1 <- NULL
  True_res$Celltype <- NULL
  True_res[is.na(True_res)==TRUE] <- 0
  True_res[nrow(True_res) + 1,] <- 0
  row.names(True_res)[nrow(True_res)] <- "otherCells"   #EPIC: otherCells
  
  Decon_res <- read_csv(paste0("~/DeconRes_ref_Best_mix_",cell_num,"_",norMeth_bulk,"_CibersortEPIC_noPro.csv"),
                        show_col_types = FALSE) %>% as.data.frame()
  
  row.names(Decon_res) <- Decon_res$...1
  Decon_res$...1 <- NULL
  Decon_res$Length <- NULL
  
  target <- row.names(True_res)
  Decon_res <- Decon_res[match(target, row.names(Decon_res)),] %>% as.data.frame()
  row.names(Decon_res) <- row.names(True_res)
  Decon_res[is.na(Decon_res)==TRUE] <- 0
  
  temp_df_a <- data.frame(Celltype = c(row.names(True_res),"SumUp"))
  temp_df_pval <- data.frame(Celltype = c(row.names(True_res),"SumUp"))
  temp_df_corr <- data.frame(Celltype = c(row.names(True_res),"SumUp"))
  
  for (ds_name in c(".", "NC", "cGBM", "^GBM")) {
    True_res_dataset <- True_res[,grep(pattern = ds_name,colnames(True_res))]   
    Decon_res_dataset <- Decon_res[,grep(pattern = ds_name,colnames(Decon_res))]   
    a <- lapply(1:nrow(True_res_dataset),function(x){
      actual    <- True_res_dataset[x,] %>% unlist() %>% as.vector()
      predicted <- Decon_res_dataset[x,] %>% unlist() %>% as.vector()
      if (length(actual) != length(predicted)) {
        print('The legnth of two array is not equal')
      }
      RMSE <- rmse(actual, predicted)
      return(RMSE)
    }) %>% unlist() %>% as.vector()
    a <- append(a,sqrt(sum((a^2)*ncol(True_res_dataset))/(ncol(True_res_dataset)*nrow(True_res_dataset))))
    temp_df_a <- temp_df_a %>% add_column(a)
    
    b <- lapply(1:nrow(True_res_dataset),function(x){
      actual    <- True_res_dataset[x,] %>% unlist() %>% as.vector()
      predicted <- Decon_res_dataset[x,] %>% unlist() %>% as.vector()
      if (length(actual) != length(predicted)) {
        print('The legnth of two array is not equal')
      }
      res <- cor.test(actual, predicted, 
                      method = "pearson")
      pval <- res$p.value
      corR <- res$estimate
      sta <- append(pval,corR)
      return(sta)
    })
    actual <- as.vector(unlist(True_res_dataset))
    predicted <- as.vector(unlist(Decon_res_dataset))
    res <- cor.test(actual, predicted, 
                    method = "pearson")
    temp_b <- c()
    b_pval <- lapply(1:length(b),function(x){
      temp_b <- append(temp_b,b[[x]][1])
      return(temp_b)
    }) %>% unlist() %>% as.vector()
    b_pval <- append(b_pval,res$p.value)
    temp_df_pval <- temp_df_pval %>% add_column(b_pval)
    temp_b <- c()
    b_corr <- lapply(1:length(b),function(x){
      temp_b <- append(temp_b,b[[x]][2])
      return(temp_b)
    }) %>% unlist() %>% as.vector()
    b_corr <- append(b_corr,res$estimate)
    temp_df_corr <- temp_df_corr %>% add_column(b_corr)
    
    remove(True_res_dataset,Decon_res_dataset,a,b,b_pval,b_corr,res,actual,predicted)
  }
  remove(target,Decon_res)
  
  temp_df_b <- cbind(temp_df_corr,temp_df_pval[,-1]) %>% as.data.frame()
  final_temp <- merge(x=temp_df_a,y=temp_df_b, by.x = "Celltype", by.y = "Celltype", all.x = TRUE) %>% as.data.frame()
  colnames(final_temp) <- c("Celltype",
                            "MergeDataset_RMSE","NC_RMSE","GSE131928_RMSE","GSE182109_RMSE",
                            "MergeDataset_PearsonR","NC_PearsonR","GSE131928_PearsonR","GSE182109_PearsonR",
                            "MergeDataset_CorP","NC_CorP","GSE131928_CorP","GSE182109_CorP")
  target <- temp_df_a$Celltype
  final_temp <- final_temp[match(target, final_temp$Celltype),] %>% as.data.frame()
  row.names(final_temp ) <- final_temp$Celltype
  final_temp$Celltype <- NULL
  
  write.csv(final_temp, file = paste0("data/RMSE_noPro/DeconStat_ref_Best_mix_",cell_num,"_",norMeth_bulk,"_CibersortEPIC_noPro.csv"))
  print(paste0("data/RMSE_noPro/DeconStat_ref_Best_mix_",cell_num,"_",norMeth_bulk,"_CibersortEPIC_noPro.csv"))
  
  remove(temp_df_a,temp_df_b,temp_df_corr, temp_df_pval,ds_name,target,final_temp,temp_b)
}

#Combine
library(readxl)
cell_num <- c(100,1000,2500,6000)
norMeth_bulk <- c("RPKM")  #TPM, TMM, vst, NormCount, RawCount, ntd, logTPM, RPKM

#SumTable_Celltype <- read_excel(paste0("~/SumTable_BulkRNA_RMSE_AllDataset_Celltype_",norMeth_bulk,".xlsx")) %>% as.data.frame()
SumTable_Celltype <- read_excel(paste0("~/SumTable_BulkRNA_Pearson_AllDataset_Celltype_",norMeth_bulk,".xlsx")) %>% as.data.frame()
row.names(SumTable_Celltype) <- SumTable_Celltype$...1
SumTable_Celltype$...1 <- NULL

elist <- list()
for (j in c(1:4)) {
  DeconStat_ref <- read_csv(paste0("~/DeconStat_ref_Best_mix_",cell_num[j],"_",norMeth_bulk,"_CibersortEPIC_noPro.csv")) %>% as.data.frame()
  row.names(DeconStat_ref) <- DeconStat_ref$...1
  DeconStat_ref$...1 <- NULL
  
  #a <- DeconStat_ref[,c(1:4)] %>% as.data.frame()    #for RMSE
  a <- DeconStat_ref[,c(5:8)] %>% as.data.frame()    #for Pearson
  b <- c()
  for (i in 1:length(row.names(a))) {
    b <- append(b,a[i,] %>% as.vector())
  }
  c <- c("MergeDataset","NC","GSE131928","GSE182109")
  b <- append(b,c) %>% unlist() %>% as.vector()
  elist[[j]] <- b
}
bb <- data.frame(CibersortEPIC_Best_mix_100 = elist[[1]],
                 CibersortEPIC_Best_mix_1000 = elist[[2]],
                 CibersortEPIC_Best_mix_2500 = elist[[3]],
                 CibersortEPIC_Best_mix_6000 = elist[[4]])
colnames(bb) <- paste0("CibersortEPIC_Best_mix_",cell_num,"_",norMeth_bulk)
bb <- bb %>% t() %>% as.data.frame()
colnames(bb) <- colnames(SumTable_Celltype)
for (i in 1:48) {
  bbb <- bb[,i]
  bbb <- round(as.double(bbb), digits = 6)
  bb[,i] <- bbb
}
SumTable_Celltype <- rbind(SumTable_Celltype,bb) %>% as.data.frame()

#View(SumTable_Celltype[grep(pattern = "CibersortEPIC_Best_mix_",row.names(SumTable_Celltype)),c("SumUp_MergeDataset","SumUp_NC")])

#write.xlsx(SumTable_Celltype, file = paste0("data/RMSE_noPro/SumTable_BulkRNA_RMSE_AllDataset_Celltype_",norMeth_bulk,".xlsx"),
#           rowNames=TRUE)
#print(paste0("data/RMSE_noPro/SumTable_BulkRNA_RMSE_AllDataset_Celltype_",norMeth_bulk,".xlsx"))
write.xlsx(SumTable_Celltype, file = paste0("~/SumTable_BulkRNA_Pearson_AllDataset_Celltype_",norMeth_bulk,".xlsx"),
           rowNames=TRUE)
print(paste0("~/SumTable_BulkRNA_Pearson_AllDataset_Celltype_",norMeth_bulk,".xlsx"))

#百分比圖-----------------------------------------------------------------------------------------------------------------------------
library(readr)
Best_mix_noPro <- read_csv("data/Deconvolution_result_TCGA/DeconRes_ref_Best_mix_CGGA_TMM_CibersortEPIC.csv") %>% as.data.frame()
row.names(Best_mix_noPro) <- Best_mix_noPro$...1
Best_mix_noPro$...1 <- NULL
Best_mix_noPro <- Best_mix_noPro %>% t() %>% as.data.frame()
#row.names(Best_mix_noPro) <- str_replace(row.names(Best_mix_noPro), "R_","W")
row.names(Best_mix_noPro) <- str_replace(row.names(Best_mix_noPro), "R","W")
row.names(Best_mix_noPro)[which(row.names(Best_mix_noPro)=="W1")] <- "W01"

CombineTable <- read_csv("data/Deconvolution_result_TCGA/CGGA_Survival_CellType_153_TMM_CombineTable.csv") %>% as.data.frame()
#row.names(CombineTable) <- CombineTable$...1
row.names(CombineTable) <- CombineTable$PATIENT_ID
CombineTable$...1 <- NULL
CombineTable$PATIENT_ID <- NULL

tempdf <- CombineTable[,c(1:10)] %>% as.data.frame()
tempdf$Patient_ID <- row.names(tempdf)
for (i in 1:(ncol(tempdf)-1)) {
  if(i==1){
    temptemp <- tempdf[,c(i,ncol(tempdf))] %>% as.data.frame()
    colnames(temptemp) <- c("Perc","Patient_ID")
    Celltype <- rep(colnames(tempdf)[i],nrow(tempdf))
  }else{
    a <- tempdf[,c(i,ncol(tempdf))] %>% as.data.frame()
    colnames(a) <- c("Perc","Patient_ID")
    temptemp <- rbind(temptemp,a) %>% as.data.frame()
    Celltype <- append(Celltype,rep(colnames(tempdf)[i],nrow(tempdf)))
    remove(a)
  }
}
temptemp$Celltype <- Celltype
remove(tempdf, Celltype,i)

factorchosen <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells","Oligodendrocytes",
                  "T cells","Tumor cells","B cells","Mural cells")

p1 <- ggplot(temptemp, aes(x = Patient_ID, y = Perc, fill = Celltype)) +
  geom_bar(position = "fill", stat = "identity", width = 0.8) +
  labs(x = "", y = "Percentage", fill = "Cell Types") +
  scale_fill_manual(values=c("Dendritic cells" = "#A6CEE3",
                             "Endothelial cells" = "#1F78B4",
                             "Macrophage-like GAMs" = "#B2DF8A",
                             "Microglia-like GAMs" = "#33A02C",
                             "NKT-like cells" = "#FB9A99",
                             "Oligodendrocytes" = "#E31A1C",
                             "T cells" = "#FDBF6F",
                             "Tumor cells" = "#FF7F00",
                             "B cells" = "#CAB2D6",
                             "Mural cells" = "#6A3D9A")) +
  # "Proliferating GAMs/tumor cells" = "#FFD92F")) +
  theme_minimal() +
  guides(fill=guide_legend(nrow=2,byrow=FALSE)) +
  theme(legend.position = "top", #bottom 底部
        axis.line.y = element_line(color = "black", size = 1),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black", size = 16),
        axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1.1, vjust = 0.5, color = "black"), 
        text = element_text(size = 14),
        legend.text = element_text(color = "black", size = 15),
        legend.title = element_text(color = "black", size = 16),
        legend.margin = margin(),
        legend.spacing.x = unit(0.25, "cm"),
        legend.key.height = unit(0.25, "cm"),legend.key.width = unit(0.4, "cm")) 

p1

#存活--------------------------------------------------------------------------
library(readr)
library(magrittr)
library(stringr)
library(tidyverse)
library(survival)
library(survminer)

#TCGA
DeconRes <- read_csv("data/Deconvolution_result_TCGA/DeconRes_ref_Best_mix_TCGA_TPM_CibersortEPIC_noPro.csv") %>% as.data.frame()
row.names(DeconRes) <- DeconRes$...1
DeconRes$...1 <- NULL
DeconRes_newly <- DeconRes[,grep("TCGA-\\w{2}-\\w{4}-01", colnames(DeconRes))] %>% as.data.frame()
#去除duplication
a <- which(duplicated(str_sub(colnames(DeconRes_newly), start = 1, end = 12))==TRUE)
DeconRes_newly <- DeconRes_newly[,-a] %>% as.data.frame()
colnames(DeconRes_newly) <- str_sub(colnames(DeconRes_newly), start = 1, end = 12)

treatment <- read_delim("data/Deconvolution_result_TCGA/data_timeline_treatment.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)
status <- read_delim("data/Deconvolution_result_TCGA/data_timeline_status.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)
clinical_patient <- read_delim("data/Deconvolution_result_TCGA/data_clinical_patient.txt", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE) %>% as.data.frame()
colnames(clinical_patient) <- clinical_patient[4,]
clinical_patient <- clinical_patient[-c(1:4),] %>% as.data.frame()
Patient_ID <- colnames(DeconRes_newly)
factorname <- c("PATIENT_ID","SUBTYPE","AGE","SEX","DAYS_LAST_FOLLOWUP","ETHNICITY","HISTORY_NEOADJUVANT_TRTYN","RACE","RADIATION_THERAPY","OS_STATUS",
                "OS_MONTHS","PFS_STATUS","PFS_MONTHS","GENETIC_ANCESTRY_LABEL")     
clinical_patient <- clinical_patient[which(clinical_patient$PATIENT_ID %in% Patient_ID),factorname] %>% as.data.frame()

DeconRes_newly <- DeconRes_newly[,clinical_patient$PATIENT_ID] %>% as.data.frame()
#資料缺失
library(TCGAbiolinks)
clinicoft <- c("submitter_id" ,"gender", "age_at_index", "vital_status", "days_to_death", "days_to_last_follow_up")
clin <- GDCquery_clinic (project = "TCGA-GBM", type = "clinical", save.csv = FALSE)
clin_target <- clin[clin$submitter_id %in% colnames(DeconRes_newly), clinicoft] %>% as.data.frame()
remove(clinicoft,clin)
clinical <- merge(clinical_patient,clin_target, by.x = "PATIENT_ID", by.y = "submitter_id", all.x = TRUE) %>% as.data.frame()
clinical <- clinical[,-which(colnames(clinical) %in% c("SEX","AGE"))] %>% as.data.frame()
remove(clin_target,clinical_patient,DeconRes,a,b,factorname,Patient_ID,sample_ID)

clinical$OS_STATUS <- str_sub(clinical$OS_STATUS, start = 1, end = 1) %>% as.numeric() 
clinical$PFS_STATUS <- str_sub(clinical$PFS_STATUS, start = 1, end = 1) %>% as.numeric() 
clinical$OS_MONTHS <- clinical$OS_MONTHS %>% as.double() 
clinical$PFS_MONTHS <- clinical$PFS_MONTHS %>% as.double() 

DeconRes_newly <- DeconRes_newly %>% t() %>% as.data.frame()
combinetable <- cbind(DeconRes_newly[order(row.names(DeconRes_newly)),],clinical[order(clinical$PATIENT_ID),]) %>% as.data.frame()

#CGMH
library(readr)
library(readxl)
#CGMH_TMM <- read_csv("data/Deconvolution_result_TCGA/DeconRes_ref_Best_mix_CGMH_TMM_CibersortEPIC.csv") %>% as.data.frame()
CGMH_TMM <- read_csv("data/Deconvolution_result_TCGA/DeconRes_ref_Best_mix_CGMH_polyA_TMM_CibersortEPIC.csv") %>% as.data.frame()
row.names(CGMH_TMM) <- CGMH_TMM$...1
CGMH_TMM$...1 <- NULL
#colnames(CGMH_TMM) <- str_replace(colnames(CGMH_TMM),"R_","W")  #exome
colnames(CGMH_TMM) <- str_replace(colnames(CGMH_TMM),"R","W")  #polyA
CGMH_TMM <- CGMH_TMM %>% t() %>% as.data.frame()
row.names(CGMH_TMM)[which(row.names(CGMH_TMM)=="W1")] <- "W01"
CGMH_TMM$PATIENT_ID <- row.names(CGMH_TMM)

CGMH_clinical <- read_excel("data/A3B_clinical_20220505.xlsx") %>% as.data.frame()
CGMH_clinical$W_sample <- str_replace(CGMH_clinical$W_sample, "w", "W")
CGMH_clinical <- CGMH_clinical[which(CGMH_clinical$grade==4),] %>% as.data.frame()
CGMH <- merge(CGMH_TMM,CGMH_clinical[,c("W_sample","PFS_status","PFS_days","OS_status","OS_days")], by.x = "PATIENT_ID", by.y = "W_sample") %>% as.data.frame()
CGMH$OS_MONTHS <- CGMH$OS_days/30.42
colnames(CGMH)[which(colnames(CGMH)=="OS_status")] <- "OS_STATUS"
CGMH$PFS_MONTHS <- CGMH$PFS_days/30.42
colnames(CGMH)[which(colnames(CGMH)=="PFS_status")] <- "PFS_STATUS"
write.csv(CGMH, file = "data/Deconvolution_result_TCGA/CGMH_Survival_CellType_polyA12GBM_TMM_CombineTable.csv")

CGMH_polyA <- CGMH
CGMH_exome$PATIENT_ID <- paste0(CGMH_exome$PATIENT_ID,"_ExomeCapture")
CGMH_polyA$PATIENT_ID <- paste0(CGMH_polyA$PATIENT_ID,"_PolyA")
target <- colnames(CGMH_exome)
CGMH_polyA <- CGMH_polyA[,match(target, colnames(CGMH_polyA))] %>% as.data.frame()
CGMH_combine <- rbind(CGMH_exome,CGMH_polyA) %>% as.data.frame()
write.csv(CGMH_combine, file = "data/Deconvolution_result_TCGA/CGMH_CellType_exome41polyA12GBM_TMM_Table.csv")

#CGGA
CGGA_TMM <- read_csv("data/Deconvolution_result_TCGA/DeconRes_ref_Best_mix_CGGA_TMM_CibersortEPIC.csv") %>% as.data.frame()
row.names(CGGA_TMM) <- CGGA_TMM$...1
CGGA_TMM$...1 <- NULL
CGGA_TMM <- CGGA_TMM %>% t() %>% as.data.frame()
CGGA_TMM$PATIENT_ID <- row.names(CGGA_TMM)

CGGA_clinical <- read_excel("data/1018_CGGA_RNASeq_Sample_Info.xlsx") %>% as.data.frame()
CGGA_clinical <- CGGA_clinical[which(CGGA_clinical$Histology=="GBM"),] %>% as.data.frame()
CGGA_clinical <- CGGA_clinical[which(CGGA_clinical$IDH_mutation_status=="Wildtype"),] %>% as.data.frame()
CGGA_clinical <- CGGA_clinical[which(CGGA_clinical$`1p19q_codeletion_status`=="Non-codel"),] %>% as.data.frame()
CGGA <- merge(CGGA_TMM,CGGA_clinical[,c("CGGA_ID","OS","Censor (alive=0; dead=1)")], by.x = "PATIENT_ID", by.y = "CGGA_ID") %>% as.data.frame()
colnames(CGGA)[13:14] <- c("OS_days","OS_status")
CGGA$OS_days <- CGGA$OS_days %>% as.numeric()
CGGA$OS_status <- CGGA$OS_status %>% as.numeric()
CGGA$OS_MONTHS <- CGGA$OS_days/30.42
colnames(CGGA)[which(colnames(CGGA)=="OS_status")] <- "OS_STATUS"
write.csv(CGGA, file = "data/Deconvolution_result_TCGA/CGGA_Survival_CellType_153_TMM_CombineTable.csv")

#best p value
oridata <- CGMH
low_25 <- ceiling(length(oridata$PATIENT_ID)*0.25)
high_25 <- floor(length(oridata$PATIENT_ID)*0.75)

oridata <- oridata[order(oridata$`Tumor cells`),] %>% as.data.frame()
pval_OS <- c()
median_H_OS <- c()
median_L_OS <- c()
pval_PFS <- c()
median_H_PFS <- c()
median_L_PFS <- c()
for (i in low_25:high_25) {
  oridata$S_group <- "H"
  oridata$S_group[1:i] <- "L"
  fit_OS <- surv_fit(Surv(OS_MONTHS,OS_STATUS) ~ S_group, data = oridata)
  pval_OS <- append(pval_OS, round(surv_pvalue(fit = fit_OS)[1,2],digits = 5))
  median_H_OS <- append(median_H_OS, surv_median(fit_OS)[1,2])
  median_L_OS <- append(median_L_OS, surv_median(fit_OS)[2,2])
  fit_PFS <- surv_fit(Surv(PFS_MONTHS,PFS_STATUS) ~ S_group, data = oridata)
  pval_PFS <- append(pval_PFS, round(surv_pvalue(fit = fit_PFS)[1,2],digits = 5))
  median_H_PFS <- append(median_H_PFS, surv_median(fit_PFS)[1,2])
  median_L_PFS <- append(median_L_PFS, surv_median(fit_PFS)[2,2])
}
sumtable <- data.frame(Number=c(low_25:high_25),
                       Pvalue_OS=pval_OS,
                       Median_surv_H_OS=median_H_OS,
                       Median_surv_L_OS=median_L_OS,
                       Pvalue_PFS=pval_PFS,
                       Median_surv_H_PFS=median_H_PFS,
                       Median_surv_L_PFS=median_L_PFS)
sumtable_Tumor <- sumtable
write.csv(sumtable_Mural, file = "data/Deconvolution_result_TCGA/Survival_result_CGGA153_MuralCells.csv")
write.csv(combinetable, file = "data/Deconvolution_result_TCGA/TCGA_Survival_CellType_153_CombineTable.csv")

#survival graph
oridata <- oridata[order(oridata$`Microglia-like GAMs`),] %>% as.data.frame()
i <- 63
oridata$S_group <- "H"
oridata$S_group[1:i] <- "L"
fit_OS <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ S_group, data = oridata)
fit_PFS <- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ S_group, data = oridata)
#break_n <- c(1000)
S_fig <- list()
pvalue_OS <- surv_pvalue(fit = fit_OS, data = oridata)
S_fig[[5]] <- ggsurvplot(fit_OS,
                         pval = paste0("p=",format(pvalue_OS$pval, digits = 4)),
                         pval.coord = c(60,0.8),
                         #risk.table = TRUE,
                         risk.table = FALSE,
                         legend.title = "Microglia-like GAMs",
                         surv.median.line = c("hv"),
                         legend = c("top"), 
                         #break.x.by = c(1000),
                         ylab="OS Probability",
                         xlab="Months after first diagnosis",
                         palette= c("#990000","#0000FF"),
                         legend.labs = c("High","Low"),
                         font.title = 16,
                         font.x =  16,
                         font.y = 16,
                         font.tickslab = 16,
                         tables.y.text = FALSE)

pvalue_PFS <- surv_pvalue(fit = fit_PFS, data = oridata)
S_fig[[6]] <- ggsurvplot(fit_PFS,
                         pval = paste0("p=",format(pvalue_PFS$pval, digits = 4)),
                         pval.coord = c(30,0.8),
                         #risk.table = TRUE,
                         risk.table = FALSE,
                         legend.title = "Oligodendrocytes",
                         surv.median.line = c("hv"),
                         legend = c("top"), 
                         #break.x.by = c(1000),
                         ylab="PFS Probability",
                         xlab="Months after first diagnosis",
                         palette= c("#990000","#0000FF"),
                         legend.labs = c("High","Low"),
                         font.title = 16,
                         font.x =  16,
                         font.y = 16,
                         font.tickslab = 16,
                         tables.y.text = FALSE)

arrange_ggsurvplots(S_fig, print = TRUE,
                    ncol = 2, nrow = 3, risk.table.height = 0.3, surv.plot.height = 0.7)

#survival graph--2 feature
factorchosen <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells","Oligodendrocytes",
                  "T cells","Tumor cells","B cells","Mural cells")
#value OS/PFS: "Dendritic cells" 50/105;"Endothelial cells" 45/86;"Macrophage-like GAMs" 105/65;"Microglia-like GAMs" 109/63;
#"Oligodendrocytes" 104/54;"Tumor cells" 100/100;"B cells" 112/39;"Mural cells" 106/75

factor1 <- "Macrophage-like GAMs"
factor2 <- "Microglia-like GAMs"
i <- 44
#j <- 51
k <- 86
#l <- 40
#oridata_temp <- oridata[,which(colnames(oridata) %in% c(factor1,factor2,"OS_MONTHS","OS_STATUS","PFS_MONTHS","PFS_STATUS"))] %>% as.data.frame()
#oridata_temp <- oridata[,which(colnames(oridata) %in% c(factor1,factor2,"OS_days","OS_status","PFS_days","PFS_status"))] %>% as.data.frame()
oridata_temp <- oridata[,which(colnames(oridata) %in% c(factor1,factor2,"OS_days","OS_status"))] %>% as.data.frame()
colnames(oridata_temp)[1:2] <- c("Celltype1","Celltype2")
oridata_temp <- oridata_temp[order(oridata_temp$Celltype1),] %>% as.data.frame()
oridata_temp$Celltype1_group_OS <- "H"
oridata_temp$Celltype1_group_OS[1:i] <- "L"
#oridata_temp$Celltype1_group_PFS <- "H"
#oridata_temp$Celltype1_group_PFS[1:j] <- "L"
oridata_temp <- oridata_temp[order(oridata_temp$Celltype2),] %>% as.data.frame()
oridata_temp$Celltype2_group_OS <- "H"
oridata_temp$Celltype2_group_OS[1:k] <- "L"
#oridata_temp$Celltype2_group_PFS <- "H"
#oridata_temp$Celltype2_group_PFS[1:l] <- "L"
oridata_temp$D_group_OS <- paste0(oridata_temp$Celltype1_group_OS,oridata_temp$Celltype2_group_OS)
#oridata_temp$D_group_PFS <- paste0(oridata_temp$Celltype1_group_PFS,oridata_temp$Celltype2_group_PFS)
#fit_OS <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ D_group_OS, data = oridata_temp)
#fit_PFS <- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ D_group_PFS, data = oridata_temp)
fit_OS <- survfit(Surv(OS_days, OS_status) ~ D_group_OS, data = oridata_temp)
#fit_PFS <- survfit(Surv(PFS_days, PFS_status) ~ D_group_PFS, data = oridata_temp)
#break_n <- c(1000)
S_fig_OS <- list()
#S_fig_PFS <- list()
S_fig_OS[[3]] <- ggsurvplot(fit_OS,
                            pval = FALSE,
                            risk.table = TRUE,
                            legend.title = "",
                            legend = c("top"), 
                            ylab="OS Probability",
                            #xlab="Months after first diagnosis",
                            xlab="Days after first diagnosis",
                            palette= c("#990000","#FF9999","#0000FF","#6699FF"),
                            legend.labs = c("High/High","High/Low","Low/High","Low/Low"),
                            font.title = 16,
                            font.x =  16,
                            font.y = 16,
                            font.tickslab = 16,
                            tables.y.text = FALSE)

S_fig_PFS[[3]] <- ggsurvplot(fit_PFS,
                             pval = FALSE,
                             risk.table = TRUE,
                             legend.title = "",
                             legend = c("top"), 
                             ylab="PFS Probability",
                             #xlab="Months after first diagnosis",
                             xlab="Days after first diagnosis",
                             palette= c("#990000","#FF9999","#0000FF","#6699FF"),
                             legend.labs = c("High/High","High/Low","Low/High","Low/Low"),
                             font.title = 16,
                             font.x =  16,
                             font.y = 16,
                             font.tickslab = 16,
                             tables.y.text = FALSE)

arrange_ggsurvplots(S_fig, print = TRUE,
                    ncol = 3, nrow = 2, risk.table.height = 0.3, surv.plot.height = 0.7)

#一個細胞種類四個分群
library("RColorBrewer")
oridata <- combinetable
low_25 <- ceiling(length(oridata$PATIENT_ID)*0.25)
high_25 <- floor(length(oridata$PATIENT_ID)*0.75)
median_50 <- floor(length(oridata$PATIENT_ID)*0.5)

oridata <- oridata[order(oridata$`Macrophage-like GAMs`),] %>% as.data.frame()
oridata$S_group <- "Q1"
oridata$S_group[1:(-1+low_25)] <- "Q4"
oridata$S_group[low_25:median_50] <- "Q3"
oridata$S_group[1+median_50:high_25] <- "Q2"
fit_OS <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ S_group, data = oridata)
fit_PFS <- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ S_group, data = oridata)
S_fig <- list()
S_fig[[1]] <- ggsurvplot(fit_OS,
                         risk.table = TRUE,
                         legend.title = "Mural cells",
                         legend = c("right"), 
                         ylab="OS Probability",
                         xlab="Months after first diagnosis",
                         palette= brewer.pal(n = 4, name = "RdYlBu"),
                         # legend.labs = c("High","Low"),
                         font.title = 16,
                         font.x =  16,
                         font.y = 16,
                         font.tickslab = 16,
                         tables.y.text = FALSE)

S_fig[[2]] <- ggsurvplot(fit_PFS,
                         risk.table = TRUE,
                         legend.title = "Mural cells",
                         legend = c("right"), 
                         ylab="PFS Probability",
                         xlab="Months after first diagnosis",
                         palette= brewer.pal(n = 4, name = "RdYlBu"),
                         #    legend.labs = c("High","Low"),
                         font.title = 16,
                         font.x =  16,
                         font.y = 16,
                         font.tickslab = 16,
                         tables.y.text = FALSE)

arrange_ggsurvplots(S_fig, print = TRUE,
                    ncol = 2, nrow = 1, risk.table.height = 0.3, surv.plot.height = 0.7)

#polyA和exomecapture的比較--------------------------------------------------------------------------------------------------------------
CombineTable <- read_csv("data/Deconvolution_result_TCGA/CGMH_CellType_exome70polyA20_TMM_Table.csv") %>% as.data.frame()
row.names(CombineTable) <- CombineTable$PATIENT_ID
CombineTable$...1 <- NULL
CombineTable$PATIENT_ID <- str_split_fixed(row.names(CombineTable),pattern = "_",2)[,1]

DB <- CombineTable[which(CombineTable$Patient_ID %in% CombineTable$Patient_ID[which(duplicated(CombineTable$Patient_ID)==TRUE)]),] %>% as.data.frame()
DB$methods <- str_split_fixed(row.names(DB),pattern = "_",2)[,2]

library(corrplot)
library(ggrepel)
library(ggpubr)
library("gridExtra")

EC <- DB[which(DB$methods=="ExomeCapture"),] %>% as.data.frame()
row.names(EC) <- EC$Patient_ID
colnames(EC)[1:11] <- paste0(colnames(EC)[1:11],"_ExomeCapture")
EC$methods <- NULL
PolyA <- DB[which(DB$methods=="PolyA"),] %>% as.data.frame()
row.names(PolyA) <- PolyA$Patient_ID
colnames(PolyA)[1:11] <- paste0(colnames(PolyA)[1:11],"_PolyA")
PolyA$methods <- NULL
DB <- merge(EC,PolyA, by.x = "Patient_ID",by.y = "Patient_ID") %>% as.data.frame()
DB <- DB[,-grep(pattern = "otherCells", colnames(DB))] %>% as.data.frame()
row.names(DB) <- DB$Patient_ID
DB$Patient_ID <- NULL

M <- cor(DB)
testRes <- cor.mtest(DB, conf.level = 0.95)
corrplot(M,
         p.mat = testRes$p,
         #         insig='blank',
         sig.level = 0.05,
         addrect = 10,
         rect.col = 'red4',
         rect.lwd = 3,
         pch.col = 'grey50',
         order = 'hclust',
         #         type = 'lower',
         tl.col = 'black',
         addCoef.col = 'white',
         number.cex = 0.6,
         #         tl.srt = 45,
         #         diag=FALSE,
         col = COL2('PuOr', 200))

factorName <- colnames(DB)
S_fig <- list()
factorName[2]
factorName[12]
S_fig[[10]] <- ggplot(data = DB, aes(x=DB$`Mural cells_ExomeCapture`, y=DB$`Mural cells_PolyA`, label = rownames(DB))) +
  geom_point() +
  geom_text_repel() +
  geom_smooth(method=lm, se=FALSE) +
  labs(y="Mural cells - PolyA", x = "Mural cells - ExomeCapture") +
  theme_bw()

grid.arrange(S_fig[[1]],S_fig[[2]],S_fig[[3]],S_fig[[4]],S_fig[[5]],
             S_fig[[6]],S_fig[[7]],S_fig[[8]],S_fig[[9]],S_fig[[10]],
             nrow = 3)                       # Number of rows


#SCTransformation_new
install.packages("remotes")
remotes::install_github("omnideconv/immunedeconv")

library(devtools)
install_github("immunogenomics/harmony", force = TRUE)
#restart

library(ggplot2)
library(magrittr)
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
library(tidyverse)
library(sctransform)

count <- Matrix::readMM("data2/GSE182109/matrix.mtx") %>% as.matrix() %>% as.data.frame()
barcodes <- read_csv("data2/GSE182109/barcodes.tsv", 
                     col_names = FALSE) %>% as.data.frame()
genes <- read_delim("data2/GSE182109/genes.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    col_names = FALSE, trim_ws = TRUE) %>% as.data.frame()
GBM.data <- cbind(genes,count)
row.names(GBM.data) <-  GBM.data$X1
GBM.data <- GBM.data[,-c(1,2)]
colnames(GBM.data) <- barcodes$X1

Final_Metadata <- read_csv("data/Meta_GBM.txt") %>% as.data.frame()
Final_Metadata <- Final_Metadata[-1,] %>% as.data.frame()

seurat_ob <- CreateSeuratObject(counts = GBM.data, project = "gbm", min.cells = 3, min.features = 200)
seurat_ob <- PercentageFeatureSet(seurat_ob, pattern = "^MT-", col.name = "percent.mt")

target <- colnames(seurat_ob)
Final_Metadata <- Final_Metadata[match(target, Final_Metadata$NAME),] %>% as.data.frame()
Patient <- factor(Final_Metadata$Patient)
names(Patient) <- colnames(seurat_ob)
Type <- factor(Final_Metadata$Type)
names(Type) <- colnames(seurat_ob)
Assignment <- factor(Final_Metadata$Assignment)
names(Assignment) <- colnames(seurat_ob)
seurat_ob <- AddMetaData(
  object = seurat_ob,
  metadata = list(Patient,Type,Assignment),
  col.name = c("Patient","Type","Assignment")
)
seurat_ob <- subset(seurat_ob, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 20)
remove(barcodes,count,genes,Assignment,Patient,Type)

seurat_ob <- subset(x = seurat_ob, subset = Patient %in% c("ndGBM-01","ndGBM-02","ndGBM-03","ndGBM-04","ndGBM-05","ndGBM-06","ndGBM-07","ndGBM-08","ndGBM-09",
                                                           "rGBM-01","rGBM-02","rGBM-03","rGBM-04","rGBM-05"))
VlnPlot(seurat_ob, features = c("nFeature_RNA"), ncol = 1, repel = TRUE, raster=FALSE)

seurat_ob_SCT <- SCTransform(seurat_ob,
                             vars.to.regress = "percent.mt",
                             ncells = 10000,
                             variable.features.n = 5000,
                             return.only.var.genes = FALSE,
                             seed.use = 42)
seurat_ob_SCT <- RunPCA(seurat_ob_SCT,npcs = 50) #npcs = 50
seurat_ob_SCT@meta.data$orig.ident <- as.character(seurat_ob_SCT@meta.data$orig.ident)
seurat_ob_SCT@meta.data$orig.ident <- as.factor(seurat_ob_SCT@meta.data$orig.ident)
seurat_ob_SCT_hmy <- RunHarmony(seurat_ob_SCT,
                                assay.use = "SCT",
                                reduction.use = "pca",
                                dims.use = 1:50,
                                group.by.vars = "orig.ident")
seurat_ob_SCT_hmy <- RunUMAP(seurat_ob_SCT_hmy,
                             assay = "SCT",
                             reduction = "harmony",
                             dims = 1:50,
                             n.neighbors = 20, min.dist = 0.15, seed.use = 42)
seurat_ob_SCT_hmy <- FindNeighbors(seurat_ob_SCT_hmy, reduction = "harmony", dims = 1:50)
seurat_ob_SCT_hmy <- FindClusters(seurat_ob_SCT_hmy, resolution = 1.0)
DimPlot(seurat_ob_SCT_hmy, label = TRUE, repel = TRUE, group.by = "seurat_clusters", raster=FALSE)
FeaturePlot(seurat_ob_SCT_hmy, features = c("MS4A1","IGFBP2"),raster=FALSE)

#The ‘corrected’ UMI counts are stored in pbmc[["SCT"]]@counts.
#We store log-normalized versions of these corrected counts in pbmc[["SCT"]]@data, which are very helpful for visualization

GSE182109_metadata <- read_csv("data/GSE182109_ndGBM_rGBM_Nond1011_metadata.txt") %>% as.data.frame()
target <- colnames(seurat_ob_SCT_hmy)
GSE182109_metadata <- GSE182109_metadata[match(target, GSE182109_metadata$cells),] %>% as.data.frame()
CellType1 <- factor(GSE182109_metadata$CellType1)
names(CellType1) <- colnames(seurat_ob_SCT_hmy)
seurat_ob_SCT_hmy <- AddMetaData(
  object = seurat_ob_SCT_hmy,
  metadata = list(CellType1),
  col.name = c("CellType1")
)
remove(target,CellType1)
DimPlot(seurat_ob_SCT_hmy, label = TRUE, repel = TRUE, group.by = "CellType1", raster=FALSE)

cc.

a <- read_csv("data/FeatureGeneSet/Cell_sampling_11celltype_GSE182109_DatasetWhole_MinDist.csv") %>% as.data.frame()
seurat_ob_SCT_hmy@meta.data$cell_ID <- colnames(seurat_ob_SCT_hmy)
sob <- subset(x = seurat_ob_SCT_hmy, subset = cell_ID %in% a$cell_ID)

Meth <- "SCT" #RawCount, LogNormalize, SCT, TMM, TPM
Sampling <- "DatasetWhole_MinDist"  #ALL,DatasetWhole_MinDist

#ALL=seurat_ob_SCT_hmy, DatasetWhole_MinDist=sob
if(Meth=="SCT" & Sampling=="ALL"){
  ex_matrix <- seurat_ob_SCT_hmy[["SCT"]]@scale.data %>% as.data.frame()
  print("SCT ALL")
}else if(Meth=="RawCount" & Sampling=="ALL"){
  ex_matrix <- seurat_ob_SCT_hmy[["SCT"]]@counts %>% as.data.frame()
  print("RawCount ALL")
}else if(Meth=="LogNormalize" & Sampling=="ALL"){
  ex_matrix <- seurat_ob_SCT_hmy[["SCT"]]@data %>% as.data.frame()
  print("LogNormalize ALL")
}else if(Meth=="TMM" & Sampling=="ALL"){
  TMM <- read_csv("data/GSE182109_norm_tmm_df.csv") %>% as.data.frame()
  row.names(TMM) <- TMM$...1
  TMM$...1 <- NULL
  ex_matrix <- TMM[,which(colnames(TMM) %in% colnames(seurat_ob_SCT_hmy))] %>% as.data.frame() #TMM for ALL
  target <- colnames(seurat_ob_SCT_hmy)
  ex_matrix <- ex_matrix[,match(target, colnames(ex_matrix))] %>% as.data.frame()
  print("TMM ALL")
}else if(Meth=="TPM" & Sampling=="ALL"){
  TPM <- read_csv("data/GSE182109_TPM_nor.csv") %>% as.data.frame()
  row.names(TPM) <- TPM$...1
  TPM$...1 <- NULL
  ex_matrix <- TPM[,which(colnames(TPM) %in% colnames(seurat_ob_SCT_hmy))] %>% as.data.frame() #TPM for ALL
  target <- colnames(seurat_ob_SCT_hmy)
  ex_matrix <- ex_matrix[,match(target, colnames(ex_matrix))] %>% as.data.frame()
  print("TPM ALL")
}else if(Meth=="SCT" & Sampling=="DatasetWhole_MinDist"){
  ex_matrix <- sob[["SCT"]]@scale.data %>% as.data.frame()
  print("SCT DatasetWhole_MinDist")
}else if(Meth=="RawCount" & Sampling=="DatasetWhole_MinDist"){
  ex_matrix <- sob[["SCT"]]@counts %>% as.data.frame()
  print("RawCount DatasetWhole_MinDist")
}else if(Meth=="LogNormalize" & Sampling=="DatasetWhole_MinDist"){
  ex_matrix <- sob[["SCT"]]@data %>% as.data.frame()
  print("LogNormalize DatasetWhole_MinDist")
}else if(Meth=="TMM" & Sampling=="DatasetWhole_MinDist"){
  TMM <- read_csv("data/GSE182109_norm_tmm_df.csv") %>% as.data.frame()
  row.names(TMM) <- TMM$...1
  TMM$...1 <- NULL
  ex_matrix <- TMM[,which(colnames(TMM) %in% colnames(sob))] %>% as.data.frame() #TMM for ALL
  target <- colnames(sob)
  ex_matrix <- ex_matrix[,match(target, colnames(ex_matrix))] %>% as.data.frame()
  print("TMM DatasetWhole_MinDist")
}else if(Meth=="TPM" & Sampling=="DatasetWhole_MinDist"){
  TPM <- read_csv("data/GSE182109_TPM_nor.csv") %>% as.data.frame()
  row.names(TPM) <- TPM$...1
  TPM$...1 <- NULL
  ex_matrix <- TPM[,which(colnames(TPM) %in% colnames(seurat_ob_SCT_hmy))] %>% as.data.frame() #TPM for ALL
  target <- colnames(seurat_ob_SCT_hmy)
  ex_matrix <- ex_matrix[,match(target, colnames(ex_matrix))] %>% as.data.frame()
  print("TPM DatasetWhole_MinDist")
}else{
  print("Something wrong!")
}

FAM <- "611"  #211,411,611
for (gene_n in c(20,50,100)) {
  FM <- read_csv(paste0("data/FeatureGeneSet/Findallmarker_GSE182109_11celltype_",FAM,"_",Sampling,"_n",gene_n,".csv")) %>% as.data.frame()
  FM <- FM[-which(FM$cluster=="Proliferating tumor cells"),] %>% as.data.frame()
  target_gene <- FM$gene
  df <- ex_matrix[which(row.names(ex_matrix) %in% target_gene),] %>% t() %>% as.data.frame()
  if(Sampling=="ALL"){
    df$celltype <- seurat_ob_SCT_hmy@meta.data$CellType1
  }else{
    df$celltype <- sob@meta.data$CellType1
  }
  df2 <- df %>% group_by(celltype) %>% summarise_all(mean) %>% t() %>% as.data.frame()
  colnames(df2) <- df2[1,]
  df2 <- df2[-1,]
  colnames(df2)[grep(pattern = "Proliferating",colnames(df2))] <- "noPro"
  df2$noPro <- NULL
  write.csv(df2, file = paste0("data/Reference_noPro/Reference_GSE182109_11celltype_",FAM,"_",Sampling,"_n",gene_n,"_",Meth,".csv"))
  remove(gene_n,FM,target_gene,df,df2)
}
#-------------------------------------------------------------------------------------------
#Cell number
seurat_ob_SCT_hmy@meta.data$cell_name <- colnames(seurat_ob_SCT_hmy)
cell_num <- c(100,1000,2500,6000)

#Tumor cell
dataset <- "GSE182109"     #NC, GSE131928, GSE182109
cell_num <- 6000  #6000,2500,1000,100
perc <- seq(from = 60, to = 100, by = 2)

RaTio <- c(0.2,0.4,0.6,0.8)
TumorCell_ID <- seurat_ob_SCT_hmy@meta.data$cell_name[seurat_ob_SCT_hmy@meta.data$CellType=="Tumor cells"]
GAMCell_ID <- seurat_ob_SCT_hmy@meta.data$cell_name[seurat_ob_SCT_hmy@meta.data$CellType %in% c("Microglia-like GAMs", "Macrophage-like GAMs")]
OtherCell_ID <- seurat_ob_SCT_hmy@meta.data$cell_name[seurat_ob_SCT_hmy@meta.data$CellType %in% c("NKT-like cells","T cells","B cells","Oligodendrocytes",
                                                                                                  "Dendritic cells","Endothelial cells","Mural cells")]

a <- data.frame(cell_ID = seq(from = 1, to = cell_num, by = 1))

if(dataset=="GSE182109"){
  Nn <- 12
}else{
  Nn <- 6
}
for (i in perc) {
  TG_num <- cell_num*i*0.01
  Other_num <- cell_num-TG_num
  print(paste0("START--",i,"% Tumor/GAM=",TG_num," Others=",Other_num))
  for (j in RaTio) {
    print(paste0("START Raio=",j))
    TumorCell_num <- round(TG_num*j, digits=0) %>% as.numeric()
    GAMCell_num <- TG_num-TumorCell_num
    print(GAMCell_num)
    for (k in 1:Nn) {
      print(k)
      Tumor_ID <- sample(TumorCell_ID,TumorCell_num,replace=TRUE)
      GAM_ID <- sample(GAMCell_ID,GAMCell_num,replace=TRUE)
      Others_ID <- sample(OtherCell_ID,Other_num,replace=TRUE)
      Total_cell_ID <- append(Tumor_ID,GAM_ID)
      Total_cell_ID <- append(Total_cell_ID,Others_ID)
      a <- a %>% add_column(Total_cell_ID)
      remove(Tumor_ID,GAM_ID,Others_ID,Total_cell_ID)
    }
    print(paste0("Completed Raio=",j))
    remove(TumorCell_num,GAMCell_num)
  }
  print(paste0("Completed--",i))
  remove(TG_num,Other_num)
}
write.csv(a, file = paste0("data/BulkRNA_noPro/",dataset,"_BulkRNA_SampleID_N",cell_num,".csv"))
#--------------------
scData <- seurat_ob_SCT_hmy

ex_matrix <- scData[["SCT"]]@counts %>% as.data.frame()

#enter here
dataset <- "GSE182109"     #NC, GSE131928, GSE182109
cell_num <- 1000  #6000,2500,1000,100

a <- read_csv(paste0("data/BulkRNA_noPro/",dataset,"_BulkRNA_SampleID_N",cell_num,".csv")) %>% as.data.frame()
a$...1 <- NULL
row.names(a) <- a$cell_ID
a$cell_ID <- NULL
if(dataset=="GSE182109"){
  NN <- 1008
}else{
  NN <- 504
}
colnames(a) <- lapply(c(1:NN), function(x)paste0("GBM",x)) %>% unlist() %>% as.vector()
temp_table <- data.frame(gene_name = row.names(ex_matrix))

for (i in c(1:length(colnames(a)))) {
  temp_table <- temp_table %>% add_column(ex_matrix[,colnames(ex_matrix) %in% a[,i]] %>% rowSums())
}
row.names(temp_table) <- temp_table$gene_name
temp_table$gene_name <- NULL
colnames(temp_table) <- lapply(c(1:NN), function(x)paste0("GBM",x)) %>% unlist() %>% as.vector()
write.csv(temp_table, file = paste0("data/BulkRNA_noPro/",dataset,"_BulkRNA_RawCount_N",cell_num,".csv"))
print(paste0("data/BulkRNA_noPro/",dataset,"_BulkRNA_RawCount_N",cell_num,".csv"))
remove(temp_table,i)
#--------------------
celln <- 100/cell_num

metadata <- data.frame(Cell_ID = scData@meta.data$cell_name,
                       Celltype = scData@meta.data$CellType1) #NOTE NEED CHANGE
temp_table <- data.frame(Celltype = levels(metadata$Celltype))

for (i in c(1:length(colnames(a)))) {
  b <- metadata[metadata$Cell_ID %in% a[,i],]
  temp_table <- temp_table %>% add_column(summary(b$Celltype)*celln)
}
row.names(temp_table) <- temp_table$Celltype
temp_table$Celltype <- NULL
temp_table <- temp_table[-grep(pattern = "Proliferating",row.names(temp_table)),] %>% as.data.frame()

colnames(temp_table) <- lapply(c(1:NN), function(x)paste0("GBM",x)) %>% unlist() %>% as.vector()
write.csv(temp_table, file = paste0("data/BulkRNA_noPro/",dataset,"_BulkRNA_Celltype_N",cell_num,".csv"))
print(paste0("data/BulkRNA_noPro/",dataset,"_BulkRNA_Celltype_N",cell_num,".csv"))
remove(b,metadata,celln,i)


a <- read_csv("data/BulkRNA_20230601/GSE131928_BulkRNA_Celltype_N100.csv") %>% as.data.frame()
row.names(a) <- a$...1
a$...1 <- NULL
colnames(a) <- lapply(c(1:length(colnames(a))), function(x)paste0("cGBM",x)) %>% unlist() %>% as.vector()
write.csv(a, file = "data/BulkRNA_20230601/GSE131928_BulkRNA_Celltype_N100.csv")

#BulkRNA transformation------------------------------------------------------------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("geneplotter")

library(geneplotter)
library(DESeq2)

dataset <- "MergeDataset" #MergeDataset, GSE182109
NN <- 100  #100,1000,2500,6000

CountMatrix <- read_csv(paste0("data/BulkRNA_20230601/",dataset,"_BulkRNA_RawCount_N",NN,".csv")) %>% as.data.frame()
row.names(CountMatrix) <- CountMatrix$geneID
#row.names(CountMatrix) <- CountMatrix$...1
CountMatrix$geneID <- NULL
CountMatrix$...1 <- NULL
coldata <- data.frame(cell_name = colnames(CountMatrix),
                      dataset_name = str_extract(colnames(CountMatrix), pattern = "[:alpha:]+"))

dds <- DESeqDataSetFromMatrix(countData = CountMatrix,
                              colData = coldata,
                              design= ~ 1)
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]
dds <- DESeq(dds, parallel = T)
vsd <- vst(dds, blind=TRUE)
vsd_matrix <- assay(vsd) %>% as.data.frame()
write.csv(vsd_matrix, file = paste0("data/BulkRNA_20230601/",dataset,"_BulkRNA_vst_N",NN,".csv"))

normcount <- counts(dds,normalized=T) %>% as.data.frame()
write.csv(normcount, file = paste0("data/BulkRNA_20230601/",dataset,"_BulkRNA_NormCount_N",NN,".csv"))

ntd <- normTransform(dds)
ntd_matrix <- assay(ntd) %>% as.data.frame()
write.csv(ntd_matrix, file = paste0("data/BulkRNA_20230601/",dataset,"_BulkRNA_ntd_N",NN,".csv"))

remove(dataset,NN,CountMatrix,coldata,dds,vsd,vsd_matrix,ntd,ntd_matrix,normcount)