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

#"seurat_ob_SCT_hmy" is GSE182109
#two dataset (ALL and DatasetWhole_MinDist)
#input dataset
seurat_ob_SCT_hmy <- LoadH5Seurat("~data/seurat_ob_SCT_hmy_GSE182109_latest.h5Seurat")


#ALL=seurat_ob_SCT_hmy===========
CellType<- c("Proliferating GAMs/tumor cells",
             "Microglia-like GAMs",
             "Macrophage-like GAMs",
             "NKT-like cells",
             "Tumor cells",
             "T cells",
             "B cells",
             "Oligodendrocytes",
             "Dendritic cells",
             "Endothelial cells",
             "Mural cells")
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

ex_matrix <- seurat_ob_SCT_hmy[["SCT"]]@data %>% as.data.frame()
target_celltype <- cluster.markers$cluster %>% unique()

#DatasetWhole_MinDist=sob===========
GSE182109_metadata <- read_csv("~data/GSE182109_ndGBM_rGBM_Nond1011_metadata.txt") %>% as.data.frame()
Celltype <- c("Microglia-like GAMs",
              "Proliferating GAMs/tumor cells",
              "Macrophage-like GAMs",
              "NKT-like cells",
              "Tumor cells",
              "T cells",
              "B cells",
              "Oligodendrocytes",
              "Dendritic cells",
              "Endothelial cells",
              "Mural cells")
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
#write.csv(a, file = "~data/Cell_sampling_11celltype_GSE182109_DatasetWhole_MinDist.csv") #only one in FeatureGeneSet 

seurat_ob_SCT_hmy@meta.data$cell_ID <- colnames(seurat_ob_SCT_hmy)
sob <- subset(x = seurat_ob_SCT_hmy, subset = cell_ID %in% a$cell_ID)


Idents(sob) <- sob@meta.data$CellType1
sob <- PrepSCTFindMarkers(sob)
cluster.markers <- FindAllMarkers(sob,
                                  assay = "SCT",
                                  only.pos = TRUE,
                                  test.use = "wilcox",
                                  recorrect_umi = FALSE,
                                  min.pct = 0.6, #0.2
                                  logfc.threshold = 0.1, #default 0.25
                                  min.diff.pct = 0.1, #0.1
                                  return.thresh = 0.05)
ex_matrix <- sob[["SCT"]]@data %>% as.data.frame() 
target_celltype <- cluster.markers$cluster %>% unique()

#============
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

FAM <- "611"  #211,411,611
Sampling <- "DatasetWhole_MinDist"  #ALL,DatasetWhole_MinDist

df <- cluster.markers %>% as.data.frame()
gene_n <- 20 #20,50,100
df$pct.diff <- df$pct.1 - df$pct.2
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
str_c(df$gene[which(df$cluster=="Proliferating GAMs/tumor cells")] %>% unlist() %>% as.vector(),collapse = ",")
str_c(df$gene[which(df$cluster=="Tumor cells")] %>% unlist() %>% as.vector(),collapse = ",")
str_c(df$gene[which(df$cluster=="Oligodendrocytes")] %>% unlist() %>% as.vector(),collapse = ",")
str_c(df$gene[which(df$cluster=="NKT-like cells")] %>% unlist() %>% as.vector(),collapse = ",")
str_c(df$gene[which(df$cluster=="B cells")] %>% unlist() %>% as.vector(),collapse = ",")
str_c(df$gene[which(df$cluster=="T cells")] %>% unlist() %>% as.vector(),collapse = ",")

write.csv(df, file = paste0("~data/Findallmarker_GSE182109_11celltype_",FAM,"_",Sampling,"_n",gene_n,".csv")) #in FeatureGeneSet