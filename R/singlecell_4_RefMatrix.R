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

