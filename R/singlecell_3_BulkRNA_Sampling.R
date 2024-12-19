#Cell number
seurat_ob_SCT_hmy@meta.data$cell_name <- colnames(seurat_ob_SCT_hmy)

#Tumor cell
dataset <- "GSE182109"  #NC, GSE131928, GSE182109
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
#extract cell sample--------------------
scData <- seurat_ob_SCT_hmy

ex_matrix <- scData[["SCT"]]@counts %>% as.data.frame()

#enter here
dataset <- "GSE182109"  #NC, GSE131928, GSE182109
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

#three datasets  merge to a big dataset


#calculate each celltype number in each sample for RMSE results============================
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

