#deconvolution method: EPIC------------------------------------------------------------------------------------------------------------
Dataset <- "MergeDataset"  #MergeDataset, GSE182109
norMeth_bulk <- "NormCount"  #TPM, TMM, vst, NormCount, RawCount, ntd, logTPM, RPKM
cell_num <- 100   #100,1000,2500,6000

exp_matrix <- read_csv("~/data.csv", show_col_types = FALSE) %>% as.data.frame()
row.names(exp_matrix) <- exp_matrix$...1
exp_matrix$...1 <- NULL

#for references
Sampling <- c("ALL")  #ALL, DatasetWhole_MinDist
FAM <- c(611)  #211,411,611
num <- c(100)   #20, 50, 100
norMeth <- c("SCT") #"TPM","TMM","LogNormalize","SCT","RawCount"


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

#deconvolution method: Cibersort----------------------------------------------------------------------------------------------------------
source("~/data/CIBERSORT_ian.R")
decon_Meth <- "Cibersort"   

#for expression
Dataset <- "MergeDataset"  #MergeDataset, GSE182109
norMeth_bulk <- "RPKM"  #TPM, TMM, vst, NormCount, RawCount, ntd, logTPM, rpkm
cell_num <- 100   #100,1000,2500,6000

exp_matrix <- read_csv("~/data.csv", show_col_types = FALSE) %>% as.data.frame()
row.names(exp_matrix) <- exp_matrix$gene_name
exp_matrix$...1 <- NULL
exp_matrix$gene_name <- NULL
exp_matrix[is.na(exp_matrix)] <- 0
exp_matrix <- exp_matrix[!is.infinite(rowSums(exp_matrix)),]

#for references
Sampling <- c("ALL")  #ALL, DatasetWhole_MinDist
FAM <- c(611)  #211,411,611
num <- c(100)   #20, 50, 100
norMeth <- c("SCT") #"TPM","TMM","LogNormalize","SCT","RawCount"

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

#Combine results
decon_Meth <- c("Cibersort","EPIC","ConsesnusTME")
FAM <- c(211,411,611)
Sampling <- c("ALL","DatasetWhole_MinDist")
num <- num <- c(20,50,100)
norMeth <- c("TPM","TMM","LogNormalize","SCT","RawCount")
Celltype <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells",
              "Oligodendrocytes","T cells","Tumor cells","B cells","Mural cells","otherCells")

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
    decon_res <- read_csv(paste0("~data/Deconvolution_noPro/DeconRes_ref_",b,c,d,e,"_BulkRNA_MergeDataset",cell_num,norMeth_bulk,"_",a,".csv"),
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
  combine_res[10,] <- c(Celltype[10],rep(0,length(colnames(combine_res))-1))
  row.names(combine_res) <- combine_res$...1
  combine_res$...1 <- NULL
  re_res <- combine_res
  for (i in 1:length(colnames(combine_res))) {
    combine_res[,i] <- combine_res[,i] %>% as.double()
    re_res[,i] <- NA
    re_res[,i] <- round(combine_res[,i]/sum(combine_res[,i]),digits = 4)*100
  }
  write.csv(re_res, file = paste0("~data/Deconvolution_noPro/DeconRes_ref_Best_mix_",cell_num,"_",norMeth_bulk,"_CibersortEPIC_noPro.csv"))
  print(paste0("~data/Deconvolution_noPro/DeconRes_ref_Best_mix_",cell_num,"_",norMeth_bulk,"_CibersortEPIC_noPro.csv"))
  remove(filelist,combine_res,re_res,i,j)
}
