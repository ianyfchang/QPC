# QPC-GBM v1.1 (last updated 05-30-2024)

library(immunedeconv)
library(tibble)
library(tidyverse)
library(readr)
library(xlsx)

QPCdecon <- function(ex_matrix,ref_list){
  celltype1 <- c("DC","MP","MGTumor","NKT","Oligo","TC","MC")
  df <- list()
  df <- lapply(celltype1, function(x) {
    dec_df <- CIBERSORT(mixture_file = ex_matrix,
                        sig_matrix = ref_list[[x]], 
                        perm = 100,
                        QN = FALSE, 
                        absolute = FALSE,
                        abs_method = "sig.score")
    dec_df_clr <- dec_df[1:10,]
    round(dec_df_clr * 100, digits = 2)
  }
  )
  
  celltype2 <- c("EC","BC")
  df2 <- list()
  df2 <- lapply(celltype2, function(x) {
    dec_df <- deconvolute_epic_custom(gene_expression_matrix = ex_matrix,
                                      signature_matrix = ref_list[[x]],
                                      signature_genes = row.names(ref_list[[x]])) %>% as.data.frame()
    dec_df_clr <- dec_df[1:10,]
    round(dec_df * 100, digits = 2)
  }
  )
  
  df <- append(df,df2)
  Celltype <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells",
                "Oligodendrocytes","T cells","Tumor cells","B cells","Mural cells")
  temp_df <- rbind(df[[1]][Celltype[1],],df2[[1]][Celltype[2],]) %>% as.data.frame()
  temp_df <- rbind(temp_df,df[[2]][Celltype[3],]) %>% as.data.frame()
  temp_df <- rbind(temp_df,df[[3]][Celltype[4],]) %>% as.data.frame()
  temp_df <- rbind(temp_df,df[[4]][Celltype[5],]) %>% as.data.frame()
  temp_df <- rbind(temp_df,df[[5]][Celltype[6],]) %>% as.data.frame()
  temp_df <- rbind(temp_df,df[[6]][Celltype[7],]) %>% as.data.frame()
  temp_df <- rbind(temp_df,df[[3]][Celltype[8],]) %>% as.data.frame()
  temp_df <- rbind(temp_df,df2[[2]][Celltype[9],]) %>% as.data.frame()
  temp_df <- rbind(temp_df,df[[7]][Celltype[10],]) %>% as.data.frame()
  
  re_res <- temp_df
  for (i in 1:length(colnames(temp_df))) {
    temp_df[,i] <- temp_df[,i] %>% as.double()
    re_res[,i] <- NA
    re_res[,i] <- round(temp_df[,i]/sum(temp_df[,i]),digits = 4)*100
  }
  write.csv(re_res, file = "~/QPC_GBM_result.csv")
  re_res
}