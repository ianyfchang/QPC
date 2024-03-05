# The Main function of QPC-GBM deconvolution

#' @param gene_expression_matrix could use TMM normalized and Rawcounts
#' @param ref_list gene expression from isolated cells, or a matrix of expression profile of cells.
#' @param df_source put your data source (Rawcounts or TMM)
#' @param Perm Set permutations for statistical analysis (≥100 permutations recommended).
#' @param QN boolean. Wheter to quantile normalize the data. Data should be normalized.
#' @param absolute Set to TRUE for CIBERSORT absolute mode.
 


#Deconvolute using CIBERSORT and Immunedeconv  
source("~/CIBERSORT_modified.R")

#input reference information================================
All_list <- list()
DC_list <- list("Cibersort", 611, "ALL", 100, "SCT")
EC_list <- list("EPIC", 611, "ALL", 20, "Raw")
MP_list <- list("Cibersort", 211, "ALL", 100, "SCT")
MG_list <- list("Cibersort", 611, "ALL", 100, "TMM")
NKT_list <- list("Cibersort", 211, "ALL", 100, "TPM")
OD_list <- list("Cibersort", 611, "DatasetWhole_MinDist", 50, "TPM")
TC_list <- list("Cibersort", 611, "ALL", 50, "Raw")
Tumor_list <- list("Cibersort", 611, "ALL", 100, "TMM")
BC_list <- list("EPIC", 611, "DatasetWhole_MinDist", 100, "TMM")
MC_list <- list("Cibersort", 411, "ALL", 100, "Raw")

All_list <- list(DC_list, EC_list, MP_list, MG_list, NKT_list, OD_list, TC_list, Tumor_list, BC_list, MC_list)


#main function
QPCdecon <- function(decon, FAM, Sampling, num, norMeth){
  ref_list <- list.files(path = paste0("~/"),pattern = paste0("Reference_GSE182109_11celltype_",FAM,"_",Sampling,"_n",num,"_",norMeth,collapse="|"),full.names = T,recursive = T)

  
  lapply(gene_expression_matrix, function(filename){
    datafile <- read_csv(filename) %>% as.data.frame() %>% column_to_rownames("gene_name")
    datafile$...1 <- NULL
    datafile$gene_name <- NULL
    
    lapply(ref_list, function(ref_file){
      ref_datafile <- read_csv(ref_file) %>% as.data.frame() %>% column_to_rownames("...1")
      print(ref_file)
      fam <- sapply(strsplit(ref_file, "_"), "[", 5)
      sampling_chr <- sapply(strsplit(ref_file, "_"), "[", 6)
      
      if (sampling_chr == "DatasetWhole") {
        sampling <- paste0(sapply(strsplit(ref_file, "_"), "[", 6),"_",sapply(strsplit(ref_file, "_"), "[", 7))
        num <- sapply(strsplit(ref_file [1], "_"), "[", 8)
        num <- gsub("n","", num)
        norMeth <- sapply(strsplit(sapply(strsplit(ref_file, "_"), "[", 9), "[.]"), "[", 1)
      } else {
        sampling <- sampling_chr
        num <- sapply(strsplit(ref_file [1], "_"), "[", 7)
        num <- gsub("n","", num)
        norMeth <- sapply(strsplit(sapply(strsplit(ref_file, "_"), "[", 8), "[.]"), "[", 1)
      }
      
      if (decon == "Cibersort") {
        dec_df <- CIBERSORT(mixture_file = datafile,
                            sig_matrix = ref_datafile,
                            perm = 100,
                            QN = FALSE, 
                            absolute = FALSE,
                            abs_method = "sig.score")
        dec_df_clr <- dec_df[1:10,]
        df <- round(dec_df_clr * 100, digits = 2)
      } else {
        # EPIC_custom
        dec_df <- deconvolute_epic_custom(gene_expression_matrix = datafile,
                                          signature_matrix = ref_datafile,
                                          signature_genes = row.names(ref_datafile))
        
        dec_df <- as.data.frame(dec_df)
        df <- round(dec_df * 100, digits = 2)
      }
      write.csv(df, paste0("~/DeconRes_ref_",fam,sampling,num,norMeth,"_",df_source,"_",decon,".csv"))
      print(paste0("~/DeconRes_ref_",fam,sampling,num,norMeth,"_",df_source,"_",decon,".csv"))
      gc()
    })
  })
}

gene_expression_matrix <- ""
df_Source <- "your data source"
for(k in 1:length(All_list)){
  decon <- All_list[[k]][[1]]
  FAM <- All_list[[k]][[2]]
  Sampling <- All_list[[k]][[3]]
  num <- All_list[[k]][[4]]
  norMeth <- All_list[[k]][[5]]
  res <- QPCdecon(decon, FAM, Sampling, num, norMeth)
}


#Merge result ==============================


MergeQPCres <- function(){
  decon_Meth <- c( "Cibersort","EPIC")
  FAM <- c(211,411,611)
  Sampling <- c("ALL","DatasetWhole_MinDist")
  num <- num <- c(20,50,100)
  norMeth <- c("TPM","TMM","LogNormalize","SCT","RawCount")
  Celltype <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells",
                "Oligodendrocytes","T cells","Tumor cells","B cells","Mural cells")
  factorlist <- list(Dendritic_cells = c(1,3,1,3,4),
                     Endothelial_cells = c(2,3,1,1,5),
                     Macrophage_GAMs = c(1,1,1,3,4),
                     Microglia_GAMs = c(1,3,1,3,2),
                     NKT_cells = c(1,1,1,3,1),
                     Oligo = c(1,3,2,2,1),
                     T_cells = c(1,3,1,2,5),
                     Tumor_cells = c(1,3,1,3,2),
                     B_cells = c(2,3,2,3,2),
                     Mural_cells = c(1,2,1,3,5))
  filelist <- c()
  for (j in 1:(length(Celltype)-1)) {
    a <- decon_Meth[factorlist[[j]][1]]
    b <- FAM[factorlist[[j]][2]]
    c <- Sampling[factorlist[[j]][3]]
    d <- num[factorlist[[j]][4]]
    e <- norMeth[factorlist[[j]][5]]
    
    decon_res <- read_csv(paste0("~/DeconRes_ref_",b,c,d,e,"_",df_source,"_",a,".csv"),
                          show_col_types = FALSE) %>% as.data.frame()
    print(paste0("~/DeconRes_ref_",b,c,d,e,"_",df_source,"_",a,".csv"))
    
    
    if(j==1){
      combine_res <- decon_res[which(decon_res$...1==Celltype[j]),]
    }else{
      combine_res <- rbind(combine_res,decon_res[which(decon_res$...1==Celltype[j]),])
    }
    filelist <- append(filelist,paste0(a,b,c,d,e))
    remove(a,b,c,d,e)
  }
  
  filelist
  
  combine_res[10,] <- c(Celltype[10],rep(0,length(colnames(combine_res))-1))
  row.names(combine_res) <- combine_res$...1
  combine_res$...1 <- NULL
  
  re_res <- combine_res
  
  for (i in 1:length(colnames(combine_res))) {
    combine_res[,i] <- combine_res[,i] %>% as.double()
    re_res[,i] <- NA
    re_res[,i] <- round(combine_res[,i]/sum(combine_res[,i]),digits = 4)*100
  }
  write.csv(re_res, file = paste0("~/DeconRes_ref_Best_mix_",df_source,"_CibersortEPIC.csv"))
  print(paste0("~/DeconRes_ref_Best_mix_",df_source,"_CibersortEPIC.csv"))
  
  
}


MerRes <- MergeQPCres()

