# QPC-GBM Deconvolution   
## Overview

QPC-GBM is a computational pipeline for the **Q**uantifying the **P**roportions of **g**lio**b**lastoma **m**ultiforme cell (GBM) Cell Type from human single cell RNA sequencing data.
It is for unified access to computational methods for estimating GBM fractions from bulk RNA sequencing data.     
                         
     
### Workflow of build QPC-GBM Deconvolution   
![image](https://github.com/ianyfchang/QPC-GBM/blob/master/Fig/fig.jpg)

                                       
1. Define 10 Cell types by Using scRNA-seq dataset: 

                                              
    |Cell Type|
      |---------|
      |Dendritic cells|
      |Endothelial cells| 
      |Macrophage-like GAMs|
      |Microglia-like GAMs| 
      |NKT-like cells|
      |Oligodendrocytes| 
      |T cells |
      |Tumor cells| 
      |B cells|
      |Mural cells|   



                                

2. Building reference matrix and bulk RNA matrix from scRNA-seq dataset and normalization:                                         
- The methods for single cell reference dataset normalization
  ```R
  Raw read counts                                     
  TPM                                
  TMM                                   
  LogNormalize                                                 
  SCT                                                        
  ```


  
- The methods for Psudo bulk RNA-seq normalization
  ```R 
  TMM
  TPM
  logTPM
  RPKM
  ntd
  vst
  NormCount
  RawCount
  ```



3. Using different deconvolution methods for estimate the proportion of the different cell types from gene expression data:
- Deconvolution methods                                    
  ```R      
  CIBERSORTx                                                                                
  EPIC
  ConsensusTME
  ```

   
(https://cibersortx.stanford.edu/)  
(https://epic.gfellerlab.org/)
(https://github.com/cansysbio/ConsensusTME)



4. Scoring with two methods that is RMSE and PearsonR, then select each cell type with the best score.
* The Root Mean Squared Error (RMSE) measures the average difference between values predicted and the actual values by a model. It provides an estimation of how well the model is able to predict the target value.
* Pearson correlation coefficient (r) is a number between –1 and 1 that measures the strength and direction of the relationship between two variables.
  


### Basic usage
You can install the released version of immunedeconv from [github](https://github.com/) with:
```R
# Loading required package
install.packages("remotes")
remotes::install_github("omnideconv/immunedeconv")

library(immunedeconv)
library(tibble)
library(tidyverse)
library(readr)
library(xlsx)                
```
                                    
### RNA-seq data normalization 
After various analysis, we recommend using Raw read counts or TMM normalized sequencing data.
An example of the input data is a gene × sample gene expression matrix in R is shown below.             
```R
dim(data)
[1] 61155    15
data[1:5,1:5]      

# Raw read counts
             TCGA-26-5133  TCGA-HT-7902  TCGA-VM-A8CF  TCGA-14-1823  TCGA-DU-6393
A1BG-AS1           77           89           38            7           72
A1CF                0            3            1            0            4
A2M             75396        36243        20502        20106        16941
A2M-AS1            26          163           19           29           31
A2ML1             170          717          290          280          351          

```
                        
### Example      
For this example, we use a dataset of GBM patients from TCGA.
```R

# make deconvolution results                
DeRes <- QPCdecon(gene_expression_matrix_path,ref_path)

# merge different cell type data
MerRes <- MergeQPCres()
```




Take a look at the results table
```R     
##                      TCGA-26-5133 TCGA-HT-7902 TCGA-VM-A8CF TCGA-14-1823 TCGA-DU-6393
## Dendritic cells              2.88         1.99         0.86         1.74         5.42
## Endothelial cells           10.06         2.87         4.10         9.72        22.64
## Macrophage-like GAMs         0.00        11.91        21.18        18.65         0.00
## Microglia-like GAMs          0.00         0.00         0.00        20.95         0.00
## NKT-like cells               0.00         0.00         0.00         0.00         0.00
## Oligodendrocytes             1.49        82.71        21.28         0.53        13.68
## T cells                      0.00         0.00         0.00         0.00         0.00
## Tumor cells                 84.00         0.50        52.57        45.03        58.26
## B cells                      0.00         0.01         0.00         0.00         0.00
## Mural cells                  1.58         0.00         0.00         3.39         0.00



```

                                                                                                                                 
If your Bulk RNA data is generated by other Normaliztion methods and RawCount or TMM-Normalizd is not available, please contact the contributor.                                                                                                                                 
                                                                                                                                 
## About reference database            
                                                                   
<div align = center><img src= "https://github.com/ianyfchang/QPC-GBM/blob/master/Fig/Fig1.png" alt = "Image" width="610px"></div>  
                                                          
To build database, the following steps will be performed:    

1. Define single cell RNA sequencing cell-type
 

                                                                              

2. Of all Reference database, we used single cell RNA sequecing in four normalized methods for test which methods are standard and widely used in scRNA-seq analysis. If you use this pipeline in your work, please cite both our paper and the method(s) you are using.                                      

         
Generation TPM normalized data          
```R
# Convert counts to TPM by EDASeq package
# Loading required package
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EDASeq")

library(EDASeq)
gene_Length <- getGeneLengthAndGCContent(gene_list, org, mode=c("biomart", "org.db"))
data <- data/ gene_Length
TPM <- as.data.frame(t( t(data) / apply(data, 2, sum, na.rm = TRUE) ) * 1e6)
```
         
Generation TMM normalized data
```R
# Convert counts to TMM by EdgeR package
# Loading required package
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR")

library(edgeR)

# make the DGEList:
dgelist <- DGEList(counts = data, group = colnames(data))
keep <- rowSums(cpm(dgelist )>1) >= 2
dgelist <- dgelist[keep, keep.lib.sizes=FALSE]
# calculate TMM normalization factors
dgelist <- calcNormFactors(dgelist,method = "TMM")
# get the normalized counts
dgelist <- cpm(dgelist)
```

Generation SCT, LogNormalize and Raw read counts. To obtain these data, we will be mainly using functions available in the Seurat package. Apart from Seurat package, we also need sctransform package which as a more accurate method of normalizing, estimating the variance of the raw filtered data.                       
```R
# Install sctransform from CRAN
install.packages("sctransform")


library(Seurat)
library(sctransform)

seurat <- CreateSeuratObject(counts = seurat_data)
# run sctransform
seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", verbose = FALSE)

# SCT normalized data
seurat[["SCT"]]$scale.data

# LogNormalize data
seurat[["SCT"]]$data

# Raw read counts
seurat[["SCT"]]$counts
```
          
3. The method for obtain feature genes. Two cell populations were used to identify the feature genes, one included the whole cell population (ALL), and the other was a 10000-cell population composed from 1000 cells per cell type that showed the smallest distance to the cell type geometric mean on UMAP (MinDist).
  
4. Find all markers for per cell types with different parameters. We used FindAllMarkers() which is find markers for every cluster compared to all remaining cells in Seurat to find genes. Therefore, we define it as 211, 411 and 611.
```R
# Use different parameters for min.pct
seurat.markers <- FindAllMarkers(seurat,
                                 only.pos = TRUE,
                                 assay = "SCT",
                                 test.use = "wilcox",
                                 recorrect_umi = FALSE,
                                 min.pct = 0.2, # 0.2, 0.4, 0.6
                                 logfc.threshold = 0.1,
                                 min.diff.pct = 0.1,
                                 return.thresh = 0.05)
```

5. Use different feature gene number including 20,50 and 100 for per cell type. Select values by smallest Adjusted p-value (p_val_adj), target expression mean, genes that show a  difference in the fraction of detection between the two groups (pct-diff) and log fold-change of the average expression (avg-log2FC).           

                                                                            



### Session info     
```R
sessionInfo()
```   

```R
## R version 4.1.3 (2022-03-10)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 22631)

## Matrix products: default     

## locale:
## [1] LC_COLLATE=Chinese (Traditional)_Taiwan.950  LC_CTYPE=Chinese (Traditional)_Taiwan.950   
## [3] LC_MONETARY=Chinese (Traditional)_Taiwan.950 LC_NUMERIC=C                                
## [5] LC_TIME=Chinese (Traditional)_Taiwan.950

## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     

## other attached packages:
## [1] edgeR_3.36.0       limma_3.50.3       immunedeconv_2.1.0 EPIC_1.1.6         xlsx_0.6.5        
## [6] lubridate_1.9.2    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.2        purrr_1.0.1       
## [11] readr_2.1.4        tidyr_1.3.0        ggplot2_3.4.4      tidyverse_2.0.0    tibble_3.2.1      

## loaded via a namespace (and not attached):
## [1] nlme_3.1-155           bitops_1.0-7           matrixStats_0.63.0     bit64_4.0.5           
## [5] filelock_1.0.2         progress_1.2.3         httr_1.4.7             GenomeInfoDb_1.30.1   
## [9] data.tree_1.1.0        tools_4.1.3            utf8_1.2.3             R6_2.5.1              
## [13] DBI_1.2.0              BiocGenerics_0.40.0    mgcv_1.8-39            colorspace_2.1-0      
## [17] mMCPcounter_1.1.0      withr_2.5.2            tidyselect_1.2.0       prettyunits_1.2.0     
## [21] preprocessCore_1.56.0  bit_4.0.5              curl_5.0.0             compiler_4.1.3        
## [25] cli_3.6.1              Biobase_2.54.0         xml2_1.3.3             ComICS_1.0.4          
## [29] scales_1.2.1           genefilter_1.76.0      rappdirs_0.3.3         digest_0.6.31         
## [33] XVector_0.34.0         pkgconfig_2.0.3        dbplyr_2.3.3           fastmap_1.1.1         
## [37] readxl_1.4.2           rlang_1.1.0            rstudioapi_0.15.0      RSQLite_2.3.1         
## [41] generics_0.1.3         testit_0.13            BiocParallel_1.28.3    RCurl_1.98-1.12       
## [45] magrittr_2.0.3         GenomeInfoDbData_1.2.7 Matrix_1.5-4           Rcpp_1.0.10           
## [49] munsell_0.5.0          S4Vectors_0.32.4       fansi_1.0.4            lifecycle_1.0.4       
## [53] stringi_1.7.12         MASS_7.3-55            zlibbioc_1.40.0        BiocFileCache_2.2.1   
## [57] grid_4.1.3             blob_1.2.4             parallel_4.1.3         crayon_1.5.2          
## [61] lattice_0.20-45        Biostrings_2.62.0      splines_4.1.3          annotate_1.72.0       
## [65] xlsxjars_0.6.1         hms_1.1.3              KEGGREST_1.34.0        locfit_1.5-9.7        
## [69] pillar_1.9.0           biomaRt_2.50.3         stats4_4.1.3           XML_3.99-0.14         
## [73] glue_1.6.2             remotes_2.4.2.1        png_0.1-8              vctrs_0.6.1           
## [77] tzdb_0.3.0             cellranger_1.1.0       gtable_0.3.4           cachem_1.0.7          
## [81] xtable_1.8-4           survival_3.2-13        rJava_1.0-6            AnnotationDbi_1.56.2  
## [85] memoise_2.0.1          IRanges_2.28.0         sva_3.42.0             timechange_0.2.0      

```





### Reference
1. Newman, A.M., Steen, C.B., Liu, C.L. et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nat Biotechnol 37, 773–782 (2019). <br> (https://doi.org/10.1038/s41587-019-0114-2)
2. Sturm, G., Finotello, F., Petitprez, F., Zhang, J. D., Baumbach, J., Fridman, W. H., List, M., Aneichyk, T. (2019).Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology.Bioinformatics, 35(14), i436-i445.  <br>(https://doi.org/10.1093/bioinformatics/btz363)
3. Neftel C, Laffy J, Filbin MG, Hara T et al. An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma. Cell 2019 Aug 8;178(4):835-849.e21. <br> (https://doi.org/10.1016/j.cell.2019.06.024.)
4. Abdelfattah N, Kumar P, Wang C, Leu JS et al. Single-cell analysis of human glioma and immune cells identifies S100A4 as an immunotherapy target. Nat Commun 2022 Feb 9;13(1):767. <br> (https://doi.org/10.1038/s41467-022-28372-y)
5. Richards, L.M., Whitley, O.K.N., MacLeod, G. et al. Gradient of Developmental and Injury Response transcriptional states defines functional vulnerabilities underpinning glioblastoma heterogeneity. Nat Cancer 2, 157–173 (2021). <br> (https://doi.org/10.1038/s43018-020-00154-9)
6. Risso D, Schwartz K, Sherlock G, Dudoit S (2011). “GC-Content Normalization for RNA-Seq Data.” BMC Bioinformatics, 12(1), 480. <br> (https://doi.org/10.1186/1471-2105-12-480)
7. Chen B, Khodadoust MS, Liu CL, Newman AM, Alizadeh AA. Profiling Tumor Infiltrating Immune Cells with CIBERSORT. Methods Mol Biol. 2018;1711:243-259.PMID: 29344893; PMCID: PMC5895181. <br> (https://doi.org/10.1007/978-1-4939-7493-1_12)
8. Racle, J., Gfeller, D. (2020). EPIC: A Tool to Estimate the Proportions of Different Cell Types from Bulk Gene Expression Data. In: Boegel, S. (eds) Bioinformatics for Cancer Immunotherapy. Methods in Molecular Biology, vol 2120. Humana, New York, NY. <br> (https://doi.org/10.1007/978-1-0716-0327-7_17)
9. Alejandro Jiménez-Sánchez, Oliver Cast, Martin L. Miller; Comprehensive Benchmarking and Integration of Tumor Microenvironment Cell Estimation Methods. Cancer Res 15 December 2019; 79 (24): 6238–6246. <br>(https://doi.org/10.1158/0008-5472.CAN-18-3560)
