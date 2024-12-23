# Optimal Deconvolution for GBM   


This is a computational deconvolution method for estimating cell type proportions in bulk RNA from GBM samples. The 10 cell types used in optimal deconvolution are defined and characterized from three public scRNA-seq datasets ([GSE182109](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182109) [1], [GSE131928](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928) [2], and the NC dataset by Richards, et al. [3]). Each cell type proportion is measured using different reference matrices with deconvolution tools of either CIBERSORT or EPIC, followed by combining each result and re-scaling to get the final cell type proportions. The method could be applied to bulk RNA-seq data generated from both poly(A) and exome capture RNA sequencing. Bulk RNA-seq data of raw count gives best deconvolution performance, whereas TMM and log normalized readcount using DESeq2 also show good results.
                         
     
![image](https://github.com/ianyfchang/QPC-GBM/blob/master/Fig/Github_fig.jpg)

                                       
## Cell types in GBM deconvolution
| Cell Type | Feature Genes
| --------- | --------------
|Dendritic cells|CD1D, FCER1A
|Endothelial cells|FLT1, PECAM1, VWF
|Macrophage-like GAMs|CD68, CD163
|Microglia-like GAMs|TMEM119, CX3CR1
|NKT-like cells|CD3E, NKG7
|Oligodendrocytes|PLP1, NKX6-2
|T cells |CD3E, IL7R
|Tumor cells|GAP43, PTPRZ1
|B cells|MS4A1, CD79A
|Mural cells|DCN, COL1A2


## Deconvolution methods used by GBM deconvolution
| Tool | Algorithm | Reference
| ---- | --------- | ----------
| [CIBERSORT](https://cibersortx.stanford.edu/) | &#965;-suppport vector regression | Newman et al. (2015) [4]
| [EPIC](https://cibersortx.stanford.edu/) | constrained least square regression | Racle et al. (2020) [5]


## Before performing GBM deconvolution
### Packages required for running deconvolution
```R
# Loading required package
# For CIBERSORT
library(e1071)
library(parallel)
library(preprocessCore)

# For GBM deconvolution
install.packages("remotes")
remotes::install_github("omnideconv/immunedeconv")

library(immunedeconv)
library(tibble)
library(tidyverse)
library(readr)
library(xlsx)                
```
                                    
### Input bulk RNA-seq normalization 
Raw read count without transformation or data transformed by TMM or log normalization are recommend to be used for GBM deconvolution.
An example of the input data is shown below: a gene × sample gene expression matrix.             
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


## Peforming GBM deconvolution      
For example, we use a bulk RNA-seq dataset of TCGA-GBM.
### Example dataset
```R
# Download example dataset
library (readr)
urlfile <- "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Sample/TCGA_Rawreadcounts.csv"
mydata <- read_csv(url(urlfile)) %>% as.data.frame() %>% column_to_rownames("...1")

# load functions
source("https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/R/CIBERSORT_modified.R")
source("https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/R/Optimal-GBM_main.R")

```
### Reference matrices
```R
# Download reference matrices
ref_path <- c("https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_DC_ori.csv",
              "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_EC_ori.csv",
              "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_MP_ori.csv",
              "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_MGTumor_ori.csv",
              "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_NKT_ori.csv",
              "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_Oligo_ori.csv",
              "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_TC_ori.csv",
              "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blobmaster/Reference_DB/Reference_BC_ori.csv",
              "https://raw.githubusercontent.com/ianyfchangOptimal-GBM/blob/master/Reference_DB/Reference_MC_ori.csv")
ref_List <- lapply(ref_path, function(x) read_csv(x) %>% as.data.frame() %>% column_to_rownames("...1"))
names(ref_List) <- c("DC","EC","MP","MGTumor","NKT","Oligo","TC","BC","MC")
remove(ref_path)

## Please do not change the order of the references
## The SCTransform was used to normalize and transform the reference matriced for deconvolution aimed to estimated the proportions
## of dendritice cells and macrophage-like GAMs. The referenc matrices contain negative values that may cause errors when performing
## deconvolution using CIBERSORT. If the errors occurred at the run of dendritic cells or macrophage-like GAMs, you can replace the
## two reference matrices with the modified ones and try again.

# Download modified reference matrices
ref_path <- c("https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_DC_mod.csv",
              "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_EC_ori.csv",
              "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_MP_mod.csv",
              "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_MGTumor_ori.csv",
              "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_NKT_ori.csv",
              "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_Oligo_ori.csv",
              "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_TC_ori.csv",
              "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_BC_ori.csv",
              "https://raw.githubusercontent.com/ianyfchang/Optimal-GBM/blob/master/Reference_DB/Reference_MC_ori.csv")
ref_List <- lapply(ref_path, function(x) read_csv(x) %>% as.data.frame() %>% column_to_rownames("...1"))
names(ref_List) <- c("DC","EC","MP","MGTumor","NKT","Oligo","TC","BC","MC")
remove(ref_path)

```

### Running GBM deconvolution
Enter the expression matrix (raw count, TMM, or Lognormalization) and the reference list.
```R
# run deconvolution                
DeRes <- GBMdecon(ex_matrix, ref_List)

## The function returns a cell type proportion table and save a "Optimal_GBM_result.csv" file
```

## Results   
```R
DeRes <- GBMdecon(mydata, ref_list)
  
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

                                                                                                                                 
If your Bulk RNA data was generated by normaliztion methods other than raw count, TMM, or lognormalization, please contact the contributor for additional information.


## Session info     
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


## References
1. Abdelfattah N, et al. Single-cell analysis of human glioma and immune cells identifies S100A4 as an immunotherapy target. Nat Commun 13, 767 (2022).
2. Neftel C, et al. An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma. Cell 178, 835-849 e821 (2019).
3. Richards LM, et al. Gradient of Developmental and Injury Response transcriptional states defines functional vulnerabilities underpinning glioblastoma heterogeneity. Nat Cancer 2, 157-173 (2021).
4. Newman AM, et al. Robust enumeration of cell subsets from tissue expression profiles. Nat Methods 12, 453-457 (2015).
5. Racle J, Gfeller D. EPIC: A Tool to Estimate the Proportions of Different Cell Types from Bulk Gene Expression Data. Methods Mol Biol 2120, 233-248 (2020).
