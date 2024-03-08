# QPC-GBM Deconvolution   
## Overview

QPC-GBM is a computational pipeline for the **Q**uantifying the **P**roportions of **g**lio**b**lastoma **m**ultiforme cell (GBM) Cell Type from human single cell RNA sequencing data.
It is for unified access to computational methods for estimating GBM fractions from bulk RNA sequencing data.

### Reference database
Of all Reference database, we used single cell RNA sequecing in four normalized methods.

[Raw read counts]


TMM (Trimmed mean of M-values)     
TPM (Transcripts per million)      
SCT



If you use this pipeline in your work, please cite both our paper and the method(s) you are using.

### RNA-seq data normalization 
After various analysis, we recommend using Raw counts or TMM normalized sequencing data.


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

An example of the input data in R is shown below.
```R
dim(data)
[1] 61155    15
data[1:5,1:5]


             TCGA-26-5133  TCGA-HT-7902  TCGA-VM-A8CF  TCGA-14-1823  TCGA-DU-6393
A1BG-AS1           77           89           38            7           72
A1CF                0            3            1            0            4
A2M             75396        36243        20502        20106        16941
A2M-AS1            26          163           19           29           31
A2ML1             170          717          290          280          351          

```


### Basic usage
You can install the released version of immunedeconv from [github](https://github.com/) with:
```R
# Loading required package
install.packages("remotes")
remotes::install_github("omnideconv/immunedeconv")

library(tibble)
library(tidyverse)
library(readr)
library(xlsx)
```



### Reference
1. Newman, A.M., Steen, C.B., Liu, C.L. et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nat Biotechnol 37, 773–782 (2019). <br> (https://doi.org/10.1038/s41587-019-0114-2)
2. Sturm, G., Finotello, F., Petitprez, F., Zhang, J. D., Baumbach, J., Fridman, W. H., List, M., Aneichyk, T. (2019).Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology.Bioinformatics, 35(14), i436-i445.  <br>(https://doi.org/10.1093/bioinformatics/btz363)
3. Neftel C, Laffy J, Filbin MG, Hara T et al. An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma. Cell 2019 Aug 8;178(4):835-849.e21. <br> (https://doi.org/10.1016/j.cell.2019.06.024.)
4. Abdelfattah N, Kumar P, Wang C, Leu JS et al. Single-cell analysis of human glioma and immune cells identifies S100A4 as an immunotherapy target. Nat Commun 2022 Feb 9;13(1):767. <br> (https://doi.org/10.1038/s41467-022-28372-y)
5. Richards, L.M., Whitley, O.K.N., MacLeod, G. et al. Gradient of Developmental and Injury Response transcriptional states defines functional vulnerabilities underpinning glioblastoma heterogeneity. Nat Cancer 2, 157–173 (2021). <br> (https://doi.org/10.1038/s43018-020-00154-9)
