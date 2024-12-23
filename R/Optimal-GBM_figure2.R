library(ggplot2)
library(magrittr)
library(stringr)
library(dplyr)
library(readr)
library(Matrix)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)
library(readxl)
library(tidyverse)

#Cell composition proportion-----------------------------------------------------------------------------------------------------------------------------
CombineTable <- read_csv("~/TCGA_Survival_CellType_145RawCount_CombineTable.csv") %>% as.data.frame()
row.names(CombineTable) <- CombineTable$...1
CombineTable$...1 <- NULL

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
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.y = element_blank(),
        #axis.title.y = element_text(color = "black", size = 16),
        axis.line.x = element_line(color = "white"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1.1, vjust = 0.5, color = "black"), 
        text = element_text(size = 14),
        legend.text = element_text(color = "black", size = 15),
        legend.title = element_text(color = "black", size = 16),
        legend.margin = margin(),
        legend.spacing.x = unit(0.25, "cm"),
        legend.key.height = unit(0.25, "cm"),legend.key.width = unit(0.4, "cm")) 

p1

#Survival analysis, best p value------------------------------------------------------------------------------------------------------
library(survival)
library(survminer)

TCGA <- read_csv("~/TCGA_Survival_CellType_145RawCount_CombineTable.csv") %>% as.data.frame()
  
oridata <- TCGA
low_25 <- ceiling(length(oridata$PATIENT_ID)*0.25)
high_25 <- floor(length(oridata$PATIENT_ID)*0.75)

oridata <- oridata[order(oridata$`Dendritic cells`),] %>% as.data.frame()
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

#Survival analysis, survival graph-------------------------------------------------------------------------------------------------------------------
TCGA <- read_csv("~/TCGA_Survival_CellType_145RawCount_CombineTable.csv") %>% as.data.frame()
row.names(TCGA) <- TCGA$...1
colnames(TCGA)[1] <- "PATIENT_ID"
CGGA <- read_csv("~/CGGA_Survival_CellType_153_RawCount_CombineTable.csv") %>% as.data.frame()
row.names(CGGA) <- CGGA$...1
colnames(CGGA)[1] <- "PATIENT_ID"
CGMH <- read_csv("~/CGMH_Survival_CellType_40GBM_RawCount_CombineTable.csv") %>% as.data.frame()
row.names(CGMH) <- CGMH$...1
colnames(CGMH)[1] <- "PATIENT_ID"

Low_numner_OS <- data.frame(TCGA = c(98,92,43,41),
                            CGMH = c(26,16,12,15),
                            CGGA = c(114,81,44,43),
                            row.names = c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs"))
Low_numner_PFS <- data.frame(TCGA = c(52,92,72,78),
                             CGMH = c(19,22,14,15),
                             row.names = c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs"))
celltype <- c("Dendritic cells","Macrophage-like GAMs","Microglia-like GAMs")
oridata <- list(TCGA,CGMH,CGGA)
dataset_name <- c("TCGA","CGMH","CGGA")

k <- 1
OS_fig <- list()
PFS_fig <- list()
pval_OS <- c()
median_H_OS <- c()
median_L_OS <- c()
pval_PFS <- c()
median_H_PFS <- c()
median_L_PFS <- c()
for (df in c(1:3)) {
  for (ctype in celltype) {
    temp_df <- oridata[[df]][order(oridata[[df]][,which(colnames(oridata[[df]])==ctype)]),] %>% as.data.frame()
    i <- Low_numner_OS[ctype,dataset_name[df]]
    temp_df$OS_group <- "H"
    temp_df$OS_group[1:i] <- "L"
    fit_OS <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ OS_group, data = temp_df)
    pvalue_OS <- surv_pvalue(fit = fit_OS, data = temp_df)
    pval_OS <- append(pval_OS, round(surv_pvalue(fit = fit_OS)[1,2],digits = 5))
    median_H_OS <- append(median_H_OS, surv_median(fit_OS)[1,2])
    median_L_OS <- append(median_L_OS, surv_median(fit_OS)[2,2])
    OS_fig[[k]] <- ggsurvplot(fit_OS,
                              pval = paste0("p=",format(pvalue_OS$pval, digits = 4)),
                              pval.coord = c(50,0.8),
                              risk.table = TRUE,
                              #risk.table = FALSE,
                              legend.title = paste0(dataset_name[df],"--",ctype),
                              surv.median.line = c("hv"),
                              legend = c("top"), 
                              ylab="OS Probability",
                              xlab="Months after first diagnosis",
                              palette= c("#B2182B","#2166AC"),
                              legend.labs = c("High","Low"),
                              font.title = 16,
                              font.x =  16,
                              font.y = 16,
                              font.tickslab = 16,
                              tables.y.text = FALSE)
    if(dataset_name[df]!="CGGA"){
      j <- Low_numner_PFS[ctype,dataset_name[df]]
      temp_df$PFS_group <- "H"
      temp_df$PFS_group[1:j] <- "L"
      fit_PFS <- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ PFS_group, data = temp_df)
      pvalue_PFS <- surv_pvalue(fit = fit_PFS, data = temp_df)
      pval_PFS <- append(pval_PFS, round(surv_pvalue(fit = fit_PFS)[1,2],digits = 5))
      median_H_PFS <- append(median_H_PFS, surv_median(fit_PFS)[1,2])
      median_L_PFS <- append(median_L_PFS, surv_median(fit_PFS)[2,2])
      PFS_fig[[k]] <- ggsurvplot(fit_PFS,
                                 pval = paste0("p=",format(pvalue_PFS$pval, digits = 4)),
                                 pval.coord = c(20,0.8),
                                 risk.table = TRUE,
                                 #risk.table = FALSE,
                                 legend.title = paste0(dataset_name[df],"--",ctype),
                                 surv.median.line = c("hv"),
                                 legend = c("top"), 
                                 ylab="PFS Probability",
                                 xlab="Months after first diagnosis",
                                 palette= c("#B2182B","#2166AC"),
                                 legend.labs = c("High","Low"),
                                 font.title = 16,
                                 font.x =  16,
                                 font.y = 16,
                                 font.tickslab = 16,
                                 tables.y.text = FALSE)
      remove(j,fit_PFS,pvalue_PFS)
    }
    remove(temp_df,i,fit_OS,pvalue_OS)
    k <- k+1
  }
}

arrange_ggsurvplots(OS_fig, print = TRUE,
                    ncol = 3, nrow = 4, risk.table.height = 0.3, surv.plot.height = 0.7)
arrange_ggsurvplots(PFS_fig, print = TRUE,
                    ncol = 2, nrow = 4, risk.table.height = 0.3, surv.plot.height = 0.7)

#3 datasets cell type percentage comparasion--------------------------------------------------------------------------------------
TCGA <- read_csv("~/TCGA_Survival_CellType_145RawCount_CombineTable.csv") %>% as.data.frame()
row.names(TCGA) <- TCGA$...1
TCGA$...1 <- NULL
TCGA <- TCGA[,c(1:10)] %>% as.data.frame()
TCGA$Dataset <- "TCGA"
TCGA <- TCGA %>% t() %>% as.data.frame()
TCGA$Celltype <- row.names(TCGA)

CGGA <- read_csv("~/CGGA_Survival_CellType_153_RawCount_CombineTable.csv") %>% as.data.frame()
row.names(CGGA) <- CGGA$...1
CGGA$...1 <- NULL
CGGA <- CGGA[,c(1:10)] %>% as.data.frame()
CGGA$Dataset <- "CGGA"
CGGA <- CGGA %>% t() %>% as.data.frame()
CGGA$Celltype <- row.names(CGGA)

CGMH <- read_csv("~/CGMH_Survival_CellType_40GBM_RawCount_CombineTable.csv") %>% as.data.frame()
row.names(CGMH) <- CGMH$...1
CGMH$...1 <- NULL
CGMH <- CGMH[,c(1:10)] %>% as.data.frame()
CGMH$Dataset <- "CGMH"
CGMH <- CGMH %>% t() %>% as.data.frame()
CGMH$Celltype <- row.names(CGMH)

MergeDataset <- merge(TCGA,CGGA, by.x = "Celltype", by.y = "Celltype", all = TRUE) %>% as.data.frame()
MergeDataset <- merge(MergeDataset, CGMH, by.x = "Celltype", by.y = "Celltype", all = TRUE) %>% as.data.frame()
targetOrder <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells","Oligodendrocytes",
                 "T cells","Tumor cells","B cells","Mural cells","Dataset")
MergeDataset <- MergeDataset[match(targetOrder, MergeDataset$Celltype),] %>% as.data.frame()
row.names(MergeDataset) <- MergeDataset$Celltype
MergeDataset$Celltype <- NULL
MergeDataset <- MergeDataset %>% t() %>% as.data.frame()

for (i in 1:(ncol(MergeDataset)-1)) {
  tempdf <- MergeDataset[,c(i,ncol(MergeDataset))] %>% as.data.frame()
  colnames(tempdf) <- c("Perc","Dataset")
  tempdf$Celltype <- targetOrder[i]
  if(i==1){
    temptemp <- tempdf
  }else{
    temptemp <- rbind(temptemp,tempdf) %>% as.data.frame()
  }
}
temptemp$Perc <- as.numeric(temptemp$Perc)
temptemp$Celltype <- factor(temptemp$Celltype, levels = targetOrder[1:10])
temptemp$Dataset <- factor(temptemp$Dataset, levels = c("TCGA","CGGA","CGMH"))

p1 <- ggplot(data = temptemp, aes(x = Celltype, y = Perc, fill = Dataset)) +
  geom_boxplot(position=position_dodge(0.9)) +
  scale_fill_manual(values=c("#B2182B","#F4A582","#2166AC")) +
  labs(x = "", y = "Percentage", fill = "Datasets") +
  theme_classic() +
  theme(legend.position = "right",
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black", size = 16),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.text.x = element_blank(),
        text = element_text(size = 14),
        legend.text = element_text(color = "black", size = 15),
        legend.title = element_text(color = "black", size = 16),
        legend.spacing.x = unit(0.2, "cm"),
        legend.key.height = unit(0.8, "cm"),legend.key.width = unit(0.8, "cm")) 
p1

#Deconvolution result summarize figures----------------------------------------------------------------------------------------------------------
library(readxl)

BulkNorm_Meth <- c("TPM","TMM","vst","NormCount","RawCount","ntd","logTPM","RPKM")

for (i in 1:length(BulkNorm_Meth)) {
  for (k in c("RMSE","Pearson")) {
    oridata <- read_excel(paste0("~/SumTable_BulkRNA_",k,"_AllDataset_Celltype_",BulkNorm_Meth[i],".xlsx")) %>% as.data.frame()
    print(paste0("open file--~/SumTable_BulkRNA_",k,"_AllDataset_Celltype_",BulkNorm_Meth[i],".xlsx"))
    row.names(oridata) <- oridata$...1
    oridata$...1 <- NULL
    if(k=="RMSE"){
      oridata <- oridata[-grep(pattern = "CibersortEPIC_Best_mix_", row.names(oridata)),] %>% as.data.frame()
    }else{
      oridata <-oridata
    }
    tempdf_all <- oridata[-grep(pattern = "DatasetWhole_MinDist",row.names(oridata)),c("SumUp_MergeDataset","SumUp_NC","SumUp_GSE131928","SumUp_GSE182109")] %>% as.data.frame()
    tempdf_dwmd <- oridata[grep(pattern = "DatasetWhole_MinDist",row.names(oridata)),c("SumUp_MergeDataset","SumUp_NC","SumUp_GSE131928","SumUp_GSE182109")] %>% as.data.frame()
    a <- str_split(row.names(tempdf_all), pattern = "_")
    b <- str_split(row.names(tempdf_dwmd), pattern = "_")
    for (j in 1:nrow(tempdf_all)) {
      tempdf_all[j,c(5:10)] <- a[[j]]
      tempdf_dwmd[j,c(5:10)] <- b[[j]][c(1,2,4,5,6,7)]
    }
    if(k=="RMSE"){
      tempdf <- rbind(tempdf_all,tempdf_dwmd) %>% as.data.frame()
      print("Done RMSE")
    }else{
      tempdf <- rbind(tempdf,tempdf_all) %>% as.data.frame()
      tempdf <- rbind(tempdf,tempdf_dwmd) %>% as.data.frame()
      print("Done PearsonR")
    }
    remove(oridata,tempdf_all,tempdf_dwmd,a,b)
  }
  colnames(tempdf)[c(5:10)] <- c("DeconMeth","Min.pct","Sampling","GeneNumber","RefNorm","BulkRNA")
  tempdf$BulkNorm <- BulkNorm_Meth[i]
  tempdf$BulkCellNumber <- c(rep(100,135),rep(1000,135),rep(2500,135),rep(6000,135),
                             rep(100,135),rep(1000,135),rep(2500,135),rep(6000,135),
                             rep(100,135),rep(1000,135),rep(2500,135),rep(6000,135),
                             rep(100,135),rep(1000,135),rep(2500,135),rep(6000,135))
  tempdf$Value <- c(rep("RMSE",1080),rep("PearsonR",1080))
  if(i==1){
    mergetable <- tempdf
  }else{
    mergetable <- rbind(mergetable,tempdf) %>% as.data.frame()
  }
  print(paste0("Mergefile--",BulkNorm_Meth[i]))
  remove(tempdf)
}
mergetable$Min.pct <- str_replace(mergetable$Min.pct,pattern = "211", replacement = "0.2")
mergetable$Min.pct <- str_replace(mergetable$Min.pct,pattern = "411", replacement = "0.4")
mergetable$Min.pct <- str_replace(mergetable$Min.pct,pattern = "611", replacement = "0.6")
mergetable$ReverseRMSE100_MergeDataset <- 100/mergetable$SumUp_MergeDataset
mergetable$ReverseRMSE100_NC <- 100/mergetable$SumUp_NC
mergetable$ReverseRMSE100_GSE131928 <- 100/mergetable$SumUp_GSE131928
mergetable$ReverseRMSE100_GSE182109 <- 100/mergetable$SumUp_GSE182109

fileinuse <- mergetable[which(mergetable$Value=="RMSE"),] %>% as.data.frame()
fileinuse$BulkNorm <- factor(fileinuse$BulkNorm, levels = c("RawCount","vst","ntd","NormCount","TMM","TPM","logTPM","RPKM"))
fileinuse$DeconMeth <- factor(fileinuse$DeconMeth, levels = c("Cibersort","EPIC","ConsesnusTME"))
fileinuse$Min.pct <- factor(fileinuse$Min.pct, levels = c("0.2","0.4","0.6"))
fileinuse$Sampling <- factor(fileinuse$Sampling, levels = c("ALL","MinDist"))
fileinuse$GeneNumber <- factor(fileinuse$GeneNumber, levels = c(20,50,100))
fileinuse$RefNorm <- factor(fileinuse$RefNorm, levels = c("RawCount","LogNormalize","SCT","TPM","TMM"))

library("RColorBrewer")
#"DeconMeth","Min.pct","Sampling","GeneNumber","RefNorm","BulkNorm","BulkCellNumber"
S_fig <- ggplot(data = fileinuse, aes(x = RefNorm, y = ReverseRMSE100_MergeDataset, fill = BulkNorm)) +
  geom_boxplot(position=position_dodge(0.9)) +
  scale_fill_manual(values = c("#B2182B","#F4A582","#2166AC","#A6CEE3","#006600","#99CC33","#6A51A3","#9E9AC8")) +
  labs(x = "", y = "100/RMSE", fill = "Bulk RNAseq\nNormalization") +
  theme_bw() +
  theme(legend.position = "right",
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.y = element_text(color = "black", size = 13),
        axis.text.x = element_text(color = "black", size = 13, vjust = -1),
        text = element_text(size = 11),
        legend.text = element_text(color = "black", size = 11),
        legend.title = element_text(color = "black", size = 13),
        legend.spacing.x = unit(0.2, "cm"),
        legend.key.height = unit(0.6, "cm"),legend.key.width = unit(0.6, "cm")) 
S_fig

#readcount and RMSE---------------------------------------------------------------------------------------------------
readcount <- data.frame(reads = c(484578,4771431,11641601,26415269,724815,6799343,15417500,30508074,718100,7154899,17779601,42080465),
                        datasets = c(rep("NC",4),rep("GSE131928",4),rep("GSE182109",4)),
                        cellnumber = c(100,1000,2500,6000,100,1000,2500,6000,100,1000,2500,6000))
fileinuse2 <- fileinuse[,c("SumUp_NC","SumUp_GSE131928","SumUp_GSE182109","ReverseRMSE100_NC","ReverseRMSE100_GSE131928","ReverseRMSE100_GSE182109","BulkCellNumber")] %>% as.data.frame()
fileinuse2 <- fileinuse2 %>% group_by(BulkCellNumber) %>% summarise_all(mean) %>% as.data.frame()
readcount$RMSE100 <- append(append(fileinuse2$ReverseRMSE100_NC,fileinuse2$ReverseRMSE100_GSE131928),fileinuse2$ReverseRMSE100_GSE182109)
readcount$datasets <- factor(readcount$datasets, levels = c("NC","GSE131928","GSE182109"))
p1 <- ggplot(readcount, aes(x = reads, y = RMSE100, group = datasets, color = datasets)) +
  geom_line(size = 1) +
  geom_point(size = 5) +
  scale_color_manual(values=c("#B2182B","#F4A582","#2166AC")) +
  #  geom_text_repel() +
  labs(y="100/RMSE", x = "Reads", color = "Datasets") +
  theme_bw() +
  theme(legend.position = "right",
        axis.line = element_line(color = "black", size = 0.5),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 14),
        axis.text.x = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 14),
        text = element_text(size = 12),
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_text(color = "black", size = 14),
        legend.spacing.x = unit(0.2, "cm"),
        legend.key.height = unit(0.6, "cm"),legend.key.width = unit(0.6, "cm")) 

p1
#Comparison between celltypes RMSE----------------------------------------------------------------------
BulkNorm_Meth <- c("TPM","TMM","vst","NormCount","RawCount","ntd","logTPM","RPKM")

for (i in 1:length(BulkNorm_Meth)) {
  for (k in c("RMSE","Pearson")) {
    oridata <- read_excel(paste0("~/SumTable_BulkRNA_",k,"_AllDataset_Celltype_",BulkNorm_Meth[i],".xlsx")) %>% as.data.frame()
    print(paste0("open file--~/SumTable_BulkRNA_",k,"_AllDataset_Celltype_",BulkNorm_Meth[i],".xlsx"))
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
    if(k=="RMSE"){
      tempdf <- rbind(tempdf_all,tempdf_dwmd) %>% as.data.frame()
      print("Done RMSE")
    }else{
      tempdf <- rbind(tempdf,tempdf_all) %>% as.data.frame()
      tempdf <- rbind(tempdf,tempdf_dwmd) %>% as.data.frame()
      print("Done PearsonR")
    }
    remove(oridata,tempdf_all,tempdf_dwmd,a,b)
  }
  colnames(tempdf)[c(49:54)] <- c("DeconMeth","Min.pct","Sampling","GeneNumber","RefNorm","BulkRNA")
  tempdf$BulkNorm <- BulkNorm_Meth[i]
  tempdf$BulkCellNumber <- c(rep(100,135),rep(1000,135),rep(2500,135),rep(6000,135),
                             rep(100,135),rep(1000,135),rep(2500,135),rep(6000,135),
                             rep(100,135),rep(1000,135),rep(2500,135),rep(6000,135),
                             rep(100,135),rep(1000,135),rep(2500,135),rep(6000,135))
  tempdf$Value <- c(rep("RMSE",1080),rep("PearsonR",1080))
  if(i==1){
    mergetable <- tempdf
  }else{
    mergetable <- rbind(mergetable,tempdf) %>% as.data.frame()
  }
  print(paste0("Mergefile--",BulkNorm_Meth[i]))
  remove(tempdf)
}
mergetable$Min.pct <- str_replace(mergetable$Min.pct,pattern = "211", replacement = "0.2")
mergetable$Min.pct <- str_replace(mergetable$Min.pct,pattern = "411", replacement = "0.4")
mergetable$Min.pct <- str_replace(mergetable$Min.pct,pattern = "611", replacement = "0.6")
mergetable$ReverseRMSE100_MergeDataset <- 100/mergetable$SumUp_MergeDataset
mergetable$ReverseRMSE100_NC <- 100/mergetable$SumUp_NC
mergetable$ReverseRMSE100_GSE131928 <- 100/mergetable$SumUp_GSE131928
mergetable$ReverseRMSE100_GSE182109 <- 100/mergetable$SumUp_GSE182109

fileinuse <- mergetable[which(mergetable$Value=="RMSE"),] %>% as.data.frame()

fileinuse$BulkNorm <- factor(fileinuse$BulkNorm, levels = c("RawCount","vst","ntd","NormCount","TMM","TPM","logTPM","RPKM"))
fileinuse$DeconMeth <- factor(fileinuse$DeconMeth, levels = c("Cibersort","EPIC","ConsesnusTME"))
fileinuse$Min.pct <- factor(fileinuse$Min.pct, levels = c("0.2","0.4","0.6"))
fileinuse$Sampling <- factor(fileinuse$Sampling, levels = c("ALL","MinDist"))
fileinuse$GeneNumber <- factor(fileinuse$GeneNumber, levels = c(20,50,100))
fileinuse$RefNorm <- factor(fileinuse$RefNorm, levels = c("RawCount","LogNormalize","SCT","TPM","TMM"))

fileinuse2 <- fileinuse[-which(fileinuse$BulkCellNumber=="100"),] %>% as.data.frame()

library("RColorBrewer")
library(patchwork)
CellType <- c("Dendritic","Endothelial","Macrophage","Microglia","NKT","Oligodendrocytes","T cells","Tumor","B cells","Mural")
CellTypetitle <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells","Oligodendrocytes","T cells",
                   "Tumor cells","B cells","Mural cells")
#"DeconMeth","Min.pct","Sampling","GeneNumber","RefNorm","BulkNorm","BulkCellNumber"

S_fig <- list()
BulkNorm <- "RawCount"  #"RawCount","vst","ntd","NormCount","TMM","TPM","logTPM","RPKM"
for (i in c(1:10)) {
  informdf <- fileinuse2[which(fileinuse2$BulkNorm==BulkNorm),c("DeconMeth","Min.pct","Sampling","GeneNumber","RefNorm","BulkRNA","BulkNorm","BulkCellNumber","Value")] %>% as.data.frame()
  tempdf <- fileinuse2[which(fileinuse2$BulkNorm==BulkNorm),grep(pattern = CellType[i], colnames(fileinuse2))] %>% as.data.frame()
  tempdf <- cbind(tempdf, 100/tempdf) %>% as.data.frame()
  tempdf <- cbind(tempdf,informdf) %>% as.data.frame()
  colnames(tempdf)[5:8] <- c("MergeDataset","NC","GSE131928","GSE182109")
  
  S_fig[[i]] <- ggplot(data = tempdf, aes(x = RefNorm, y = MergeDataset, fill = DeconMeth)) +
    geom_boxplot(position=position_dodge(0.9)) +
    scale_fill_manual(values=c("#B2182B","#F4A582","#2166AC")) +
    ggtitle(CellTypetitle[[i]]) +
    labs(x = "", y = "100/RMSE", fill = "Deconvolution\nMethod") +
    theme_bw() +
    theme(legend.position = "right",
          panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"),
          axis.text.y = element_text(color = "black"),
          axis.title.y = element_text(color = "black", size = 13),
          axis.text.x = element_blank(),
          text = element_text(size = 13),
          legend.text = element_text(color = "black", size = 11),
          legend.title = element_text(color = "black", size = 13),
          legend.spacing.x = unit(0.2, "cm"),
          legend.key.height = unit(0.6, "cm"),legend.key.width = unit(0.6, "cm"),
          plot.title = element_text(hjust = 0.5, size = 13))
}
combined <- S_fig[[1]] + S_fig[[2]] + S_fig[[3]] + S_fig[[4]] + S_fig[[5]] + S_fig[[6]] + S_fig[[7]] + S_fig[[8]] +
  S_fig[[9]] + S_fig[[10]] & theme(legend.position = "right")
combined + plot_layout(ncol = 5, nrow = 2, guides = "collect")


#geom_point----------------------------------------------------------------------------------------------------------
library(readxl)
CellType <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells","Oligodendrocytes","T cells",
              "Tumor cells","B cells","Mural cells", "Sumup_merge")
Best_deconvolution_res <- read_excel("~/Best_deconvolution_res.xlsx") %>% as.data.frame()
Best_deconvolution_res <- Best_deconvolution_res[,c(1:8)] %>% as.data.frame()
colnames(Best_deconvolution_res) <- c("CellType","DeconMeth","Min.pct","Sampling","GeneNumber","RefNorm","BulkNorm","RMSE")
Best_deconvolution_res$RefNorm <- factor(Best_deconvolution_res$RefNorm, levels = c("RawCount","LogNormalize","SCT","TPM","TMM"))
Best_deconvolution_res$BulkNorm <- factor(Best_deconvolution_res$BulkNorm, levels = c("RawCount","vst","ntd","NormCount","TMM","TPM","logTPM","RPKM"))
Best_deconvolution_res$CellType <- factor(Best_deconvolution_res$CellType,
                                          levels = c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells",
                                                     "Oligodendrocytes","T cells","Tumor cells","B cells","Mural cells", "Sumup_merge","Best100",
                                                     "Best1000","Best2500","Best6000"))
Best_deconvolution_res$DeconMeth <- factor(Best_deconvolution_res$DeconMeth, levels = c("Cibersort","EPIC","ConsesnusTME"))

tempdf <- Best_deconvolution_res[which(Best_deconvolution_res$CellType %in% CellType),] %>% as.data.frame()
tempdf$RMSE100 <- 100/tempdf$RMSE
tempdf$CellType <- factor(tempdf$CellType,
                          levels = rev(c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells",
                                         "Oligodendrocytes","T cells","Tumor cells","B cells","Mural cells", "Sumup_merge")))

ggplot(data = tempdf, aes(x = BulkNorm, y = CellType)) +
  geom_point(aes(shape=DeconMeth, color=RefNorm, size=RMSE100)) +
  scale_color_manual(values = c("#B2182B","#F4A582","#2166AC","#A6CEE3","#006600")) +
  labs(size = "100/RMSE", color = "Reference\nNormalization", shape = "Deconvolution\nMethod") +
  theme_bw() + 
  theme(legend.position = "right",
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"),
        axis.text.y = element_text(color = "black", size = 13),
        axis.title.y = element_text(color = "black", size = 14),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title.x = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 13),
        legend.title = element_text(color = "black", size = 14),
        legend.spacing.x = unit(0.2, "cm"),
        legend.key.height = unit(0.5, "cm"),legend.key.width = unit(0.5, "cm")) +
  scale_x_discrete(guide = guide_axis(angle = 90))

#Best
Best_deconvolution_res <- Best_deconvolution_res[c(grep(pattern = "Best", Best_deconvolution_res$CellType),grep(pattern = "Single", Best_deconvolution_res$CellType)),c("CellType","BulkNorm","RMSE")] %>% as.data.frame()
Best_deconvolution_res$Method <- c(rep("Model",32),rep("Single",32))
Best_deconvolution_res <- Best_deconvolution_res[order(Best_deconvolution_res$CellType),] %>% as.data.frame()
Best_deconvolution_res$Pseudo_BulkRNAseq <- rep(c(rep(1,8),rep(10,8),rep(25,8),rep(60,8)),2)
Best_deconvolution_res$RMSE100 <- 100/Best_deconvolution_res$RMSE
p1 <- ggplot(Best_deconvolution_res, aes(x = Pseudo_BulkRNAseq, y = RMSE100, color = Method)) +
  geom_line(size = 0.5) +
  geom_point(size = 2.5) +
  scale_color_manual(values=c("#B2182B","#2166AC")) +
  facet_wrap(vars(BulkNorm), nrow = 2) +
  labs(y="100/RMSE", x = "Cell Number(10^2)", color = "Deconvolution Method") +
  xlim(-5,65) +
  ylim(10,40) +
  theme_bw() +
  theme(legend.position = "top",
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 14),
        axis.text.x = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 14),
        text = element_text(size = 12),
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_text(color = "black", size = 14),
        legend.spacing.x = unit(0.2, "cm"),
        legend.key.height = unit(0.6, "cm"),legend.key.width = unit(0.6, "cm"),
        strip.text.x = element_text(size = 14)) 

p1

#True proportions and estimated proportions--------------------------------------------------------------------------------------------
True_res <- read_csv("~/MergeDataset_BulkRNA_Celltype_N2500.csv",
                     show_col_types = FALSE) %>% as.data.frame()
True_res <- True_res[-which(True_res$Celltype=="Others"),] %>% as.data.frame()
row.names(True_res) <- True_res$Celltype
True_res$...1 <- NULL
True_res$Celltype <- NULL
True_res[is.na(True_res)==TRUE] <- 0
True_res[nrow(True_res) + 1,] <- 0
row.names(True_res)[nrow(True_res)] <- "otherCells"
Decon_res <- read_csv("~/DeconRes_ref_Best_mix_2500_RawCount_CibersortEPIC_noPro.csv",
                      show_col_types = FALSE) %>% as.data.frame()
row.names(Decon_res) <- Decon_res$...1
Decon_res$...1 <- NULL
Decon_res$Length <- NULL
target <- row.names(True_res)
Decon_res <- Decon_res[match(target, row.names(Decon_res)),] %>% as.data.frame()
row.names(Decon_res) <- row.names(True_res)
Decon_res[is.na(Decon_res)==TRUE] <- 0

tempdf <- True_res
tempdf <- tempdf %>% t() %>% as.data.frame()
tempdf$Dataset[grep(pattern = "cGBM", row.names(tempdf))] <- "GSE131928"
tempdf$Dataset[grep(pattern = "^GBM", row.names(tempdf))] <- "GSE182109"
tempdf$Dataset[grep(pattern = "NC", row.names(tempdf))] <- "NC"
tempdf$Sample_ID <- row.names(tempdf)

for (i in 1:(ncol(tempdf)-3)) {
  if(i==1){
    temptemp <- tempdf[,c(i,ncol(tempdf)-1,ncol(tempdf))] %>% as.data.frame()
    colnames(temptemp) <- c("Perc","Dataset","Sample_ID")
    Celltype <- rep(colnames(tempdf)[i],nrow(tempdf))
  }else{
    a <- tempdf[,c(i,ncol(tempdf)-1,ncol(tempdf))] %>% as.data.frame()
    colnames(a) <- c("Perc","Dataset","Sample_ID")
    temptemp <- rbind(temptemp,a) %>% as.data.frame()
    Celltype <- append(Celltype,rep(colnames(tempdf)[i],nrow(tempdf)))
    remove(a)
  }
}
temptemp$Celltype <- Celltype
row.names(temptemp) <- NULL
#temptemp_estimate <- temptemp
temptemp_true <- temptemp
remove(tempdf,Celltype,i,temptemp)

factorchosen <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells","Oligodendrocytes",
                  "T cells","Tumor cells","B cells","Mural cells")

temp_true_GSE182109 <- temptemp_true[which(temptemp_true$Dataset=="GSE182109"),] %>% as.data.frame()
Sample_order_GSE182109 <- paste0("GBM",seq(from = 1, to = 1008, by =1))
temp_estm_GSE182109 <- temptemp_estimate[which(temptemp_estimate$Dataset=="GSE182109"),] %>% as.data.frame()
temp_true_2data <- temptemp_true[-which(temptemp_true$Dataset=="GSE182109"),] %>% as.data.frame()
Sample_order_2data <- c(paste0("cGBM",seq(from = 1, to = 504, by =1)),paste0("NC",seq(from = 1, to = 504, by =1)))
temp_estm_2data <- temptemp_estimate[-which(temptemp_estimate$Dataset=="GSE182109"),] %>% as.data.frame()

#True-GSE182109, estm-GSE182109, true-2data, estimate-2data
p2 <- ggplot(temp_estm_GSE182109, aes(x = factor(Sample_ID,levels = Sample_order_GSE182109), y = Perc, fill = Celltype)) +
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
  scale_y_continuous(sec.axis = sec_axis(~., name="Percentage")) +
  theme_minimal() +
  guides(fill=guide_legend(nrow=2,byrow=FALSE)) +
  theme(legend.position = "top", #bottom 底部
        axis.line.y = element_line(color = "black", size = 1),
        axis.text.y = element_text(color = "black",size = 11),
        axis.title.y = element_blank(),
        axis.line.x = element_line(color = "white"),
        axis.text.x = element_blank(), 
        legend.text = element_text(color = "black", size = 13),
        legend.title = element_text(color = "black", size = 16),
        legend.margin = margin(),
        legend.spacing.x = unit(0.25, "cm"),
        legend.key.height = unit(0.25, "cm"),legend.key.width = unit(0.4, "cm")) 

ggarrange(p1,p2,p3,p4, ncol = 1, nrow = 4, common.legend = TRUE)

#RawCount only Comparison between celltypes RMSE/Pearson----------------------------------------------------------------------
#From Comparison between celltypes RMSE
#mergetable should be made first

RMSE <- read_excel(paste0("~/SumTable_BulkRNA_RMSE_AllDataset_Celltype_RawCount.xlsx")) %>% as.data.frame()
row.names(RMSE) <- RMSE$...1
RMSE$...1 <- NULL
Pear <- read_excel(paste0("~/SumTable_BulkRNA_Pearson_AllDataset_Celltype_RawCount.xlsx")) %>% as.data.frame()
row.names(Pear) <- Pear$...1
Pear$...1 <- NULL

best_df <- RMSE[grep(pattern = "Best", row.names(RMSE)),-c(49:52)] %>% as.data.frame()
best_df <- rbind(best_df,Pear[grep(pattern = "Best", row.names(Pear)),-c(49:52)] %>% as.data.frame()) %>% as.data.frame()

best_df$DeconMeth <- "OptGBM"
best_df$Min.pct <- "Mix"
best_df$Sampling <- "Mix"
best_df$GeneNumber <- "Mix"
best_df$RefNorm <- "Mix"
best_df$BulkRNA <- rep(c("RawCount100","RawCount1000","RawCount2500","RawCount6000"),2)    #change
best_df$BulkNorm <- "RawCount"    #change
best_df$BulkCellNumber <- rep(c(100,1000,2500,6000))
best_df$Value <- c(rep("RMSE",4),rep("Pearson",4))
mergetable <- rbind(mergetable,best_df) %>% as.data.frame()

mergetable$Min.pct <- str_replace(mergetable$Min.pct,pattern = "211", replacement = "0.2")
mergetable$Min.pct <- str_replace(mergetable$Min.pct,pattern = "411", replacement = "0.4")
mergetable$Min.pct <- str_replace(mergetable$Min.pct,pattern = "611", replacement = "0.6")
mergetable$ReverseRMSE100_MergeDataset <- 100/mergetable$SumUp_MergeDataset
mergetable$ReverseRMSE100_NC <- 100/mergetable$SumUp_NC
mergetable$ReverseRMSE100_GSE131928 <- 100/mergetable$SumUp_GSE131928
mergetable$ReverseRMSE100_GSE182109 <- 100/mergetable$SumUp_GSE182109

#RawCount only single and Optimal-GBM RMSE curves
fileinuse <- mergetable[which(mergetable$Value=="RMSE" & mergetable$BulkNorm=="RawCount"),] %>% as.data.frame()
fileinuse$DeconMeth_new <- "Others"
fileinuse$DeconMeth_new[which(fileinuse$DeconMeth=="OptGBM")] <- "OptGBM"

readcount <- data.frame(Reads = c(661398,6470143,15654576,35271068,484578,4771431,11641601,26415269,724815,6799343,15417500,30508074,718100,7154899,17779601,42080465),
                        Datasets = c(rep("Merge",4),rep("NC",4),rep("GSE131928",4),rep("GSE182109",4)),
                        BulkCellNumber = c(100,1000,2500,6000,100,1000,2500,6000,100,1000,2500,6000,100,1000,2500,6000))
readcount <- rbind(readcount,readcount) %>% as.data.frame()
readcount$DeconMeth_new <- c(rep("Others",16),rep("OptGBM",16))

fileinuse2 <- fileinuse[,c("ReverseRMSE100_MergeDataset","ReverseRMSE100_NC","ReverseRMSE100_GSE131928","ReverseRMSE100_GSE182109","BulkCellNumber","DeconMeth_new")] %>% as.data.frame()
fileinuse2 <- fileinuse2 %>% group_by(BulkCellNumber, DeconMeth_new) %>% summarise_all(mean) %>% as.data.frame()
for (i in 1:4) {
  if(i==1){
    fileinuse2_t <- fileinuse2[,c(1,2,(2+i))] %>% as.data.frame()
    colnames(fileinuse2_t) <- c("BulkCellNumber","DeconMeth_new","RMSE100")
    fileinuse2_t$Datasets <- str_split_fixed(colnames(fileinuse2)[2+i], pattern = "_",2)[,2]
  }else{
    a <- fileinuse2[,c(1,2,(2+i))] %>% as.data.frame()
    colnames(a) <- c("BulkCellNumber","DeconMeth_new","RMSE100")
    a$Datasets <- str_split_fixed(colnames(fileinuse2)[2+i], pattern = "_",2)[,2]
    fileinuse2_t <- rbind(fileinuse2_t,a) %>% as.data.frame()
    remove(a)
  }
}
fileinuse2_t$Datasets <- str_replace(fileinuse2_t$Datasets, pattern = "MergeDataset", replacement = "Merge")
readcount2 <- merge(x=readcount,y=fileinuse2_t,
                    by.x=c("BulkCellNumber","DeconMeth_new","Datasets"),
                    by.y=c("BulkCellNumber","DeconMeth_new","Datasets"),
                    all = TRUE) %>% as.data.frame()
readcount2$DataMeth <- paste(readcount2$Datasets,readcount2$DeconMeth_new, sep = "--")

p1 <- ggplot(readcount2[-which(readcount2$Datasets=="Merge"),], aes(x = Reads, y = RMSE100, group = DataMeth, color = DataMeth)) +
  geom_line(aes(linetype=DataMeth), size = 0.5) +
  geom_point(size = 2.5) +
  scale_color_manual(values=c("#B2182B","#F4A582","#2166AC","#A6CEE3","#006600","#99CC33","#6A51A3","#9E9AC8")) +
  scale_linetype_manual(values=c("solid", "dashed","solid", "dashed","solid", "dashed")) +
  labs(y="100/RMSE", x = "Reads", color = "Dataset--\nDeconvolution method") +
  theme_bw() + 
  theme(legend.position = "right",
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 14),
        axis.text.x = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_text(color = "black", size = 14)) 

p1
#RawCount only single and Optimal-GBM and each cell type RMSE/PearsonR--------------------------------------------------------------------------------------
#mergetable from above
fileinuse <- mergetable[which(mergetable$BulkNorm=="RawCount"),
                        c(grep(pattern = "MergeDataset",colnames(mergetable)),49:57)] %>% as.data.frame()
colnames(fileinuse) <- str_remove(colnames(fileinuse),pattern = "_MergeDataset")
fileinuse <- fileinuse[,-which(colnames(fileinuse) %in% c("otherCells","ReverseRMSE100","BulkRNA"))] %>% as.data.frame()

OptGBM_temp <- fileinuse[which(fileinuse$DeconMeth=="OptGBM"),] %>% as.data.frame()
Others_temp_RMSE <- fileinuse[-which(fileinuse$DeconMeth=="OptGBM"),] %>% as.data.frame()
Others_temp_RMSE <- Others_temp_RMSE[which(Others_temp_RMSE$Value=="RMSE"),] %>% as.data.frame()
Others_temp_Pearson <- fileinuse[-which(fileinuse$DeconMeth=="OptGBM"),] %>% as.data.frame()
Others_temp_Pearson <- Others_temp_Pearson[which(Others_temp_Pearson$Value=="PearsonR"),] %>% as.data.frame()

target_by <- c("DeconMeth","Min.pct","Sampling","GeneNumber","RefNorm","BulkNorm","BulkCellNumber")
RMSE_mean <- c()
PearsonR_mean <- c()
RMSE100_mean <- c()
Celltype <- c()
Method <- c()
for (i in 1:11) {
  tempdf <- merge(x=Others_temp_RMSE[,c(i,12:18)],y=Others_temp_Pearson[,c(i,12:18)], by.x=target_by, by.y=target_by, all=TRUE)
  colnames(tempdf)[8:9] <- c("RMSE","PearsonR")
  tempdf$RMSE100 <- 100/tempdf$RMSE
  for (j in c(100,1000,2500,6000)) {
    a <- tempdf[which(tempdf$BulkCellNumber==j),]
    a$Celltype <- colnames(Others_temp_RMSE)[i]
    a$Method <- "Best"
    a <- a[order(a$RMSE),]
    if(i==1 & j==100){
      best_df <- a[1,c(6:12)] %>% as.data.frame()
    }else{
      best_df <- rbind(best_df,a[1,c(6:12)] %>% as.data.frame()) %>% as.data.frame()
    }
    RMSE_mean <- append(RMSE_mean,mean(a$RMSE))
    PearsonR_mean <- append(PearsonR_mean,mean(a$PearsonR))
    RMSE100_mean <- append(RMSE100_mean,mean(a$RMSE100))
    Celltype <- append(Celltype,colnames(Others_temp_RMSE)[i])
    Method <- append(Method,"Mean")
    remove(a)
  }
  remove(tempdf)
}
mean_df <- data.frame(BulkNorm = rep("RawCount",44),
                      BulkCellNumber = rep(c(100,1000,2500,6000),11),
                      RMSE = RMSE_mean,
                      PearsonR = PearsonR_mean,
                      RMSE100 = RMSE100_mean,
                      Celltype = Celltype,
                      Method = Method)

a_index <- OptGBM_temp[c(1:4),c("BulkNorm","BulkCellNumber","DeconMeth")] %>% as.data.frame()
row.names(a_index) <- NULL
for (i in 1:11) {
  a <- cbind(OptGBM_temp[c(1:4),i],OptGBM_temp[c(5:8),i]) %>% as.data.frame()
  colnames(a) <- c("RMSE","PearsonR")
  a$RMSE100 <- 100/a$RMSE
  a$Celltype <- rep(colnames(OptGBM_temp)[i],4)
  a <- cbind(a,a_index) %>% as.data.frame()
  if(i==1){
    OptGBM_df <- a
  }else{
    OptGBM_df <- rbind(OptGBM_df,a) %>% as.data.frame()
  }
  remove(a)
}
colnames(OptGBM_df)[which(colnames(OptGBM_df)=="DeconMeth")] <- "Method"

tempdf <- merge(x=Others_temp_RMSE[,c(1:11,12:18)],y=Others_temp_Pearson[,c(1:11,12:18)], by.x=target_by, by.y=target_by, all=TRUE)
for (i in c(100,1000,2500,6000)) {
  a <- tempdf[which(tempdf$BulkCellNumber==i),]
  a <- a[order(a$SumUp.x),]
  if(i==100){
    singlebest_df <- a[1,c(6:29)] %>% as.data.frame()
  }else{
    singlebest_df <- rbind(singlebest_df,a[1,c(6:29)] %>% as.data.frame()) %>% as.data.frame()
  }
  remove(a)
}
remove(tempdf)
for (i in 1:11) {
  a <- singlebest_df[,c(1,2,(2+i),(13+i))] %>% as.data.frame()
  colnames(a)[3:4] <- c("RMSE","PearsonR")
  a$RMSE100 <- 100/a$RMSE
  a$Celltype <- str_remove(colnames(singlebest_df)[2+i], pattern = ".x")
  a$Method <- "OneRefBest"
  if(i==1){
    oneref_df <- a
  }else{
    oneref_df <- rbind(oneref_df,a) %>% as.data.frame()
  }
  remove(a)
}

target_by <- colnames(best_df)
mean_df <- mean_df[,sort(target_by)] %>% as.data.frame()
OptGBM_df <- OptGBM_df[,sort(target_by)] %>% as.data.frame()
oneref_df <- oneref_df[,sort(target_by)] %>% as.data.frame()
total_df <- rbind(best_df,mean_df) %>% as.data.frame()
total_df <- rbind(total_df,OptGBM_df) %>% as.data.frame()
total_df <- rbind(total_df,oneref_df) %>% as.data.frame()
total_df$group <- paste(total_df$BulkCellNumber,total_df$Method,sep = "_")
group_order <- sort(unique(total_df$group), decreasing = FALSE)
celltype_order <- rev(c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells",
                        "Oligodendrocytes","T cells","Tumor cells","B cells","Mural cells", "SumUp"))

ggplot(data = total_df, aes(x = factor(group, levels = group_order), y = factor(Celltype, levels = celltype_order))) +
  geom_point(aes(color=PearsonR, size=RMSE100)) +
  scale_color_gradient2(midpoint=0, low="#2166AC", mid="#CCCCCC", high="#B2182B", limits = c(-1,1)) +
  labs(size = "100/RMSE", color = "Pearson R") +
  theme_bw() + 
  theme(legend.position = "right",
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"),
        axis.text.y = element_text(color = "black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        #  axis.text.x = element_text(color = "black", size = 13),
        axis.title.x = element_blank(),
        legend.text = element_text(color = "black", size = 13),
        legend.title = element_text(color = "black", size = 14),
        legend.spacing.x = unit(0.2, "cm"),
        legend.key.height = unit(0.5, "cm"),legend.key.width = unit(0.5, "cm")) +
  scale_x_discrete(guide = guide_axis(angle = 90))

#N2500 for example，differences of RMSE100 and pearsonR in BulkNorm and Celltype-----------------------------------------------------------------------------
BulkNorm_order <- c("RawCount","vst","ntd","NormCount","TMM","TPM","logTPM","RPKM")

for (i in 1:8) {
  RMSE <- read_excel(paste0("~/SumTable_BulkRNA_RMSE_AllDataset_Celltype_",BulkNorm_order[i],".xlsx")) %>% as.data.frame()
  row.names(RMSE) <- RMSE$...1
  RMSE$...1 <- NULL
  RMSE <- RMSE[grep(pattern = "Best",row.names(RMSE)),c(1,5,9,13,17,21,25,29,33,37,45)] %>% as.data.frame()
  RMSE <- RMSE[grep(pattern = "2500",row.names(RMSE)),] %>% as.data.frame()
  colnames(RMSE) <- str_remove(colnames(RMSE), pattern = "_MergeDataset")
  row.names(RMSE) <- NULL
  Pear <- read_excel(paste0("~/SumTable_BulkRNA_Pearson_AllDataset_Celltype_",BulkNorm_order[i],".xlsx")) %>% as.data.frame()
  row.names(Pear) <- Pear$...1
  Pear$...1 <- NULL
  Pear <- Pear[grep(pattern = "Best",row.names(Pear)),c(1,5,9,13,17,21,25,29,33,37,45)] %>% as.data.frame()
  Pear <- Pear[grep(pattern = "2500",row.names(Pear)),] %>% as.data.frame()
  colnames(Pear) <- str_remove(colnames(Pear), pattern = "_MergeDataset")
  row.names(Pear) <- NULL
  for(j in 1:11){
    temp_df <- cbind(RMSE[,j],Pear[,j]) %>% as.data.frame()
    colnames(temp_df) <- c("RMSE","PearsonR")
    temp_df$Celltype <- colnames(RMSE)[j]
    temp_df$BulkNorm <- BulkNorm_order[i]
    if(i==1 & j==1){
      combinetable <- temp_df
    }else{
      combinetable <- rbind(combinetable,temp_df) %>% as.data.frame()
    }
    remove(temp_df)
    print(paste0("Done--",BulkNorm_order[i],"--",colnames(RMSE)[j]))
  }
  remove(RMSE,Pear)
  print(paste0("Done--",BulkNorm_order[i]))
}
combinetable$RMSE100 <- 100/combinetable$RMSE

celltype_order <- rev(c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells",
                        "Oligodendrocytes","T cells","Tumor cells","B cells","Mural cells", "SumUp"))
ggplot(data = combinetable, aes(x = factor(BulkNorm, levels = BulkNorm_order), y = factor(Celltype, levels = celltype_order))) +
  geom_point(aes(color=PearsonR, size=RMSE100)) +
  scale_color_gradient2(midpoint=0, low="#2166AC", mid="#CCCCCC", high="#B2182B", limits = c(-1,1)) +
  labs(size = "100/RMSE", color = "Pearson R") +
  theme_bw() + 
  theme(legend.position = "right",
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"),
        axis.text.y = element_text(color = "black", size = 13),
        axis.title.y = element_blank(),
        #  axis.text.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 13),
        axis.title.x = element_blank(),
        legend.text = element_text(color = "black", size = 13),
        legend.title = element_text(color = "black", size = 14),
        legend.spacing.x = unit(0.2, "cm"),
        legend.key.height = unit(0.5, "cm"),legend.key.width = unit(0.5, "cm")) +
  scale_x_discrete(guide = guide_axis(angle = 90))

#N2500_celltype_Optimal-GBM_x-y correlation plot-------------------------------------------------------------------------
True_res <- read_csv(paste0("~/MergeDataset_BulkRNA_Celltype_N2500.csv"),
                     show_col_types = FALSE) %>% as.data.frame()
True_res <- True_res[-which(True_res$Celltype=="Others"),] %>% as.data.frame()
row.names(True_res) <- True_res$Celltype
True_res$...1 <- NULL
True_res$Celltype <- NULL
True_res[is.na(True_res)==TRUE] <- 0
True_res[nrow(True_res) + 1,] <- 0
row.names(True_res)[nrow(True_res)] <- "otherCells"   #EPIC: otherCells
True_res <- True_res %>% t() %>% as.data.frame()

Decon_res <- read_csv(paste0("~/DeconRes_ref_Best_mix_2500_RawCount_CibersortEPIC_noPro.csv"),
                      show_col_types = FALSE) %>% as.data.frame()
row.names(Decon_res) <- Decon_res$...1
Decon_res$...1 <- NULL
Decon_res <- Decon_res %>% t() %>% as.data.frame()
target <- row.names(True_res)
Decon_res <- Decon_res[match(target, row.names(Decon_res)),] %>% as.data.frame()
target <- colnames(True_res)
Decon_res <- Decon_res[,match(target, colnames(Decon_res))] %>% as.data.frame()
Decon_res[is.na(Decon_res)==TRUE] <- 0

for (i in 1:10) {
  temp_df <- cbind(True_res[,i], Decon_res[,i]) %>% as.data.frame()
  colnames(temp_df) <- c("True","Estimate")
  temp_df$Celltype <- colnames(True_res)[i]
  temp_df$Sample <- row.names(True_res)
  temp_df$Dataset <- "GSE182109"
  temp_df$Dataset[grep(pattern = "NC", temp_df$Sample)] <- "NC"
  temp_df$Dataset[grep(pattern = "cGBM", temp_df$Sample)] <- "GSE131928"
  if(i==1){
    combinetable <- temp_df
  }else{
    combinetable <- rbind(combinetable,temp_df) %>% as.data.frame()
  }
  remove(temp_df)
}
combinetable$Celltype <- factor(combinetable$Celltype, levels = c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells",
                                                                  "Oligodendrocytes","T cells","Tumor cells","B cells","Mural cells"))
combinetable$Dataset <- factor(combinetable$Dataset, levels = c("GSE131928","GSE182109","NC"))
celltype <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells","Oligodendrocytes","T cells",
              "Tumor cells","B cells","Mural cells")

cell.labs <- c("Dendritic cells","Endothelial cells","Macrophage-like\nGAMs","Microglia-like\nGAMs","NKT-like cells","Oligodendrocytes","T cells",
               "Tumor cells","B cells","Mural cells")
names(cell.labs) <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells","Oligodendrocytes","T cells",
                      "Tumor cells","B cells","Mural cells")
p1 <- ggplot(combinetable, aes(x = True, y = Estimate, color = factor(paste0(Dataset,"    ")))) +
  geom_point(size = 0.1) +
  scale_color_manual(values=c("#B2182B","#F4A582","#2166AC")) +
  geom_smooth(method=lm, se=FALSE) +
  facet_wrap(vars(Celltype), nrow = 2, scales="free",
             labeller = labeller(Celltype = cell.labs)) +
  labs(y="Estimated proportion (%)", x = "True proportion (%)", color = "Dataset") +
  theme_bw() +
  theme(legend.position = "top",
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 15),
        axis.title.x = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black", size = 18),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.1, "cm"),
        legend.key.width = unit(0.6, "cm"),
        strip.text.x = element_text(size = 18),
        strip.background = element_rect(color = "black", size = 1)) 

p1

i <- 6
a131928 <- combinetable[which(combinetable$Celltype==celltype[i] & combinetable$Dataset=="GSE131928"),]
a182109 <- combinetable[which(combinetable$Celltype==celltype[i] & combinetable$Dataset=="GSE182109"),]
aNC <- combinetable[which(combinetable$Celltype==celltype[i] & combinetable$Dataset=="NC"),]
ba131928 <- cor.test(a131928$True, a131928$Estimate)
ba182109 <- cor.test(a182109$True, a182109$Estimate)
baNC <- cor.test(aNC$True, aNC$Estimate)
ba131928$p.value
ba131928$estimate
ba182109$p.value
ba182109$estimate
baNC$p.value
baNC$estimate

#Correlation between cell type abundance-----------------------------------------------------------------------------------------------------------------------
library(corrplot)
library(ggrepel)
library(ggpubr)
library(gridExtra)

TCGA <- read_csv("~/TCGA_Survival_CellType_145RawCount_CombineTable.csv") %>% as.data.frame()
row.names(TCGA) <- TCGA$...1
TCGA$...1 <- NULL
TCGA <- TCGA[,c(1:10)] %>% as.data.frame()
TCGA <- TCGA[,order(colnames(TCGA))] %>% as.data.frame()
CGGA <- read_csv("~/CGGA_Survival_CellType_153_RawCount_CombineTable.csv") %>% as.data.frame()
row.names(CGGA) <- CGGA$...1
CGGA$...1 <- NULL
CGGA <- CGGA[,c(1:10)] %>% as.data.frame()
CGGA <- CGGA[,order(colnames(CGGA))] %>% as.data.frame()
CGMH <- read_csv("d~/CGMH_Survival_CellType_40GBM_RawCount_CombineTable.csv") %>% as.data.frame()
row.names(CGMH) <- CGMH$...1
CGMH$...1 <- NULL
CGMH <- CGMH[,c(1:10)] %>% as.data.frame()
CGMH <- CGMH[,order(colnames(CGMH))] %>% as.data.frame()
mergedata <- rbind(TCGA,CGGA,CGMH) %>% as.data.frame()

tempdf <- mergedata
if(length(which(colSums(tempdf)==0))!=0){
  tempdf <- tempdf[,-which(colSums(tempdf)==0)] %>% as.data.frame()
}
M <- cor(tempdf)
testRes <- cor.mtest(tempdf, conf.level = 0.95)
corrplot(M,
         p.mat = testRes$p,
         #insig='blank',
         sig.level = 0.05,
         addrect = 3,
         rect.col = 'black',
         rect.lwd = 3,
         pch.col = 'grey50',
         order = 'hclust',
         #type = 'lower',
         tl.col = 'black',
         tl.cex = 1.5,
         addCoef.col = 'black',
         number.cex = 0.8,
         #tl.srt = 45,
         #diag=FALSE,
         col = rev(COL2('RdBu', 200)))

#Newly diagnosis and recurrence---------------------------------------------------------------------------------------------------------------------
TCGA_GBM_recur <- read_csv("~/Raw_Counts_HGNC_TCGA_GBM_total_duplicateRemoved_tumor_NewlyRecur.csv") %>% as.data.frame()
row.names(TCGA_GBM_recur) <- TCGA_GBM_recur$...1
TCGA_GBM_recur$...1 <- NULL
colnames(TCGA_GBM_recur) <- str_sub(colnames(TCGA_GBM_recur), start = 1,end = 15)
#OptGBM-GBM deconvolution
OptGBM_table <- OptGBM_table %>% t() %>% as.data.frame()
OptGBM_table$TCGA_ID <- str_sub(row.names(OptGBM_table), start = 1, end = 12)
OptGBM_table$PRS_type <- "Newly"
OptGBM_table$PRS_type[grep("02", str_split_fixed(row.names(OptGBM_table),pattern = "-",4)[,4])] <- "Recurrent"
OptGBM_table$Endo_Tumor <- OptGBM_table$`Endothelial cells` + OptGBM_table$`Tumor cells`
OptGBM_table$Immune <- OptGBM_table$`B cells`+ OptGBM_table$`NKT-like cells` + OptGBM_table$`T cells`
OptGBM_table$GAM_Mural <- OptGBM_table$`Macrophage-like GAMs` + OptGBM_table$`Microglia-like GAMs` + OptGBM_table$`Mural cells`

celltype <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells","Oligodendrocytes","T cells",
              "Tumor cells","B cells","Mural cells")
for (i in 1:10) {
  temp_df <- OptGBM_table[,c("TCGA_ID","PRS_type",celltype[i])]
  colnames(temp_df) <- c("TCGA_ID","PRS_type","OptGBM")
  temp_df$Celltype <- celltype[i]
  if(i==1){
    combinetable <- temp_df
  }else{
    combinetable <- rbind(combinetable,temp_df) %>% as.data.frame()
  }
  remove(temp_df)
}
combinetable_TCGA <- combinetable

celltype2 <- c("Dendritic cells","Endo_Tumor","Immune","GAM_Mural","Oligodendrocytes")
for (i in 1:5) {
  temp_df <- OptGBM_table[,c("TCGA_ID","PRS_type",celltype2[i])]
  colnames(temp_df) <- c("TCGA_ID","PRS_type","OptGBM")
  temp_df$Celltype <- celltype2[i]
  if(i==1){
    combinetable <- temp_df
  }else{
    combinetable <- rbind(combinetable,temp_df) %>% as.data.frame()
  }
  remove(temp_df)
}
combinetable_TCGA_C5 <- combinetable

combinetable_TCGA$Celltype <- factor(combinetable_TCGA$Celltype, levels = c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells",
                                                                            "Oligodendrocytes","T cells","Tumor cells","B cells","Mural cells"))
combinetable_TCGA$TCGA_ID <- factor(combinetable_TCGA$TCGA_ID, levels = c("TCGA-06-0210","TCGA-06-0190","TCGA-06-0211","TCGA-19-4065","TCGA-14-1034","TCGA-06-0125"))
cell.labs <- c("Dendritic cells","Endothelial cells","Macrophage-like\nGAMs","Microglia-like\nGAMs","NKT-like cells","Oligodendrocytes","T cells",
               "Tumor cells","B cells","Mural cells")
names(cell.labs) <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells","Oligodendrocytes","T cells",
                      "Tumor cells","B cells","Mural cells")
blank_data <- data.frame(Celltype = factor(c("Dendritic cells","Dendritic cells","Endothelial cells","Endothelial cells","Macrophage-like GAMs",
                                             "Macrophage-like GAMs","Microglia-like GAMs","Microglia-like GAMs","NKT-like cells","NKT-like cells",
                                             "Oligodendrocytes","Oligodendrocytes","T cells","T cells","Tumor cells","Tumor cells","B cells","B cells",
                                             "Mural cells","Mural cells"), levels = c("Dendritic cells","Endothelial cells","Macrophage-like GAMs",
                                                                                      "Microglia-like GAMs","NKT-like cells",
                                                                                      "Oligodendrocytes","T cells","Tumor cells","B cells","Mural cells")),
                         TCGA_ID = c(rep("TCGA-06-0125",20)),
                         x = c("Newly","Recurrent","Newly","Recurrent","Newly","Recurrent","Newly","Recurrent","Newly","Recurrent",
                               "Newly","Recurrent","Newly","Recurrent","Newly","Recurrent","Newly","Recurrent","Newly","Recurrent"),
                         y = c(0, 4, 0, 30, 0, 35, 0, 30, -0.05, 0.05, 0, 25, 0, 0.6, 0, 90, 0, 0.015, 0, 85))

p1 <- ggplot(combinetable_TCGA, aes(x = PRS_type, y = OptGBM, group = TCGA_ID)) +
  geom_point(aes(color = factor(paste0(TCGA_ID,"           "))), size = 3) +
  geom_line(aes(color = factor(paste0(TCGA_ID,"           "))), size=1) +
  scale_color_manual(values=c("#B2182B","#F4A582","#2166AC","#A6CEE3","#006600","#99CC33")) +
  geom_blank(data = blank_data, aes(x = x, y = y, group = TCGA_ID)) +
  facet_wrap(vars(Celltype), nrow = 2, scales="free_y",
             labeller = labeller(Celltype = cell.labs)) +
  labs(y="OptGBM proportion (%)", x = "") +
  theme_bw() +
  theme(legend.position = "right",
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 18, vjust = 2),
        axis.text.x = element_text(color = "black", size = 18, angle = 40, vjust = 0.7, hjust = 0.7),
        axis.title.x = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black", size = 18),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.1, "cm"),
        legend.key.width = unit(1, "cm"),
        strip.text.x = element_text(size = 18),
        strip.background = element_rect(color = "black", size = 1)) +
  guides(color = guide_legend(byrow = TRUE))
p1

celltype2 <- c("Dendritic cells","Endo_Tumor","Immune","GAM_Mural","Oligodendrocytes")
combinetable_TCGA_C5$Celltype <- factor(combinetable_TCGA_C5$Celltype, levels = c("Dendritic cells","Endo_Tumor","Immune","GAM_Mural","Oligodendrocytes"))
combinetable_TCGA_C5$TCGA_ID <- factor(combinetable_TCGA_C5$TCGA_ID, levels = c("TCGA-06-0210","TCGA-06-0190","TCGA-06-0211","TCGA-19-4065","TCGA-14-1034","TCGA-06-0125"))
cell.labs <- c("Dendritic cells","Tumor/Endothelial\ncells","Immune cells","GAMs/Mural cells","Oligodendrocytes")
names(cell.labs) <- c("Dendritic cells","Endo_Tumor","Immune","GAM_Mural","Oligodendrocytes")
blank_data <- data.frame(Celltype = factor(c("Dendritic cells","Dendritic cells","Endo_Tumor","Endo_Tumor","Immune","Immune","GAM_Mural",
                                             "GAM_Mural","Oligodendrocytes","Oligodendrocytes"),
                                           levels = c("Dendritic cells","Endo_Tumor","Immune","GAM_Mural","Oligodendrocytes")),
                         TCGA_ID = c(rep("TCGA-06-0125",10)),
                         x = c("Newly","Recurrent","Newly","Recurrent","Newly","Recurrent","Newly","Recurrent","Newly","Recurrent"),
                         y = c(0, 4, 0, 110, 0, 0.5, 0, 110, 0, 20))
p1 <- ggplot(combinetable_TCGA_C5, aes(x = PRS_type, y = OptGBM, group = TCGA_ID)) +
  geom_point(aes(color = factor(paste0(TCGA_ID,"           "))), size = 3) +
  geom_line(aes(color = factor(paste0(TCGA_ID,"           "))), size=1) +
  scale_color_manual(values=c("#B2182B","#F4A582","#2166AC","#A6CEE3","#006600","#99CC33")) +
  geom_blank(data = blank_data, aes(x = x, y = y, group = TCGA_ID)) +
  facet_wrap(vars(Celltype), nrow = 1, scales="free",
             labeller = labeller(Celltype = cell.labs)) +
  labs(y="OptGBM proportion (%)", x = "") +
  theme_bw() +
  theme(legend.position = "top",
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 18, vjust = 2),
        axis.text.x = element_text(color = "black", size = 18, angle = 40, vjust = 0.7, hjust = 0.7),
        axis.title.x = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black", size = 18),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.1, "cm"),
        legend.key.width = unit(1, "cm"),
        strip.text.x = element_text(size = 18),
        strip.background = element_rect(color = "black", size = 1)) +
  guides(color = guide_legend(byrow = TRUE))
p1

#CGGA-----------------------------------------------
#建立recurrent data-----
library(readr)
CGGA693 <- read_delim("~/CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt",  #downloaded from CGGA website
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE) %>% as.data.frame()
row.names(CGGA693) <- CGGA693$gene_name
CGGA693$gene_name <- NULL
CGGA325 <- read_delim("~/CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt", #downloaded from CGGA website
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE) %>% as.data.frame()
row.names(CGGA325) <- CGGA325$gene_name
CGGA325$gene_name <- NULL
CGGA <- cbind(CGGA325,CGGA693) %>% as.data.frame()
write.csv(CGGA,file = "~/CGGA1018_RawReadCount.csv")

library(readxl)
library(readr)
Sample_Info <- read_excel("~/1018_CGGA_RNASeq_Sample_Info.xlsx") %>% as.data.frame() #downloaded from CGGA website
Sample_Info <- Sample_Info[which(Sample_Info$IDH_mutation_status=="Wildtype" & Sample_Info$Histology %in% c("GBM","rGBM")),] %>% as.data.frame()
CGGA <- read_csv("~/CGGA1018_RawReadCount.csv") %>% as.data.frame()
row.names(CGGA) <- CGGA$...1
CGGA$...1 <- NULL
CGGA <- CGGA[,which(colnames(CGGA) %in% Sample_Info$CGGA_ID)] %>% as.data.frame()
#Optimal-GBM deconvolution
OptGBM_table <- OptGBM_table %>% t() %>% as.data.frame()
OptGBM_table$CGGA_ID <- row.names(OptGBM_table)
write.csv(OptGBM_table, file = "~/OptGBM_table_CGGA_RawCount_Recur.csv")
#建立recurrent data-----

Sample_Info <- read_excel("~/1018_CGGA_RNASeq_Sample_Info.xlsx") %>% as.data.frame() #downloaded from CGGA website
Sample_Info <- Sample_Info[which(Sample_Info$IDH_mutation_status=="Wildtype" & Sample_Info$Histology %in% c("rGBM")),] %>% as.data.frame()
CGGA_Recur <- read_csv("~/OptGBM_table_CGGA_RawCount_Recur.csv") %>% as.data.frame()
row.names(CGGA_Recur) <- CGGA_Recur$...1
CGGA_Recur$...1 <- NULL
CGGA_Recur$CGGA_ID <- NULL
CGGA_Recur <- CGGA_Recur[which(row.names(CGGA_Recur) %in% Sample_Info$CGGA_ID),] %>% as.data.frame()
CGGA_Recur <- CGGA_Recur[,order(colnames(CGGA_Recur))] %>% as.data.frame()
CGGA_Recur$PRS_type <- "Recurrent"

CGGA_newly <- read_csv("~/CGGA_Survival_CellType_153_RawCount_CombineTable.csv") %>% as.data.frame()
row.names(CGGA_newly) <- CGGA_newly$...1
CGGA_newly$...1 <- NULL
CGGA_newly <- CGGA_newly[,c(1:10)] %>% as.data.frame()
CGGA_newly <- CGGA_newly[,order(colnames(CGGA_newly))] %>% as.data.frame()
CGGA_newly$PRS_type <- "Newly"

CGGA <- rbind(CGGA_newly,CGGA_Recur) %>% as.data.frame()
CGGA$CGGA_ID　<- row.names(CGGA)

celltype <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells","Oligodendrocytes","T cells",
              "Tumor cells","B cells","Mural cells")
for (i in 1:10) {
  temp_df <- CGGA[,c("CGGA_ID","PRS_type",celltype[i])]
  colnames(temp_df) <- c("CGGA_ID","PRS_type","OptGBM")
  temp_df$Celltype <- celltype[i]
  if(i==1){
    combinetable <- temp_df
  }else{
    combinetable <- rbind(combinetable,temp_df) %>% as.data.frame()
  }
  remove(temp_df)
}

combinetable$Celltype <- factor(combinetable$Celltype, levels = c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells",
                                                                  "Oligodendrocytes","T cells","Tumor cells","B cells","Mural cells"))
cell.labs <- c("Dendritic cells","Endothelial cells","Macrophage-like\nGAMs","Microglia-like\nGAMs","NKT-like cells","Oligodendrocytes","T cells",
               "Tumor cells","B cells","Mural cells")
names(cell.labs) <- c("Dendritic cells","Endothelial cells","Macrophage-like GAMs","Microglia-like GAMs","NKT-like cells","Oligodendrocytes","T cells",
                      "Tumor cells","B cells","Mural cells")
blank_data <- data.frame(Celltype = factor(c("Dendritic cells","Dendritic cells","Endothelial cells","Endothelial cells","Macrophage-like GAMs",
                                             "Macrophage-like GAMs","Microglia-like GAMs","Microglia-like GAMs","NKT-like cells","NKT-like cells",
                                             "Oligodendrocytes","Oligodendrocytes","T cells","T cells","Tumor cells","Tumor cells","B cells","B cells",
                                             "Mural cells","Mural cells"), levels = c("Dendritic cells","Endothelial cells","Macrophage-like GAMs",
                                                                                      "Microglia-like GAMs","NKT-like cells",
                                                                                      "Oligodendrocytes","T cells","Tumor cells","B cells","Mural cells")),
                         PRS_type = c(rep("Newly",20)),
                         x = c("Newly","Recurrent","Newly","Recurrent","Newly","Recurrent","Newly","Recurrent","Newly","Recurrent",
                               "Newly","Recurrent","Newly","Recurrent","Newly","Recurrent","Newly","Recurrent","Newly","Recurrent"),
                         y = c(0, 5.5, 0, 50, 0, 60, 0, 50, 0, 1.75, 0, 100, 0, 7, 0, 110, 0, 1.5, 0, 85))

p1 <- ggplot(combinetable, aes(x = PRS_type, y = OptGBM, color = factor(paste0(PRS_type,"    ")))) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_color_manual(values=c("#B2182B","#2166AC")) +
  geom_blank(data = blank_data, aes(x = x, y = y, group = PRS_type)) +
  facet_wrap(vars(Celltype), nrow = 2, scales="free_y",
             labeller = labeller(Celltype = cell.labs)) +
  labs(y="OptGBM proportion (%)", x = "") +
  theme_bw() +
  theme(legend.position = "top",
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 18, vjust = 2),
        axis.text.x = element_text(color = "black", size = 18, angle = 40, vjust = 0.7, hjust = 0.7),
        axis.title.x = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black", size = 18),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.1, "cm"),
        legend.key.width = unit(1, "cm"),
        strip.text.x = element_text(size = 18),
        strip.background = element_rect(color = "black", size = 1)) 

p1

CGGA$Endo_Tumor <- CGGA$`Endothelial cells` + CGGA$`Tumor cells`
CGGA$Immune <- CGGA$`B cells`+ CGGA$`NKT-like cells` + CGGA$`T cells`
CGGA$GAM_Mural <- CGGA$`Macrophage-like GAMs` + CGGA$`Microglia-like GAMs` + CGGA$`Mural cells`
celltype2 <- c("Dendritic cells","Endo_Tumor","Immune","GAM_Mural","Oligodendrocytes")
for (i in 1:5) {
  temp_df <- CGGA[,c("CGGA_ID","PRS_type",celltype2[i])]
  colnames(temp_df) <- c("CGGA_ID","PRS_type","OptGBM")
  temp_df$Celltype <- celltype2[i]
  if(i==1){
    combinetable <- temp_df
  }else{
    combinetable <- rbind(combinetable,temp_df) %>% as.data.frame()
  }
  remove(temp_df)
}
combinetable_CGGA_C5 <- combinetable

celltype2 <- c("Dendritic cells","Endo_Tumor","Immune","GAM_Mural","Oligodendrocytes")
combinetable_CGGA_C5$Celltype <- factor(combinetable_CGGA_C5$Celltype, levels = c("Dendritic cells","Endo_Tumor","Immune","GAM_Mural","Oligodendrocytes"))
cell.labs <- c("Dendritic cells","Tumor/Endothelial\ncells","Immune cells","GAMs/Mural cells","Oligodendrocytes")
names(cell.labs) <- c("Dendritic cells","Endo_Tumor","Immune","GAM_Mural","Oligodendrocytes")
blank_data <- data.frame(Celltype = factor(c("Dendritic cells","Dendritic cells","Endo_Tumor","Endo_Tumor","Immune","Immune","GAM_Mural",
                                             "GAM_Mural","Oligodendrocytes","Oligodendrocytes"),
                                           levels = c("Dendritic cells","Endo_Tumor","Immune","GAM_Mural","Oligodendrocytes")),
                         PRS_type = c(rep("Newly",10)),
                         x = c("Newly","Recurrent","Newly","Recurrent","Newly","Recurrent","Newly","Recurrent","Newly","Recurrent"),
                         y = c(0, 5, 0, 110, 0, 7.5, 0, 110, 0, 90))
p1 <- ggplot(combinetable_CGGA_C5, aes(x = PRS_type, y = OptGBM, color = factor(paste0(PRS_type,"    ")))) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_color_manual(values=c("#B2182B","#2166AC")) +
  geom_blank(data = blank_data, aes(x = x, y = y, group = PRS_type)) +
  facet_wrap(vars(Celltype), nrow = 1, scales="free_y",
             labeller = labeller(Celltype = cell.labs)) +
  labs(y="OptGBM proportion (%)", x = "") +
  theme_bw() +
  theme(legend.position = "top",
        panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 18, vjust = 2),
        axis.text.x = element_text(color = "black", size = 18, angle = 40, vjust = 0.7, hjust = 0.7),
        axis.title.x = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black", size = 18),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.1, "cm"),
        legend.key.width = unit(1, "cm"),
        strip.text.x = element_text(size = 18),
        strip.background = element_rect(color = "black", size = 1)) 
p1

pairwise.t.test(OptGBM_table$`Tumor cells`, OptGBM_table$PRS_type, p.adjust.method = "bonferroni", paired = FALSE, alternative = c("two.sided"))
median(oridata$`Tumor cells`[which(oridata$PRS_type=="Primary")])
median(oridata$`Tumor cells`[which(oridata$PRS_type=="Recurrent")])
