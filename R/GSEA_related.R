library(magrittr)
library(stringr)
library(tidyverse)
library(readxl)

TMM <- read_excel("D:/A3_glioma/GBM_deconvolution_paper/20231221_ThreeGSEDataset_GSEA_TMM.xlsx") %>% as.data.frame()
TMM_005_025 <- TMM[which(TMM$`NOM p-val`<0.05 & TMM$`FDR q-val`<0.25),] %>% as.data.frame()
TMM_NESpos <- TMM[which(TMM$NES>0),] %>% as.data.frame()
TMM_005_025_NESpos <- TMM_005_025[which(TMM_005_025$NES>0),] %>% as.data.frame()
a <- TMM_005_025_NESpos %>% group_by(NAME) %>% count()
a <- TMM_005_025_NESpos[which(TMM_005_025_NESpos$NAME %in% a$NAME[which(a$n==2)]),] %>% as.data.frame()
a <- a[which(a$celltype %in% c("MG","MP")),] %>% as.data.frame()
aa <- a %>% group_by(NAME) %>% count()
target_pathway1 <-unique(aa$NAME[which(aa$n==2)])

TMM_005_NESpos <- TMM[which(TMM$`NOM p-val`<0.05 & TMM$NES>0),] %>% as.data.frame()
TMM_005_NESpos_minMGMP <- TMM_005_NESpos[-which(TMM_005_NESpos$celltype %in% c("MP","MG")),] %>% as.data.frame()
b <- unique(TMM_005_NESpos_minMGMP$NAME)
target_pathway2 <- setdiff(target_pathway1,b)

library(RColorBrewer)
library(wesanderson)
library(stringr)
TMM <- read_excel("D:/A3_glioma/GBM_deconvolution_paper/20231221_ThreeGSEDataset_GSEA_TMM.xlsx") %>% as.data.frame()
pathway_selected <- target_pathway2
df <- TMM[which(TMM$NAME %in% pathway_selected),] %>% as.data.frame()
df <- separate(data=df, col = tag, sep = "[=]", into = c("tag","num"), remove = FALSE)
df$Generatio <- (as.numeric(gsub( "%", "", df$num)))/ 100
df$value <- "<0.05"
df$value[which(df$`NOM p-val`>=0.05)] <- ">=0.05"
df$Pvalue005 <- df$`NOM p-val`
df$Pvalue005[which(df$`NOM p-val`>=0.05)] <- 0.05
df$Symbol_1 <- ifelse(df$Symbol== "REACTOME", "RCT", df$Symbol)
df$Pathway_1 <- tolower(df$Pathway) 
df$Pathway_1 <- str_to_title(df$Pathway_1)
df$Pathway_2 <- paste0("(",df$Symbol_1,") ",df$Pathway_1)

df_arrange <- df[order(df$NAME,df$celltype),] %>% as.data.frame()
df_arrange <- data.frame(NAME = unique(df_arrange$NAME),
                         MP_NES = df_arrange$NES[which(df_arrange$celltype=="MP")],
                         MG_NES = df_arrange$NES[which(df_arrange$celltype=="MG")],
                         MP_pval = df_arrange$`NOM p-val`[which(df_arrange$celltype=="MP")],
                         MG_pval = df_arrange$`NOM p-val`[which(df_arrange$celltype=="MG")])
df_arrange$ave_NES <- (df_arrange$MP_NES+df_arrange$MG_NES)/2
df_arrange$ave_pval <- (df_arrange$MP_pval+df_arrange$MG_pval)/2
df_arrange <- df_arrange[order(df_arrange$ave_pval, decreasing = FALSE),] %>% as.data.frame()

df$celltype <- factor(df$celltype, levels = c("BC", "TC", "NKT","EC", "MC", "OD", "Tumor","MP", "MG","DC" ))
df$Pathway_2 <- factor(df$Pathway_2, levels = rev(c("(WP) Tlr4_signaling_and_tolerance",
                                                    "(WP) Ltf_danger_signal_response_pathway",
                                                    "(WP) Ldl_influence_on_cd14_and_tlr4",
                                                    "(WP) Quercetin_and_nfkb_ap1_induced_apoptosis",
                                                    "(WP) Pparalpha_pathway",
                                                    "(WP) Agerage_pathway",
                                                    "(WP) Degradation_pathway_of_sphingolipids_including_diseases",
                                                    "(WP) Transcription_factor_regulation_in_adipogenesis",
                                                    "(PID) Pdgfrb_pathway", 
                                                    "(PID) Il1_pathway",
                                                    "(WP) Il1_signaling_pathway",
                                                    "(WP) Il3_signaling_pathway",
                                                    "(PID) Il6_7_pathway",
                                                    "(PID) Fcer1_pathway",
                                                    "(PID) Bcr_5pathway",
                                                    "(WP) Antiviral_and_antiinflammatory_effects_of_nrf2_on_sarscov2_pathway",
                                                    "(WP) Mirna_role_in_immune_response_in_sepsis",
                                                    "(WP) Hostpathogen_interaction_of_human_coronaviruses_interferon_induction",
                                                    "(WP) Mitochondrial_immune_response_to_sarscov2",
                                                    "(WP) Hostpathogen_interaction_of_human_coronaviruses_autophagy",
                                                    "(WP) Type_i_interferon_induction_and_signaling_during_sarscov2_infection",
                                                    "(WP) 4hydroxytamoxifen_dexamethasone_and_retinoic_acids_regulation_of_p27_expression",
                                                    "(WP) Glycosaminoglycan_degradation",
                                                    "(WP) Glycosylation_and_related_congenital_defects",
                                                    "(WP) Follicle_stimulating_hormone_fsh_signaling_pathway")))

p <- ggplot(df,aes(x=celltype ,y=Pathway_2,colour= NES,size= Pvalue005,shape = value))+
  geom_point(alpha = 0.8)+ 
  scale_size_continuous(range = c(1,5),breaks = seq(0,0.05,by = 0.01),
                        trans = "reverse")+
  theme_bw() + 
  xlab("")+ ylab("")+
  scale_x_discrete(labels =  c("B cells",
                               "T cells",
                               "NKT-like cells",
                               "Endothelial cells",
                               "Mural cells",
                               "Oligodendrocytes",
                               "Tumor cells",
                               "Macrophage-like GAMs",
                               "Microglia-like GAMs",
                               "Dendritic cells"))+
  # scale_y_discrete(labels = rev(c("(WP) TLR4 signaling and tolerance",
  #                                 "(WP) LTF danger signal response pathway",
  #                                 "(WP) LDL influence on CD14 and TLR4",
  #                                 paste0("(WP) Quercetin and NF",expression(kappa),"B/AP1 induced apoptosis"),
  #                                 paste0("(WP) PPAR",expression(alpha)," pathway"),
  #                                 "(WP) AGE/RAGE pathway",
  #                                 "(WP) Degradation pathway of sphingolipids including diseases",
  #                                 "(WP) Transcription factor regulation in adipogenesis",
  #                                 "(PID) PDGFRB pathway", 
  #                                 "(PID) IL1 pathway",
  #                                 "(WP) IL1 signaling pathway",
  #                                 "(WP) IL3 signaling pathway",
  #                                 "(PID) IL6/IL7 pathway",
  #                                 "(PID) FCER1 pathway",
  #                                 "(PID) BCR signaling pathway",
  #                                 "(WP) Effects of NRF2 on SARSCoV-2 pathway",
  #                                 "(WP) MiRNA role in immune response in sepsis",
  #                                 "(WP) Hostpathogen interaction-coronaviruses interferon induction",
  #                                 "(WP) Mitochondrial immune response to SARSCoV-2",
  #                                 "(WP) Hostpathogen interaction-coronaviruses autophagy",
  #                                 "(WP) IFN induction and signaling during SARSCoV-2 infection",
  #                                 "(WP) Regulation of p27 expression",
  #                                 "(WP) Glycosaminoglycan degradation",
  #                                 "(WP) Glycosylation and related congenital defects",
  #                                 "(WP) Follicle stimulating hormone FSH signaling pathway")))+
  labs(size = "P value", shape = "Value")+
  scale_color_gradient2(midpoint=0, low="#2166AC", mid="#CCCCCC", high="#B2182B", limits = c(-1,1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 11,color = "black"),
        axis.text.y = element_text(size = 11,color = "black")) 
p

library(writexl)
library(readr)
library(readxl)
write_xlsx(df, "D:/A3_glioma/GBM_deconvolution_paper/Figure7_MPMG_specific_signaling.xlsx")

#dat.gct <- read.delim(file="D:/ssGSEA/2024-01-17/20231222_TCGA_GBM_TMM_ssGSEA-scores.gct", skip=2) %>% as.data.frame()
dat.gct <- read.delim(file="D:/A3_glioma//GBM_deconvolution_paper/20240217_CGGA_GBM_TMM_ssGSEA-scores.gct", skip=2) %>% as.data.frame()
score <- dat.gct[,c(1083:1621)] %>% as.data.frame()
row.names(score) <- dat.gct[,1]
Figure7_MPMG_specific_signaling <- read_excel("D:/A3_glioma/GBM_deconvolution_paper/Figure7_MPMG_specific_signaling.xlsx") %>% as.data.frame()
#TCGA_Survival <- read_csv("D:/A3_glioma/GBM_deconvolution_paper/TCGA_Survival_CellType_145RawCount_CombineTable.csv") %>% as.data.frame()
CGGA_Survival <- read_csv("D:/A3_glioma/GBM_deconvolution_paper/CGGA_Survival_CellType_153_RawCount_CombineTable.csv") %>% as.data.frame()
colnames(score) <- str_replace_all(colnames(score), pattern = "\\.", replacement = "-")
score <- score %>% t() %>% as.data.frame()
score$PATIENT_ID <- row.names(score)
#OS_data <- TCGA_Survival[,c(13,20,21,22,23)] %>% as.data.frame()
OS_data <- CGGA_Survival[,c(2,17,18,19,20)] %>% as.data.frame()
#OS_score <- merge(x=OS_data, y=score, by.x = "PATIENT_ID", by.y = "PATIENT_ID", all = TRUE) %>% as.data.frame()
OS_score <- merge(x=OS_data, y=score, by.x = "PATIENT_ID", by.y = "PATIENT_ID", all.x = TRUE) %>% as.data.frame()
row.names(OS_score) <- OS_score$PATIENT_ID

library(survival)
library(survminer)

oridata <- OS_score
low_25 <- ceiling(length(oridata$PATIENT_ID)*0.25)
high_25 <- floor(length(oridata$PATIENT_ID)*0.75)

signal_name <- colnames(oridata)[6:30]
cname_list <- c("L_Number")
for (i in 1:25) {
  oridata_short <- oridata[,c(2,3,4,5,i+5)] %>% as.data.frame()
  colnames(oridata_short)[5] <- "target_sig"
  oridata_short <- oridata_short[order(oridata_short$target_sig),] %>% as.data.frame()
  pval_OS <- c()
  median_H_OS <- c()
  median_L_OS <- c()
  #  pval_PFS <- c()
  #  median_H_PFS <- c()
  #  median_L_PFS <- c()
  for (j in low_25:high_25) {
    oridata_short$S_group <- "H"
    oridata_short$S_group[1:j] <- "L"
    fit_OS <- surv_fit(Surv(OS_MONTHS,OS_STATUS) ~ S_group, data = oridata_short)
    pval_OS <- append(pval_OS, round(surv_pvalue(fit = fit_OS)[1,2],digits = 5))
    median_H_OS <- append(median_H_OS, surv_median(fit_OS)[1,2])
    median_L_OS <- append(median_L_OS, surv_median(fit_OS)[2,2])
    #    fit_PFS <- surv_fit(Surv(PFS_MONTHS,PFS_STATUS) ~ S_group, data = oridata_short)
    #    pval_PFS <- append(pval_PFS, round(surv_pvalue(fit = fit_PFS)[1,2],digits = 5))
    #    median_H_PFS <- append(median_H_PFS, surv_median(fit_PFS)[1,2])
    #    median_L_PFS <- append(median_L_PFS, surv_median(fit_PFS)[2,2])
  }
  if(i==1){
    sumtable_final <- data.frame(Number=c(low_25:high_25),
                                 Pvalue_OS=pval_OS,
                                 Median_surv_H_OS=median_H_OS,
                                 Median_surv_L_OS=median_L_OS)
    #                                 Pvalue_PFS=pval_PFS,
    #                                 Median_surv_H_PFS=median_H_PFS,
    #                                 Median_surv_L_PFS=median_L_PFS)
  }else{
    sumstable <- data.frame(Pvalue_OS=pval_OS,
                            Median_surv_H_OS=median_H_OS,
                            Median_surv_L_OS=median_L_OS)
    #                            Pvalue_PFS=pval_PFS,
    #                            Median_surv_H_PFS=median_H_PFS,
    #                            Median_surv_L_PFS=median_L_PFS)
    sumtable_final <- cbind(sumtable_final,sumstable) %>% as.data.frame()
  }
  cname_list_temp <- c(paste0(signal_name[i],"_Pvalue_OS"),paste0(signal_name[i],"_Median_surv_H_OS"),paste0(signal_name[i],"_Median_surv_L_OS"))
  #                       paste0(signal_name[i],"_Pvalue_PFS"),paste0(signal_name[i],"_Median_surv_H_PFS"),paste0(signal_name[i],"_Median_surv_L_PFS"))
  cname_list <- append(cname_list,cname_list_temp)
  #  remove(sumstable,cname_list_temp,oridata_short,fit_OS,fit_PFS)
  remove(sumstable,cname_list_temp,oridata_short,fit_OS)
}

colnames(sumtable_final) <- cname_list
library(writexl)
#write_xlsx(sumtable_final, "D:/A3_glioma/GBM_deconvolution_paper/TCGA_ssGSEA_25signal_OSPFS.xlsx")
write_xlsx(sumtable_final, "D:/A3_glioma/GBM_deconvolution_paper/CGGA_ssGSEA_25signal_OS.xlsx")


oridata <- OS_score
oridata <- oridata[order(oridata$WP_LTF_DANGER_SIGNAL_RESPONSE_PATHWAY),] %>% as.data.frame()
oridata$OS_group <- "H"
oridata$OS_group[1:73] <- "L"
fit_OS <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ OS_group, data = oridata)
pvalue_OS <- surv_pvalue(fit = fit_OS, data = oridata)
OS_fig <- ggsurvplot(fit_OS,
                     pval = paste0("p=",format(pvalue_OS$pval, digits = 4)),
                     pval.coord = c(50,0.8),
                     risk.table = TRUE,
                     #risk.table = FALSE,
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
OS_fig

oridata$PFS_group <- "H"
oridata$PFS_group[1:73] <- "L"
fit_PFS <- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ PFS_group, data = oridata)
pvalue_PFS <- surv_pvalue(fit = fit_PFS, data = oridata)
PFS_fig <- ggsurvplot(fit_PFS,
                      pval = paste0("p=",format(pvalue_PFS$pval, digits = 4)),
                      pval.coord = c(20,0.8),
                      risk.table = TRUE,
                      #risk.table = FALSE,
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
PFS_fig

library(readr)
sig25 <- read_delim("D:/ssGSEA/ssGSEA2.0-master/db/msigdb/c2.20240117.symbols.gmt", 
                    delim = "\t", escape_double = FALSE, 
                    col_names = FALSE, trim_ws = TRUE) %>% as.data.frame()
#OS_sig_graph
os_sig_name <- c("(WP) Glycosaminoglycan degradation",
                 "(WP) LTF danger signal response pathway",
                 "(WP) LDL influence on CD14 and TLR4",
                 "(WP) Transcription factor regulation in adipogenesis",
                 "(WP) AGE/RAGE pathway",
                 "(WP) Hostpathogen interaction-coronaviruses autophagy",
                 "(WP) Mitochondrial immune response to SARSCoV-2",
                 "(WP) Degradation pathway of sphingolipids including diseases")
Low_numner <- data.frame(OS = c(74,108,108,94,51,53,105,99),
                         PFS = c(75,106,108,101,104,67,54,107),
                         row.names = factor(os_sig_name, levels = c("(WP) Glycosaminoglycan degradation",
                                                                    "(WP) LTF danger signal response pathway",
                                                                    "(WP) LDL influence on CD14 and TLR4",
                                                                    "(WP) Transcription factor regulation in adipogenesis",
                                                                    "(WP) AGE/RAGE pathway",
                                                                    "(WP) Hostpathogen interaction-coronaviruses autophagy",
                                                                    "(WP) Mitochondrial immune response to SARSCoV-2",
                                                                    "(WP) Degradation pathway of sphingolipids including diseases")))
os_sig_name_cap <- c("WP_GLYCOSAMINOGLYCAN_DEGRADATION",
                     "WP_LTF_DANGER_SIGNAL_RESPONSE_PATHWAY",
                     "WP_LDL_INFLUENCE_ON_CD14_AND_TLR4",
                     "WP_TRANSCRIPTION_FACTOR_REGULATION_IN_ADIPOGENESIS",
                     "WP_AGE_RAGE_PATHWAY",
                     "WP_HOST_PATHOGEN_INTERACTION_OF_HUMAN_CORONAVIRUSES_AUTOPHAGY",
                     "WP_MITOCHONDRIAL_IMMUNE_RESPONSE_TO_SARS_COV_2",
                     "WP_DEGRADATION_PATHWAY_OF_SPHINGOLIPIDS_INCLUDING_DISEASES")
k <- 1
OS_fig <- list()
for (j in 1:length(os_sig_name)) {
  oridata <- OS_score[,append(c(2,3,4,5),which(colnames(OS_score)==os_sig_name_cap[j]))] %>% as.data.frame()
  colnames(oridata)[5] <- "target_sig"
  oridata <- oridata[order(oridata$target_sig),] %>% as.data.frame()
  i <- Low_numner$OS[j]
  oridata$OS_group <- "H"
  oridata$OS_group[1:i] <- "L"
  fit_OS <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ OS_group, data = oridata)
  pvalue_OS <- surv_pvalue(fit = fit_OS, data = oridata)
  OS_fig[[k]] <- ggsurvplot(fit_OS,
                            pval = paste0("p=",format(pvalue_OS$pval, digits = 4)),
                            pval.coord = c(50,0.2),
                            #risk.table = TRUE,
                            risk.table = FALSE,
                            legend.title = os_sig_name[j],
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
  remove(oridata,i,fit_OS,pvalue_OS)
  k <- k+1
}
gg <- arrange_ggsurvplots(OS_fig, print = FALSE,
                          ncol = 2, nrow = 4, risk.table.height = 0.3, surv.plot.height = 0.7)
ggsave("D:/A3_glioma/GBM_deconvolution_paper/ssGSEA_25signal_OS.pdf", gg,
       width = 8, height = 15, units = "in")

k <- 1
PFS_fig <- list()
for (j in 1:length(os_sig_name)) {
  oridata <- OS_score[,append(c(2,3,4,5),which(colnames(OS_score)==os_sig_name_cap[j]))] %>% as.data.frame()
  colnames(oridata)[5] <- "target_sig"
  oridata <- oridata[order(oridata$target_sig),] %>% as.data.frame()
  i <- Low_numner$PFS[j]
  oridata$PFS_group <- "H"
  oridata$PFS_group[1:i] <- "L"
  fit_PFS <- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ PFS_group, data = oridata)
  pvalue_PFS <- surv_pvalue(fit = fit_PFS, data = oridata)
  PFS_fig[[k]] <- ggsurvplot(fit_PFS,
                             pval = paste0("p=",format(pvalue_PFS$pval, digits = 4)),
                             pval.coord = c(28,0.3),
                             #risk.table = TRUE,
                             risk.table = FALSE,
                             legend.title = os_sig_name[j],
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
  remove(oridata,i,fit_PFS,pvalue_PFS)
  k <- k+1
}
gg <- arrange_ggsurvplots(PFS_fig, print = FALSE,
                          ncol = 2, nrow = 4, risk.table.height = 0.3, surv.plot.height = 0.7)
ggsave("D:/A3_glioma/GBM_deconvolution_paper/ssGSEA_25signal_PFS.pdf", gg,
       width = 8, height = 15, units = "in")

WP_GLYCOSAMINOGLYCAN_DEGRADATION <- c("ARSB","GALNS","GLB1","GNS","GUSB","HEXA","HEXB","HGSNAT","HPSE","HPSE2","HYAL1","HYAL2","HYAL4","IDS","IDUA","NAGLU","SGSH")
WP_LTF_DANGER_SIGNAL_RESPONSE_PATHWAY <- c("AGER","CD14","CXCL8","IFNA1","IFNB1","IL1A","IL1B","IL6","IRAK1","IRAK4","LTF","MAPK1","MYD88","NFKB1","TLR2","TLR4",
                                           "TNF","TRAF6","TREM1")
WP_LDL_INFLUENCE_ON_CD14_AND_TLR4 <- c("AKT1","ATF4","CCL2","CD14","CREB1","CREB3","CREB3L1","CREB3L2","CREB3L3","CREB3L4","CREB5","IL10","IL1B","IL6","MAPK11",
                                       "MAPK12","MAPK13","MAPK14","NFKB1","NFKB2","REL","RELA","RELB","TLR4")
WP_TRANSCRIPTION_FACTOR_REGULATION_IN_ADIPOGENESIS <- c("ADIPOQ","CEBPA","CEBPB","CEBPD","CREB1","FOXO1","IL6","INSR","IRS1","IRS2","LEP","LPIN1","MAPK8","NR3C1",
                                                        "NRIP1","PCK2","PPARG","PPARGC1A","RXRA","SLC2A4","TNF","TWIST1")
WP_AGE_RAGE_PATHWAY <- c("AGER","AKT1","ALPL","ATF2","CASP3","CASP8","CASP9","CDC42","CHUK","CYCS","DDOST","DIAPH1","EGFR","EZR","FOXO1","FOXO4","HIF1A","IKBKB",
                         "INHBB","INS","INSR","IRAK4","IRS1","JAK2","JUN","LGALS3","MAP2K1","MAPK1","MAPK14","MAPK3","MAPK8","MAPK9","MMP13","MMP14","MMP2","MMP7",
                         "MMP9","MSN","MSR1","MYD88","NCF1","NFKB1","NFKBIA","NOS2","NOS3","PLA2G4A","PRKCA","PRKCB","PRKCD","PRKCZ","RAC1","RAF1","RELA","RHOA",
                         "ROCK1","SHC1","SMAD2","SMAD3","SOD1","SP1","SRC","STAT1","STAT3","STAT5A","STAT5B","TIRAP")
WP_HOST_PATHOGEN_INTERACTION_OF_HUMAN_CORONAVIRUSES_AUTOPHAGY <- c("ATG10","ATG12","ATG13","ATG16L1","ATG16L2","ATG3","ATG4A","ATG5","ATG7","BECN1","MAP1LC3A",
                                                                   "MTOR","PIK3C3","PIK3R4","RB1CC1","ULK1","ULK2","WIPI1","ZFYVE1")
WP_MITOCHONDRIAL_IMMUNE_RESPONSE_TO_SARS_COV_2 <- c("ACAD9","ACE","ACE2","AGTR1","AGTR2","BCS1L","CGAS","CTSL","ECSIT","IFIH1","IKBKE","IRF3","IRF7","MAVS",
                                                    "NDUFAF1","NDUFB9","NFKB1","NFKB2","NLRX1","NOX1","PHB1","PHB2","REN","RIGI","STING1","TBK1","TICAM1","TLR3",
                                                    "TLR7","TMPRSS2","TOMM70","TRAF3","TRAF6")
WP_DEGRADATION_PATHWAY_OF_SPHINGOLIPIDS_INCLUDING_DISEASES <- c("ARSA","GALC","GBA1","GLA","GLB1","GM2A","HEXA","HEXB","LIPA","NEU1","NEU2","NEU3","NEU4","NPC1",
                                                                "NPC2","PSAP","SCARB2")
# sig_list <- list(WP_GLYCOSAMINOGLYCAN_DEGRADATION,WP_LTF_DANGER_SIGNAL_RESPONSE_PATHWAY,WP_LDL_INFLUENCE_ON_CD14_AND_TLR4,
#                  WP_TRANSCRIPTION_FACTOR_REGULATION_IN_ADIPOGENESIS,WP_AGE_RAGE_PATHWAY,WP_HOST_PATHOGEN_INTERACTION_OF_HUMAN_CORONAVIRUSES_AUTOPHAGY,
#                  WP_MITOCHONDRIAL_IMMUNE_RESPONSE_TO_SARS_COV_2,WP_DEGRADATION_PATHWAY_OF_SPHINGOLIPIDS_INCLUDING_DISEASES)
sig_list <- list(WP_LTF_DANGER_SIGNAL_RESPONSE_PATHWAY,WP_LDL_INFLUENCE_ON_CD14_AND_TLR4,
                 WP_TRANSCRIPTION_FACTOR_REGULATION_IN_ADIPOGENESIS,WP_AGE_RAGE_PATHWAY)
sig_list_total <- c()
for (i in 1:4) {
  sig_list_total <- append(sig_list_total,sig_list[[i]])
}
sig_list_total <- unique(sig_list_total)
# temp_df <- data.frame(row.names = sig_list_total,
#                       WP_GLYCOSAMINOGLYCAN_DEGRADATION = ifelse(sig_list_total %in% WP_GLYCOSAMINOGLYCAN_DEGRADATION, 1, 0),
#                       WP_LTF_DANGER_SIGNAL_RESPONSE_PATHWAY = ifelse(sig_list_total %in% WP_LTF_DANGER_SIGNAL_RESPONSE_PATHWAY, 1, 0),
#                       WP_LDL_INFLUENCE_ON_CD14_AND_TLR4 = ifelse(sig_list_total %in% WP_LDL_INFLUENCE_ON_CD14_AND_TLR4, 1, 0),
#                       WP_TRANSCRIPTION_FACTOR_REGULATION_IN_ADIPOGENESIS = ifelse(sig_list_total %in% WP_TRANSCRIPTION_FACTOR_REGULATION_IN_ADIPOGENESIS, 1, 0),
#                       WP_AGE_RAGE_PATHWAY = ifelse(sig_list_total %in% WP_AGE_RAGE_PATHWAY, 1, 0),
#                       WP_HOST_PATHOGEN_INTERACTION_OF_HUMAN_CORONAVIRUSES_AUTOPHAGY = ifelse(sig_list_total %in% WP_HOST_PATHOGEN_INTERACTION_OF_HUMAN_CORONAVIRUSES_AUTOPHAGY, 1, 0),
#                       WP_MITOCHONDRIAL_IMMUNE_RESPONSE_TO_SARS_COV_2 = ifelse(sig_list_total %in% WP_MITOCHONDRIAL_IMMUNE_RESPONSE_TO_SARS_COV_2, 1, 0),
#                       WP_DEGRADATION_PATHWAY_OF_SPHINGOLIPIDS_INCLUDING_DISEASES = ifelse(sig_list_total %in% WP_DEGRADATION_PATHWAY_OF_SPHINGOLIPIDS_INCLUDING_DISEASES, 1, 0))
temp_df <- data.frame(row.names = sig_list_total,
                      WP_LTF_DANGER_SIGNAL_RESPONSE_PATHWAY = ifelse(sig_list_total %in% WP_LTF_DANGER_SIGNAL_RESPONSE_PATHWAY, 1, 0),
                      WP_LDL_INFLUENCE_ON_CD14_AND_TLR4 = ifelse(sig_list_total %in% WP_LDL_INFLUENCE_ON_CD14_AND_TLR4, 1, 0),
                      WP_TRANSCRIPTION_FACTOR_REGULATION_IN_ADIPOGENESIS = ifelse(sig_list_total %in% WP_TRANSCRIPTION_FACTOR_REGULATION_IN_ADIPOGENESIS, 1, 0),
                      WP_AGE_RAGE_PATHWAY = ifelse(sig_list_total %in% WP_AGE_RAGE_PATHWAY, 1, 0)
)

temp_df <- temp_df %>% t() %>% as.data.frame()
library(pheatmap)
php <- pheatmap(temp_df,
                color = c("lightgrey","black"),
                cellwidth = 10,
                cellheight = 10,
                border_color = "black",
                clustering_method = "ward.D2")
ggsave("D:/A3_glioma/GBM_deconvolution_paper/genename_8signal_fig.pdf", php,
       width = 4, height = 35, units = "in")

#IHC score and survival--------------------------------------------------------------------------------------------------------------------
library(readxl)
library(stringr)
library(magrittr)
library(tidyverse)

TCGG_clinical <- read_excel("D:/A3_glioma/A3B_clinical_20240215.xlsx") %>% as.data.frame()
IHC_total <- read_excel("D:/A3_glioma/GBM_deconvolution_paper/IHC/IHC_total_20240215.xlsx") %>% as.data.frame()
a <- lapply(str_split(IHC_total$labnumber...1, pattern = "-"), function(x)x[1]) %>% unlist() %>% as.vector()
IHC_total$SAMPLE_ID <- a
IHC_total$SAMPLE_ID <- str_replace_all(IHC_total$SAMPLE_ID, pattern = "w", replacement = "W")
TCGG_clinical$W_sample <- str_replace_all(TCGG_clinical$W_sample, pattern = "w", replacement = "W")
IHC_total$pathology[which(is.na(IHC_total$pathology)==TRUE)] <- "Recurrent"
tempdf <- TCGG_clinical[,c("W_sample","PFS_status","PFS_days","OS_status","OS_days")] %>% as.data.frame()
IHC_total <- merge(IHC_total, tempdf, by.x = "SAMPLE_ID", by.y = "W_sample", all.x = TRUE) %>% as.data.frame()
IHC_total$labnumber...1 <- NULL
IHC_total$labnumber...12 <- NULL
IHC_total$BVD_1mm2[which(IHC_total$BVD_1mm2=="X")] <- NA
IHC_total$EC_area[which(IHC_total$EC_area=="X")] <- NA
IHC_total$BVD_1mm2 <- as.numeric(IHC_total$BVD_1mm2)
IHC_total$EC_area <- as.numeric(IHC_total$EC_area)
IHC_total$PFS_MONTHS <- IHC_total$PFS_days/30.43
IHC_total$OS_MONTHS <- IHC_total$OS_days/30.43
IHC_total$PFS_STATUS <- IHC_total$PFS_status
IHC_total$OS_STATUS <- IHC_total$OS_status
library(writexl)
write_xlsx(IHC_total, "D:/A3_glioma/GBM_deconvolution_paper/TCGG_IHC_score_clinical.xlsx")

library(survival)
library(survminer)

oridata <- IHC_total[which(is.na(IHC_total$OS_status)==FALSE),] %>% as.data.frame()
oridata <- oridata[which(is.na(oridata$CD68_score3)==FALSE),] %>% as.data.frame()
oridata <- oridata[which(is.na(oridata$EC_area)==FALSE),] %>% as.data.frame()
low_25 <- ceiling(length(oridata$SAMPLE_ID)*0.25)
high_25 <- floor(length(oridata$SAMPLE_ID)*0.75)

oridata <- oridata[order(oridata$BVD_1mm2),] %>% as.data.frame()
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
remove(pval_OS,median_H_OS,median_L_OS,pval_PFS,median_H_PFS,median_L_PFS,sumtable,fit_OS,fit_PFS,i)


oridata_ori <- IHC_total[which(is.na(IHC_total$OS_status)==FALSE),] %>% as.data.frame()
oridata <- oridata_ori[which(is.na(oridata_ori$CD68_score3)==FALSE),] %>% as.data.frame()
a <- sort((oridata$CD68_score3 %>% unique()), decreasing = FALSE)
pval_OS <- c()
median_H_OS <- c()
median_L_OS <- c()
pval_PFS <- c()
median_H_PFS <- c()
median_L_PFS <- c()
for (i in 1:length(a)) {
  oridata$S_group <- "H"
  oridata$S_group[which(oridata$CD68_score3 %in% a[1:i])] <- "L"
  fit_OS <- surv_fit(Surv(OS_MONTHS,OS_STATUS) ~ S_group, data = oridata)
  pval_OS <- append(pval_OS, round(surv_pvalue(fit = fit_OS)[1,2],digits = 5))
  median_H_OS <- append(median_H_OS, surv_median(fit_OS)[1,2])
  median_L_OS <- append(median_L_OS, surv_median(fit_OS)[2,2])
  fit_PFS <- surv_fit(Surv(PFS_MONTHS,PFS_STATUS) ~ S_group, data = oridata)
  pval_PFS <- append(pval_PFS, round(surv_pvalue(fit = fit_PFS)[1,2],digits = 5))
  median_H_PFS <- append(median_H_PFS, surv_median(fit_PFS)[1,2])
  median_L_PFS <- append(median_L_PFS, surv_median(fit_PFS)[2,2])
}
sumtable <- data.frame(Number=c(1:length(a)),
                       Pvalue_OS=pval_OS,
                       Median_surv_H_OS=median_H_OS,
                       Median_surv_L_OS=median_L_OS,
                       Pvalue_PFS=pval_PFS,
                       Median_surv_H_PFS=median_H_PFS,
                       Median_surv_L_PFS=median_L_PFS)

oridata_ori <- IHC_total[which(is.na(IHC_total$OS_status)==FALSE),] %>% as.data.frame()
oridata_BVA <- oridata_ori[which(is.na(oridata_ori$EC_area)==FALSE),] %>% as.data.frame()
oridata_CD68 <- oridata_ori[which(is.na(oridata_ori$CD68_score3)==FALSE),] %>% as.data.frame()
a <- sort((oridata_CD68$CD68_score3 %>% unique()), decreasing = FALSE)
#oridata$OS_group <- "H"
#oridata$OS_group[1:39] <- "L"
#oridata$OS_group[which(oridata$CD68_score3 %in% a[1:6])] <- "L"
oridata_CD68$GAM_group <- "H"
oridata_CD68$GAM_group[which(oridata_CD68$CD68_score3 %in% a[1:6])] <- "L"
oridata_BVA <- oridata_BVA[order(oridata_BVA$EC_area),] %>% as.data.frame()
oridata_BVA$OS_BVA_group <- "H"
oridata_BVA$OS_BVA_group[1:39] <- "L"
oridata_BVA$PFS_BVA_group <- "H"
oridata_BVA$PFS_BVA_group[1:34] <- "L"
#oridata_BVD <- oridata_BVD[order(oridata_BVD$EC_area),] %>% as.data.frame()
#oridata_BVD$BVD_group <- "H"
#oridata_BVD$BVD_group[1:33] <- "L"
#oridata <- merge(oridata_BVD[,c("SAMPLE_ID","BVD_group")],oridata_CD68,by.x = "SAMPLE_ID", by.y = "SAMPLE_ID", all = FALSE) %>% as.data.frame()
oridata <- merge(oridata_BVA[,c("SAMPLE_ID","OS_BVA_group","PFS_BVA_group")],oridata_CD68,by.x = "SAMPLE_ID", by.y = "SAMPLE_ID", all = FALSE) %>% as.data.frame()

fit_OS <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ GAM_group+OS_BVA_group, data = oridata)
#pvalue_OS <- surv_pvalue(fit = fit_OS, data = oridata)
OS_fig <- ggsurvplot(fit_OS,
                     #pval = paste0("p=",format(pvalue_OS$pval, digits = 4)),
                     #pval.coord = c(50,0.8),
                     #risk.table = TRUE,
                     risk.table = FALSE,
                     #surv.median.line = c("hv"),
                     legend = c("top"), 
                     ylab="OS Probability",
                     xlab="Months after first diagnosis",
                     #palette= c("#B2182B","#2166AC"),
                     palette= c("#B2182B","#F4A582","#2166AC","#A6CEE3"),
                     #legend.labs = c("High","Low"),
                     font.title = 16,
                     font.x =  16,
                     font.y = 16,
                     font.tickslab = 16,
                     tables.y.text = FALSE)
OS_fig

oridata$Two_group_OS <- paste0(oridata$GAM_group,oridata$OS_BVA_group)
oridata$Two_group_PFS <- paste0(oridata$GAM_group,oridata$PFS_BVA_group)

group_num <- c("LH","LL") 
fit_OS <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ Two_group_OS, data = oridata[which(oridata$Two_group_OS %in% group_num),])
#fit_OS
fit_PFS <- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ Two_group_PFS, data = oridata[which(oridata$Two_group_PFS %in% group_num),])
#fit_PFS
round(surv_pvalue(fit = fit_OS)[1,2],digits = 5)
round(surv_pvalue(fit = fit_PFS)[1,2],digits = 5)

#oridata$PFS_group <- "H"
#oridata$PFS_group[1:34] <- "L"
#oridata$PFS_group[which(oridata$CD68_score3 %in% a[1:6])] <- "L"
oridata_CD68$GAM_group <- "H"
oridata_CD68$GAM_group[which(oridata_CD68$CD68_score3 %in% a[1:6])] <- "L"
oridata_BVA <- oridata_BVA[order(oridata_BVA$EC_area),] %>% as.data.frame()
oridata_BVA$PFS_BVA_group <- "H"
oridata_BVA$PFS_BVA_group[1:34] <- "L"
oridata <- merge(oridata_BVA[,c("SAMPLE_ID","PFS_BVA_group")],oridata_CD68,by.x = "SAMPLE_ID", by.y = "SAMPLE_ID", all = FALSE) %>% as.data.frame()
fit_PFS <- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ GAM_group+PFS_BVA_group, data = oridata)
#pvalue_PFS <- surv_pvalue(fit = fit_PFS, data = oridata)
PFS_fig <- ggsurvplot(fit_PFS,
                      #pval = paste0("p=",format(pvalue_PFS$pval, digits = 4)),
                      #pval.coord = c(20,0.8),
                      risk.table = TRUE,
                      #risk.table = FALSE,
                      surv.median.line = c("hv"),
                      legend = c("top"), 
                      ylab="PFS Probability",
                      xlab="Months after first diagnosis",
                      #palette= c("#B2182B","#2166AC"),
                      palette= c("#B2182B","#F4A582","#2166AC","#A6CEE3"),
                      #legend.labs = c("High","Low"),
                      xlim = c(0,60),
                      font.title = 16,
                      font.x =  16,
                      font.y = 16,
                      font.tickslab = 16,
                      tables.y.text = FALSE)
PFS_fig

arrange_ggsurvplots(S_fig, print = TRUE,
                    ncol = 2, nrow = 2, risk.table.height = 0.3, surv.plot.height = 0.7)

#blood vessel signaling--------------------------------------------------------------------------------------
library(magrittr)
library(stringr)
library(tidyverse)
library(readxl)

TMM <- read_excel("D:/A3_glioma/GBM_deconvolution_paper/20231221_ThreeGSEDataset_GSEA_TMM.xlsx") %>% as.data.frame()
TMM_005_025 <- TMM[which(TMM$`NOM p-val`<0.05 & TMM$`FDR q-val`<0.25),] %>% as.data.frame()
#TMM_NESpos <- TMM[which(TMM$NES>0),] %>% as.data.frame()
#TMM_005_025_NESpos <- TMM_005_025[which(TMM_005_025$NES>0),] %>% as.data.frame()
#a <- TMM_005_025_NESpos %>% group_by(NAME) %>% count()
a <- TMM_005_025 %>% group_by(NAME) %>% count()
#a <- TMM_005_025_NESpos[which(TMM_005_025_NESpos$NAME %in% a$NAME[which(a$n==1)]),] %>% as.data.frame()
a <- TMM_005_025[which(TMM_005_025$NAME %in% a$NAME[which(a$n==1)]),] %>% as.data.frame()
a <- a[which(a$celltype=="EC"),] %>% as.data.frame()
target_pathway2 <-unique(a$NAME)

#GSEA TCGA GAM----------------------------------------------------------------------------------------------------------------------------------------------
library(readr)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(stringr)
library(magrittr)

TCGA <- read_csv("D:/A3_glioma/GBM_deconvolution_paper/TCGA_Survival_CellType_145RawCount_CombineTable.csv") %>% as.data.frame()
row.names(TCGA) <- TCGA$...1
colnames(TCGA)[1] <- "PATIENT_ID"
TCGA <- TCGA[,-c(13)] %>% as.data.frame()
TCGA$GAMs <- TCGA$`Macrophage-like GAMs`+TCGA$`Microglia-like GAMs`
TCGA <- TCGA[order(TCGA$GAMs),] %>% as.data.frame()
TCGA$GAM_group <- "High"
TCGA$GAM_group[c(1:44)] <- "Low"

TCGA_TMM <- read_delim("D:/ssGSEA/TCGA_TMM_20240117.txt", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE) %>% as.data.frame()
TCGA$GAM_group[c(1:44)] <- "Low"
row.names(TCGA_TMM) <- TCGA_TMM$NAME
TCGA_TMM$NAME <- NULL

genelsit <- row.names(TCGA_TMM)
High <- rowMeans(TCGA_TMM[,which(colnames(TCGA_TMM) %in% TCGA$PATIENT_ID[which(TCGA$GAM_group=="High")])])
Low <- rowMeans(TCGA_TMM[,which(colnames(TCGA_TMM) %in% TCGA$PATIENT_ID[which(TCGA$GAM_group=="Low")])])
Ratio_HighLow <- High/Low
names(Ratio_HighLow) <- genelsit
Ratio_HighLow <- sort(Ratio_HighLow, decreasing = TRUE)

a <- bitr(names(Ratio_HighLow), fromType="SYMBOL", toType="ENTREZID", OrgDb = "org.Hs.eg.db", drop = FALSE)
if(length(a$SYMBOL[which(duplicated(a$SYMBOL))])!=0){
  b <- a[-which(duplicated(a$SYMBOL)),]
}else{
  b <- a
}
names(Ratio_HighLow) <- b$ENTREZID

# m_t2g <- read.gmt("D:/ssGSEA/db/msigdb/c2.cp.wikipathways.v2023.2.Hs.entrez.gmt") %>% as.data.frame()
# targeted_term <- c("WP_LDL_INFLUENCE_ON_CD14_AND_TLR4",
#                    "WP_LTF_DANGER_SIGNAL_RESPONSE_PATHWAY",
#                    "WP_TRANSCRIPTION_FACTOR_REGULATION_IN_ADIPOGENESIS",
#                    "WP_AGE_RAGE_PATHWAY")
# targeted_m_t2g <- m_t2g[which(m_t2g$term %in% targeted_term),]

y <- gseWP(Ratio_HighLow,
           organism = "Homo sapiens",
           pAdjustMethod = "BH",
           pvalueCutoff = 1.0)

y_trans <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID') 
temp_data <- y_trans@result %>% as.data.frame()
temp_data$GeneRatio <- as.numeric(str_remove(str_remove(str_split_fixed(temp_data$leading_edge, pattern = ",",n=3)[,1], pattern = "tags="),"%"))/100
temp_data <- temp_data[which(temp_data$ID %in% c("WP2324","WP3599","WP4478","WP5272","WP4815","WP4863","WP5038")),] %>% as.data.frame()
y_trans@result <- temp_data
geneSet <- y_trans@geneSets[names(y_trans@geneSets) %in% c("WP2324","WP3599","WP4478","WP5272","WP4815","WP4863","WP5038")]
y_trans@geneSets <- geneSet

targeted_gene_dataset <- temp_data
a <- c()
for (i in c(1:length(targeted_gene_dataset$core_enrichment))) {
  a <- paste(a,targeted_gene_dataset$core_enrichment[i], sep = "/")
}
aa <- str_split(a,pattern = "/") %>% unlist()
aa <- aa[-1]
aaa <- unique(aa)
aaaa <- bitr(aaa, fromType="SYMBOL", toType="ENTREZID", OrgDb = "org.Hs.eg.db", drop = FALSE) #wiki_microglia pathways中有差別的基因
FC <- Ratio_HighLow[names(Ratio_HighLow) %in% aaaa$ENTREZID & is.na(names(Ratio_HighLow))==FALSE]
FC[which(FC>6)] <- 6

options(ggrepel.max.overlaps = Inf)
library(ggnewscale)
cne <- cnetplot(y_trans,
                showCategory = 7,
                node_label="gene",
                #cex_label_category = "none",
                cex_label_gene = 0.8,
                color_category = "#FFCC00",
                foldChange=FC,
                geneSize = 20) +
  scale_size_continuous(name = "Size",
                        range = c(3,10),
                        limits = c(0,40),
                        breaks=c(10,20,30,40)) +
  scale_color_gradient2(name='Fold Changes', low='darkgreen', high='firebrick',
                        limits= c(-6,6),
                        midpoint = 0,
                        breaks=c(-6, -3, 0, 3, 6))
cne

#Figure 7
TMM <- read_excel("F:/Glioma/single cell (download)/OurPaperData/ssGSEA/20231221_ThreeGSEDataset_GSEA_TMM.xlsx") %>% as.data.frame()
TMM$pathway <- paste0("(", TMM$Symbol, ") ", TMM$Pathway)

target_pathway2 <- c("(WP) Tlr4_signaling_and_tolerance",
                     "(WP) Ltf_danger_signal_response_pathway",
                     "(WP) Ldl_influence_on_cd14_and_tlr4",
                     "(WP) Quercetin_and_nfkb_ap1_induced_apoptosis",
                     "(WP) Pparalpha_pathway",
                     "(WP) Agerage_pathway",
                     "(WP) Degradation_pathway_of_sphingolipids_including_diseases",
                     "(WP) Transcription_factor_regulation_in_adipogenesis",
                     "(PID) Pdgfrb_pathway", 
                     "(PID) Il1_pathway",
                     "(WP) Il1_signaling_pathway",
                     "(WP) Il3_signaling_pathway",
                     "(PID) Il6_7_pathway",
                     "(PID) Fcer1_pathway",
                     "(PID) Bcr_5pathway",
                     "(WP) Antiviral_and_antiinflammatory_effects_of_nrf2_on_sarscov2_pathway",
                     "(WP) Mirna_role_in_immune_response_in_sepsis",
                     "(WP) Hostpathogen_interaction_of_human_coronaviruses_interferon_induction",
                     "(WP) Mitochondrial_immune_response_to_sarscov2",
                     "(WP) Hostpathogen_interaction_of_human_coronaviruses_autophagy",
                     "(WP) Type_i_interferon_induction_and_signaling_during_sarscov2_infection",
                     "(WP) 4hydroxytamoxifen_dexamethasone_and_retinoic_acids_regulation_of_p27_expression",
                     "(WP) Glycosaminoglycan_degradation",
                     "(WP) Glycosylation_and_related_congenital_defects",
                     "(WP) Follicle_stimulating_hormone_fsh_signaling_pathway")
pathway_selected <- toupper(target_pathway2)
df <- TMM[which(TMM$pathway %in% pathway_selected),] %>% as.data.frame()
df <- separate(data=df, col = tag, sep = "[=]", into = c("tag","num"), remove = FALSE)
df$Generatio <- (as.numeric(gsub( "%", "", df$num)))/ 100
df$value <- "<0.05"
df$value[which(df$`NOM p-val`>=0.05)] <- ">=0.05"
df$Pvalue005 <- df$`NOM p-val`
df$Pvalue005[which(df$`NOM p-val`>=0.05)] <- 0.05
df$Symbol_1 <- ifelse(df$Symbol== "REACTOME", "RCT", df$Symbol)
df$Pathway_1 <- tolower(df$Pathway) 
df$Pathway_1 <- str_to_title(df$Pathway_1)
df$Pathway_2 <- paste0("(",df$Symbol_1,") ",df$Pathway_1)

# df_arrange <- df[order(df$NAME,df$celltype),] %>% as.data.frame()
# df_arrange <- data.frame(NAME = unique(df_arrange$NAME),
#                          MP_NES = df_arrange$NES[which(df_arrange$celltype=="MP")],
#                          MG_NES = df_arrange$NES[which(df_arrange$celltype=="MG")],
#                          MP_pval = df_arrange$`NOM p-val`[which(df_arrange$celltype=="MP")],
#                          MG_pval = df_arrange$`NOM p-val`[which(df_arrange$celltype=="MG")])
# df_arrange$ave_NES <- (df_arrange$MP_NES+df_arrange$MG_NES)/2
# df_arrange$ave_pval <- (df_arrange$MP_pval+df_arrange$MG_pval)/2
# df_arrange <- df_arrange[order(df_arrange$ave_pval, decreasing = FALSE),] %>% as.data.frame()

df$celltype <- factor(df$celltype, levels = c("BC", "TC", "NKT","EC", "MC", "OD", "Tumor","MP", "MG","DC" ))
df$Pathway_2 <- factor(df$Pathway_2, levels = rev(c("(WP) Tlr4_signaling_and_tolerance",
                                                    "(WP) Ltf_danger_signal_response_pathway",
                                                    "(WP) Ldl_influence_on_cd14_and_tlr4",
                                                    "(WP) Quercetin_and_nfkb_ap1_induced_apoptosis",
                                                    "(WP) Pparalpha_pathway",
                                                    "(WP) Agerage_pathway",
                                                    "(WP) Degradation_pathway_of_sphingolipids_including_diseases",
                                                    "(WP) Transcription_factor_regulation_in_adipogenesis",
                                                    "(PID) Pdgfrb_pathway", 
                                                    "(PID) Il1_pathway",
                                                    "(WP) Il1_signaling_pathway",
                                                    "(WP) Il3_signaling_pathway",
                                                    "(PID) Il6_7_pathway",
                                                    "(PID) Fcer1_pathway",
                                                    "(PID) Bcr_5pathway",
                                                    "(WP) Antiviral_and_antiinflammatory_effects_of_nrf2_on_sarscov2_pathway",
                                                    "(WP) Mirna_role_in_immune_response_in_sepsis",
                                                    "(WP) Hostpathogen_interaction_of_human_coronaviruses_interferon_induction",
                                                    "(WP) Mitochondrial_immune_response_to_sarscov2",
                                                    "(WP) Hostpathogen_interaction_of_human_coronaviruses_autophagy",
                                                    "(WP) Type_i_interferon_induction_and_signaling_during_sarscov2_infection",
                                                    "(WP) 4hydroxytamoxifen_dexamethasone_and_retinoic_acids_regulation_of_p27_expression",
                                                    "(WP) Glycosaminoglycan_degradation",
                                                    "(WP) Glycosylation_and_related_congenital_defects",
                                                    "(WP) Follicle_stimulating_hormone_fsh_signaling_pathway")))

p <- ggplot(df,aes(x=celltype ,y=Pathway_2,colour= NES,size= Pvalue005,shape = value))+
  geom_point(alpha = 0.8)+ 
  scale_size_continuous(range = c(1,5),breaks = seq(0,0.05,by = 0.01),
                        trans = "reverse")+
  theme_bw() + 
  xlab("")+ ylab("")+
  scale_x_discrete(labels =  c("B cells",
                               "T cells",
                               "NKT-like cells",
                               "Endothelial cells",
                               "Mural cells",
                               "Oligodendrocytes",
                               "Tumor cells",
                               "Macrophage-like GAMs",
                               "Microglia-like GAMs",
                               "Dendritic cells"))+
  scale_y_discrete(labels = rev(c("(WP) TLR4 signaling and tolerance",
                                  "(WP) LTF danger signal response pathway",
                                  "(WP) LDL influence on CD14 and TLR4",
                                  paste0("(WP) Quercetin and NF",expression(kappa),"B/AP1 induced apoptosis"),
                                  paste0("(WP) PPAR",expression(alpha)," pathway"),
                                  "(WP) AGE/RAGE pathway",
                                  "(WP) Degradation pathway of sphingolipids including diseases",
                                  "(WP) Transcription factor regulation in adipogenesis",
                                  "(PID) PDGFRB pathway", 
                                  "(PID) IL1 pathway",
                                  "(WP) IL1 signaling pathway",
                                  "(WP) IL3 signaling pathway",
                                  "(PID) IL6/IL7 pathway",
                                  "(PID) FCER1 pathway",
                                  "(PID) BCR signaling pathway",
                                  "(WP) Effects of NRF2 on SARSCoV-2 pathway",
                                  "(WP) MiRNA role in immune response in sepsis",
                                  "(WP) Hostpathogen interaction-coronaviruses interferon induction",
                                  "(WP) Mitochondrial immune response to SARSCoV-2",
                                  "(WP) Hostpathogen interaction-coronaviruses autophagy",
                                  "(WP) IFN induction and signaling during SARSCoV-2 infection",
                                  "(WP) Regulation of p27 expression",
                                  "(WP) Glycosaminoglycan degradation",
                                  "(WP) Glycosylation and related congenital defects",
                                  "(WP) Follicle stimulating hormone FSH signaling pathway")))+
  labs(size = "P value", shape = "Value")+ 
  scale_color_gradient2(midpoint=0, low="#2166AC", mid="#CCCCCC", high="#B2182B", limits = c(-1,1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 11,color = "black"),
        axis.text.y = element_text(size = 11,color = "black")) +guides(shape = FALSE, size = FALSE)

p 


TCGA <- read_xlsx("F:/Glioma/single cell (download)/OurPaperData/ssGSEA/TCGA_ssGSEA_25signal_OSPFS.xlsx")
CGGA <- read_xlsx("F:/Glioma/single cell (download)/OurPaperData/ssGSEA/CGGA_ssGSEA_25signal_OS.xlsx")
TCGA_pw <- c("WP_4_HYDROXYTAMOXIFEN_DEXAMETHASONE_AND_RETINOIC_ACIDS_REGULATION_OF_P27_EXPRESSION",
             "WP_AGE_RAGE_PATHWAY",                                                                
             "WP_ANTIVIRAL_AND_ANTI_INFLAMMATORY_EFFECTS_OF_NRF2_ON_SARS_COV_2_PATHWAY",           
             "WP_DEGRADATION_PATHWAY_OF_SPHINGOLIPIDS_INCLUDING_DISEASES",                         
             "WP_FOLLICLE_STIMULATING_HORMONE_FSH_SIGNALING_PATHWAY",                              
             "WP_GLYCOSAMINOGLYCAN_DEGRADATION",                                                   
             "WP_GLYCOSYLATION_AND_RELATED_CONGENITAL_DEFECTS",                                    
             "WP_HOST_PATHOGEN_INTERACTION_OF_HUMAN_CORONAVIRUSES_AUTOPHAGY",                      
             "WP_HOST_PATHOGEN_INTERACTION_OF_HUMAN_CORONAVIRUSES_INTERFERON_INDUCTION",           
             "WP_IL_1_SIGNALING_PATHWAY",                                                          
             "WP_IL_3_SIGNALING_PATHWAY",                                                          
             "WP_LDL_INFLUENCE_ON_CD14_AND_TLR4",                                                  
             "WP_LTF_DANGER_SIGNAL_RESPONSE_PATHWAY",                                              
             "WP_MIRNA_ROLE_IN_IMMUNE_RESPONSE_IN_SEPSIS",                                         
             "WP_MITOCHONDRIAL_IMMUNE_RESPONSE_TO_SARS_COV_2",                                     
             "WP_PPAR_ALPHA_PATHWAY",                                                              
             "WP_QUERCETIN_AND_NF_KB_AP_1_INDUCED_APOPTOSIS",                                      
             "WP_TLR4_SIGNALING_AND_TOLERANCE",                                                    
             "WP_TRANSCRIPTION_FACTOR_REGULATION_IN_ADIPOGENESIS",                                 
             "WP_TYPE_I_INTERFERON_INDUCTION_AND_SIGNALING_DURING_SARS_COV_2_INFECTION",           
             "PID_BCR_5PATHWAY",                                                                   
             "PID_FCER1_PATHWAY",                                                                  
             "PID_IL1_PATHWAY",                                                                    
             "PID_IL6_7_PATHWAY",                                                                  
             "PID_PDGFRB_PATHWAY" )
#for TCGA
df <- TCGA[,1]
df_OS <- data.frame()
df_PFS <- data.frame()

for (i in 1:length(TCGA_pw)) {
  data <- TCGA[str_subset(colnames(TCGA), TCGA_pw[i])]
  total_df <- cbind(df,data)
  total_df$Pathway <- TCGA_pw[i]
  total_df_OS <- total_df[-str_which(colnames(total_df), "PFS$")]
  total_df_OS$Label <- "TCGA-OS"
  total_df_PFS <- total_df[-str_which(colnames(total_df), "OS$")]
  total_df_PFS$Label <- "TCGA-PFS"
  OS <- paste0(TCGA_pw[i],"_Pvalue_OS")
  PFS <- paste0(TCGA_pw[i],"_Pvalue_PFS")
  
  total_df_OS_min <- total_df_OS[which.min(total_df_OS[,OS]),]
  total_df_PFS_min <- total_df_PFS[which.min(total_df_PFS[,PFS]),]
  
  colnames(total_df_OS_min) <- c("L_Number", "Pvalue", "Median_surv_H", "Median_surv_L", "pathway", "Label")
  colnames(total_df_PFS_min) <- c("L_Number", "Pvalue", "Median_surv_H", "Median_surv_L", "pathway", "Label")
  df_OS <- rbind(df_OS, total_df_OS_min)
  df_PFS <- rbind(df_PFS, total_df_PFS_min)
  print(TCGA_pw[i])
}

TCGA_df <- rbind(df_OS,df_PFS)



#for CGGA (NO PFS)
df <- CGGA[,1]
df_OS <- data.frame()
for (i in 1:length(TCGA_pw)) {
  data <- CGGA[str_subset(colnames(CGGA), TCGA_pw[i])]
  total_df <- cbind(df,data)
  total_df$Pathway <- TCGA_pw[i]
  total_df_OS <- total_df
  total_df_OS$Label <- "CGGA-OS"
  OS <- paste0(TCGA_pw[i],"_Pvalue_OS")
  
  
  total_df_OS_min <- total_df_OS[which.min(total_df_OS[,OS]),]
  
  
  colnames(total_df_OS_min) <- c("L_Number", "Pvalue", "Median_surv_H", "Median_surv_L", "pathway", "Label")
  df_OS <- rbind(df_OS, total_df_OS_min)
  print(TCGA_pw[i])
}
CGGA_df <- df_OS


merdata <- rbind(TCGA_df,CGGA_df)
merdata$Months <- merdata$Median_surv_H - merdata$Median_surv_L

merdata$value <- "<0.05"
merdata$value[which(merdata$Pvalue>=0.05)] <- ">=0.05"
merdata$Pvalue005 <- merdata$Pvalue
merdata$Pvalue005[which(merdata$Pvalue>=0.05)] <- 0.05
merdata$Label <- factor(merdata$Label, levels = c("TCGA-OS", "TCGA-PFS", "CGGA-OS"))
merdata$pathway <- factor(merdata$pathway, levels = rev(c("WP_TLR4_SIGNALING_AND_TOLERANCE",
                                                          "WP_LTF_DANGER_SIGNAL_RESPONSE_PATHWAY",
                                                          "WP_LDL_INFLUENCE_ON_CD14_AND_TLR4",
                                                          "WP_QUERCETIN_AND_NF_KB_AP_1_INDUCED_APOPTOSIS",
                                                          "WP_PPAR_ALPHA_PATHWAY", 
                                                          "WP_AGE_RAGE_PATHWAY",
                                                          "WP_DEGRADATION_PATHWAY_OF_SPHINGOLIPIDS_INCLUDING_DISEASES",
                                                          "WP_TRANSCRIPTION_FACTOR_REGULATION_IN_ADIPOGENESIS",
                                                          "PID_PDGFRB_PATHWAY",
                                                          "PID_IL1_PATHWAY",
                                                          "WP_IL_1_SIGNALING_PATHWAY",                                                          
                                                          "WP_IL_3_SIGNALING_PATHWAY",
                                                          "PID_IL6_7_PATHWAY",
                                                          "PID_FCER1_PATHWAY",
                                                          "PID_BCR_5PATHWAY",
                                                          "WP_ANTIVIRAL_AND_ANTI_INFLAMMATORY_EFFECTS_OF_NRF2_ON_SARS_COV_2_PATHWAY",
                                                          "WP_MIRNA_ROLE_IN_IMMUNE_RESPONSE_IN_SEPSIS",
                                                          "WP_HOST_PATHOGEN_INTERACTION_OF_HUMAN_CORONAVIRUSES_INTERFERON_INDUCTION",
                                                          "WP_MITOCHONDRIAL_IMMUNE_RESPONSE_TO_SARS_COV_2",
                                                          "WP_HOST_PATHOGEN_INTERACTION_OF_HUMAN_CORONAVIRUSES_AUTOPHAGY",
                                                          "WP_TYPE_I_INTERFERON_INDUCTION_AND_SIGNALING_DURING_SARS_COV_2_INFECTION",
                                                          "WP_4_HYDROXYTAMOXIFEN_DEXAMETHASONE_AND_RETINOIC_ACIDS_REGULATION_OF_P27_EXPRESSION",
                                                          "WP_GLYCOSAMINOGLYCAN_DEGRADATION",
                                                          "WP_GLYCOSYLATION_AND_RELATED_CONGENITAL_DEFECTS",
                                                          "WP_FOLLICLE_STIMULATING_HORMONE_FSH_SIGNALING_PATHWAY"
)))


p1 <- ggplot(merdata,aes(x=Label ,y=pathway,colour= Months,size= Pvalue005,shape = value))+
  geom_point(alpha = 0.8)+ 
  scale_size_continuous(range = c(1,5),breaks = seq(0,0.05,by = 0.01),
                        trans = "reverse")+
  theme_bw() + 
  xlab("")+ ylab("")+
  scale_x_discrete(labels = c("TCGA-OS", "TCGA-PFS", "CGGA-OS"))+
  scale_y_discrete(labels = rev(c("(WP) TLR4 signaling and tolerance",
                                  "(WP) LTF danger signal response pathway",
                                  "(WP) LDL influence on CD14 and TLR4",
                                  paste0("(WP) Quercetin and NF",expression(kappa),"B/AP1 induced apoptosis"),
                                  paste0("(WP) PPAR",expression(alpha)," pathway"),
                                  "(WP) AGE/RAGE pathway",
                                  "(WP) Degradation pathway of sphingolipids including diseases",
                                  "(WP) Transcription factor regulation in adipogenesis",
                                  "(PID) PDGFRB pathway", 
                                  "(PID) IL1 pathway",
                                  "(WP) IL1 signaling pathway",
                                  "(WP) IL3 signaling pathway",
                                  "(PID) IL6/IL7 pathway",
                                  "(PID) FCER1 pathway",
                                  "(PID) BCR signaling pathway",
                                  "(WP) Effects of NRF2 on SARSCoV-2 pathway",
                                  "(WP) MiRNA role in immune response in sepsis",
                                  "(WP) Hostpathogen interaction-coronaviruses interferon induction",
                                  "(WP) Mitochondrial immune response to SARSCoV-2",
                                  "(WP) Hostpathogen interaction-coronaviruses autophagy",
                                  "(WP) IFN induction and signaling during SARSCoV-2 infection",
                                  "(WP) Regulation of p27 expression",
                                  "(WP) Glycosaminoglycan degradation",
                                  "(WP) Glycosylation and related congenital defects",
                                  "(WP) Follicle stimulating hormone FSH signaling pathway")))+
  labs(size = "P value", shape = "Value")+
  scale_color_gradient2(midpoint=0, low="#CC6600", mid="white", high="#2D004B", limits = c(-10,10))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 11,color = "black"),
        #axis.text.y = element_text(size = 11,color = "black"
        axis.text.y = element_blank()
  ) 
p1
library(patchwork)
#patch plots a and b
patchwork1 <- p+p1
# Remove title and legend from second subplot
patchwork1[[2]] <- p