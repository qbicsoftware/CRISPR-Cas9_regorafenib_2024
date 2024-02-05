rm(list = ls(all = TRUE)); graphics.off()
#setwd("/Volumes/groups/Neuro-Oncology/Clinical Neuro-Oncology/data backup/GBM_Screens/Broad screening/Knockout_synthetic_lethal_raw_data")
setwd("//fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/data backup/GBM_Screens/Broad screening/Activation_resistance_raw_data")
library(MAGeCKFlute)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
source("ScatterViewCustom.R")
options(ggrepel.max.overlaps = Inf)
#essential.genes = read.delim("CEGv2.txt", header = T)


####please put one of the following : Abema   CCNU  Evero   Rego  TMZ  VAL
drug <- "VAL"
gene_list <- c("BCL2", "BCL2L1")

####please put one of the following : LN229_GS.9_pilote   LN229_GS.9_3  LN229_GS.9    LN18_T98G    
#cell_lines 

##########
#batch <- paste0(
 # "/Volumes/groups/Neuro-Oncology/Clinical Neuro-Oncology/data backup/GBM_Screens/Broad screening/Knockout_synthetic_lethal_raw_data/grouped_cell_lines/mageck/mle/")
batch <- paste0(
  "//fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/data backup/GBM_Screens/Broad screening/Activation_resistance_raw_data/results_grouped_activation/mageck/mle/")
Drugvsplasmid <- ReadBeta(paste0(batch,drug,"_vs_plasmid/",drug,"_vs_plasmid.gene_summary.txt"))
DrugvsDMSO <- ReadBeta(paste0(batch,drug,"_vs_DMSO/",drug,"_vs_DMSO.gene_summary.txt"))
DMSOvsplasmid <- ReadBeta(paste0(batch,"DMSO_vs_plasmid_",str_to_lower(drug),"/DMSO_vs_plasmid_",str_to_lower(drug),".gene_summary.txt"))


activation.Scatter = Drugvsplasmid  %>% left_join(DMSOvsplasmid)

####Rename columns to remove the cell_line name
colnames(activation.Scatter)[2] <- drug
colnames(activation.Scatter)[3] = "DMSO"
colnames(DrugvsDMSO)[2] <- drug
#######################################################################################

#######################################################################################
#This is to get the data for the rank plots in an object that RankView is happy with
genelist.DrugvsDMSO <- DrugvsDMSO[[drug]]
names(genelist.DrugvsDMSO) <- DrugvsDMSO$Gene
#######################################################################################

###############GRAPHS##################################################################


#Scatter View######
Scatter.Activation = ScatterViewCustom(activation.Scatter, "DMSO", drug, 
                                       label = "Gene",size =2, groups = c("top", "bottom"), top=0,
                                       toplabels = gene_list,
                                       auto_cut_diag  = 2, group_col = c("#e41a1c", "#377eb8"),
                                       alpha = 0.5, display_cut = TRUE)
print(Scatter.Activation + theme_classic() + theme(legend.position="none"))


#Scatter View without genes
Scatter.Activation.blank = ScatterViewCustom(activation.Scatter, "DMSO", drug, 
                                             label = "Gene",size =2, groups = c("top", "bottom"), top=0,
                                             auto_cut_diag  = 2, group_col = c("#e41a1c", "#377eb8"),
                                             alpha = 0.5, display_cut = TRUE)
print(Scatter.Activation.blank + theme_classic() + theme(legend.position="none"))
####################

#9square
Scatter.9sq = ScatterViewCustom(activation.Scatter, "DMSO", drug, 
                                label = "Gene",size =2,groups = c("top", "bottom"), top=0,
                                model = c("ninesquare"),
                                auto_cut_diag  = 2, group_col = c("#e41a1c", "#377eb8"),
                                alpha = 0.5, display_cut = TRUE)
print(Scatter.9sq + theme_classic() + theme(legend.position="none"))
####################

#Rank
Rank.Activation.DrugvsDMSO = RankView(genelist.DrugvsDMSO,
                                      top = 0, bottom = 0,
                                      genelist = gene_list,
                                      main = paste0("Rank Activation ", " ", drug, "vs DMSO"),
                                      width = 2, height = 4)
print(Rank.Activation.DrugvsDMSO + theme_classic() +theme(legend.position="none"), vp=grid::viewport(gp=grid::gpar(cex=2)))


#Rank without genes
Rank.Activation.DrugvsDMSO.blank = RankView(genelist.DrugvsDMSO,
                                            top = 0, bottom = 0,
                                            main = paste0("Rank Activation ", " ", drug, "vs DMSO"),
                                            width = 2, height = 4)
print(Rank.Activation.DrugvsDMSO.blank + theme_classic() +theme(legend.position="none"), vp=grid::viewport(gp=grid::gpar(cex=2)))
####################

###########################################################################################

