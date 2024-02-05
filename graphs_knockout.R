rm(list = ls(all = TRUE)); graphics.off()
setwd("//fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/data backup/GBM_Screens/Broad screening/Knockout_synthetic_lethal_raw_data")
#setwd("/Volumes/groups/Neuro-Oncology/Clinical Neuro-Oncology/data backup/GBM_Screens/Broad screening/Knockout_synthetic_lethal_raw_data")
library(MAGeCKFlute)
library(ggplot2)
library(tidyverse)
library(dplyr)
source("ScatterViewCustom.R")
options(ggrepel.max.overlaps = Inf)



#################PARAMETERS YOU CAN MODIFY#############################################
#please put "2. batch", "1. batch" "3. batch"
batch_number <- "1. batch"

#please put "GS.9", "LN229", "LN18", "T98G"
cell_line <- "T98G"
#
drug <- "TMZ"
gene_list <- c("NF1", "MAP2K1")


#######################################################################################
essential.genes = read.delim("CEGv2.txt", header = T)

#####I am going in the folder of the cell line and drug you specified##################

batch <- paste0("//fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/data backup/GBM_Screens/Broad screening/Knockout_synthetic_lethal_raw_data/", batch_number,"/Results/")



if(batch_number != "1. batch") {
setwd(paste0(batch,cell_line,"_",drug))
  } else if (batch_number != "grouped_cell_lines")
  {
    setwd(paste0(batch,cell_line,"_",drug, "_vs_DMSO"))
  } else
  {
    setwd(batch)
  }
########################################################################################



#Reading the 4 tables#################################################################
if(batch_number != "1. batch") {
Drugvsplasmid <- ReadBeta(paste0("plasmid/",cell_line, "_",drug,"_vs_plasmid.gene_summary.txt"))
DrugvsDMSO <- ReadBeta(paste0(cell_line, "_",drug,"_vs_DMSO.gene_summary.txt"))
DMSOvsplasmid <- ReadBeta(paste0(batch,cell_line, "_DMSO/",cell_line,"_DMSO_vs_plasmid.gene_summary.txt"))
DMSOvsplasmid.volcano = read.delim(paste0(batch,cell_line, "_DMSO/",cell_line,"_DMSO_vs_plasmid.gene_summary.txt"))
} else if(batch_number != "grouped_cell_lines")
{ 
  DMSOvsplasmid <- ReadBeta(paste0(batch,cell_line, "_DMSO_vs_plasmid/",cell_line,"_DMSO_vs_plasmid.gene_summary.txt"))
  DrugvsDMSO <- ReadBeta(paste0(cell_line, "_",drug,"_vs_DMSO.gene_summary.txt"))
  Drugvsplasmid <- ReadBeta(paste0(batch, cell_line, "_",drug,"_vs_plasmid/",cell_line, "_",drug,"_vs_plasmid.gene_summary.txt"))
  DMSOvsplasmid.volcano = read.delim(paste0(batch,cell_line, "_DMSO_vs_plasmid/",cell_line,"_DMSO_vs_plasmid.gene_summary.txt"))
  
} else{
  DMSOvsplasmid <- ReadBeta(paste0(batch,cell_line, "_DMSO_vs_plasmid/",cell_line,"_DMSO_vs_plasmid.gene_summary.txt"))
  DrugvsDMSO <- ReadBeta(paste0(cell_line, "_",drug,"_vs_DMSO.gene_summary.txt"))
  Drugvsplasmid <- ReadBeta(paste0(batch, cell_line, "_",drug,"_vs_plasmid/",cell_line, "_",drug,"_vs_plasmid.gene_summary.txt"))
  DMSOvsplasmid.volcano = read.delim(paste0(batch,cell_line, "_DMSO_vs_plasmid/",cell_line,"_DMSO_vs_plasmid.gene_summary.txt"))
}

#Merge Drug vs Plasmid and DMSO vs Plasmid
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
x_var <- colnames(DMSOvsplasmid.volcano)[3]
y_var <- colnames(DMSOvsplasmid.volcano)[7]

colnames(essential.genes)[1] = "Gene"
DMSOvsplasmid.volcano %>%
  ggplot(aes_string(x = x_var, y = paste0("-log10(", y_var, ")"), label = "Gene")) +
  geom_point()+
  geom_point(data = DMSOvsplasmid.volcano %>% inner_join(essential.genes, by = "Gene"), 
             aes(label = Gene), color = "red") +
  theme_classic()



#Scatter View######
Scatter.Activation = ScatterViewCustom(activation.Scatter, "DMSO", drug, 
                                           label = "Gene",size =2, groups = c("top", "bottom"), top=0,
                                           toplabels = gene_list,
                                           auto_cut_diag  = 2, group_col = c("#e41a1c", "#377eb8"),
                                           alpha = 0.5, display_cut = TRUE)
print(Scatter.Activation + theme_classic() + theme(legend.position="none"))


write.table(Scatter.Activation[["data"]], 
            file = paste0(drug, "_", cell_line, "_scatter_data.tsv"), 
            sep="\t", quote = F, col.names = T, row.names = F)


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
                                           main = paste0("Rank Activation ", cell_line, " ", drug, "vs DMSO"),
                                           width = 2, height = 4)
print(Rank.Activation.DrugvsDMSO + theme_classic() +theme(legend.position="none"), vp=grid::viewport(gp=grid::gpar(cex=2)))


#Rank without genes
Rank.Activation.DrugvsDMSO.blank = RankView(genelist.DrugvsDMSO,
                                                 top = 0, bottom = 0,
                                                 main = paste0("Rank Activation ", cell_line, " ", drug, "vs DMSO"),
                                                 width = 2, height = 4)
print(Rank.Activation.DrugvsDMSO.blank + theme_classic() +theme(legend.position="none"), vp=grid::viewport(gp=grid::gpar(cex=2)))
####################

###########################################################################################


