library(readxl)
library(pheatmap)
library(tidyverse)

#load RaceID and custom functions
source(file.path("R", "RaceID3_StemID2_class.R"))
source(file.path("R", "sankowski-et-al-functions.R"))

#load sc data
load(file.path("data", "sc_gbm_microglia_nn.RData"))
load(file.path("data", "order_clusters_gbm_micr.RData"))
load(file.path("data", "retain_cl_gbm_micr.RData"))
load(file.path("data", "metadata_ctrl_gbm.RData"))

# Figure 6e 
#load data
cytof_data <- read_excel(file.path("data", "cytof_expression.xlsx"))[,-27] %>% as.data.frame()

cytof_data <- cytof_data[cytof_data$`t-test` <.05,]
colnames(cytof_data) <- gsub(".fcs", "", colnames(cytof_data))
rownames(cytof_data) <- cytof_data[,1]
cytof_data[,26] <- NULL

#plot
pheatmap(cytof_data[,-1],
         scale = 'row', 
         color = colorRampPalette(rev(toupper(c('#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b'))))(100),
         cellwidth = 12,
         cellheight = 12)

# Supplementary Figure 8
# t-SNE plots
genes <- c("APOE", "FCGR1A", "FCGR3A", "TREM2", "CD86", "ADGRG1", "HLA-DRA", "CSF1R", "CX3CR1", "IRF8", "P2RY12", "TMEM119", "CD14", "ITGAX", "FCGR2A", "ADGRE1",
           "^CSF2", "PTPRC", "CD68", "CD44", "CSF2RA", "CD47", "^IL6", "^LRP1", "CD101", "^ICAM1", "CCNB1", "ITGAM", "IL10", "CD33", "CD274", "CCL4", "CCL2", "FLT3",
           "CCR7", "TGFB1", "^FAS", "CCR5", "TNF", "SIRPA", "PDCD1", "MRC1", "MKI67", "ITGB2", "CD163", "CD74", "SLC2A5")

walk(genes, function(i) print(plotexptsne2( .gene=i, point_size = 5)))

#line plots
genes2 <- name2id(genes, rownames(sc@ndata))
data_long <- make_data_long(.genes = genes2)
walk(gsub("(-|\\|)", ".",genes2), function(i) print(gene_line_plot(.gene=i)))
