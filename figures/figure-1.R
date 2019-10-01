#figure 1
source("bin/RaceID3_StemID2_class.R")
source("bin/sankowski-et-al-functions.R")
library(tidyverse)
library(viridis)
library(reshape2)

#1b
    load("data/df_control_all_cells.RData")
    tsne_plot_no_outline(df_all, FILL = df$Cluster, fill_colors = colors_fig, point_outline = colors_fig, point_size = 3, line_width = 0.25)
    

#1c
    #load data
    load("data/counts_control.RData")
    load("data/df_control.RData")
    
    #plot heatmap
    plot_heatmap(.sc = counts_control, .df=df_control, .up_genes= up_genes_control, .colors = viridis(100))
    
#1d
    