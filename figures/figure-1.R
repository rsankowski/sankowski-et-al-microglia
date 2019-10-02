#figure 1
source("bin/RaceID3_StemID2_class.R")
source("bin/sankowski-et-al-functions.R")
library(tidyverse)
library(viridis)
library(reshape2)

#1b
    load("data/df_control_all_cells.RData")
    tsne_plot_no_outline(df_all, FILL = df_all$Cluster, fill_colors = colors_fig, point_outline = colors_fig, point_size = 3, line_width = 0.25)
    

#1c
    #load data
    load("data/counts_control.RData")
    load("data/df_control.RData")
    
    #plot heatmap
    plot_heatmap(.sc = counts_control, .df=df_control, .up_genes= up_genes_control, .colors = viridis(100))
    
#1d
    #marimekko plot - patients
    mosaicGG2(df_control, "Cluster", "anon_ID", c(colors_pat, colors_many), rect_col = 'black', line_width = 0.1)
    
#1e
    tsne_plot_no_outline(df_control, FILL = df_control$Cluster, fill_colors = colors_fig, point_outline = colors_fig, point_size = 3, line_width = 0.25)

#1f
    pie_chart(.df=df_control)

#1g
    df <- list()
    index=1
    datasets <- c("training_set", "test_set", "masuda-et-al", "velmeshev-et-al", "schirmer-et-al", "mathys-et-al", "jackel-et-al") #,"aut_asd", "tirosh", "venteicher", "darmanis"
    files <- grep("predictions", list.files("data/deep-learning-predictor-control/"), value = T)
    
    for (i in files) {
        load(file.path("data/deep-learning-predictor-control",i))
    }
    
    df[datasets] <- list(predictions_train, predictions_test, predictions_masuda, predictions_aut_healthy, predictions_schirmer, predictions_mathys, predictions_jackel) #, predictions_aut_asd, predictions_oligo, predictions_astro, predictions_gbm
    
    df_long <- bind_rows(df, .id="id")
    df_long$id <- factor(df_long$id, levels = datasets)
    df_long$Cluster <- factor(df_long$Cluster, levels = c("8", "1", "9", "5", "6", "3", "2", "7"))
    
    #marimekko plot
    mosaicGG2(df_long, "id", "Cluster", colors_fig, rect_col = 'black', line_width = 0.1) 
    
    #marimekko stat plot
    mosaicGG(df_long,"id", "Cluster", rect_col = 'black', line_width = 0.1)

#1h
    genes <- (name2id(c("CX3CR1", "HLA-DRA", "^SPP1", "CCL2"), rownames(counts_control)))
    genes <- as.list(genes)
    lapply(genes, function(x) plotexptsne2(gene=x, .df = df_control, .sc = counts_control, point_size=3, line_width = 0))

#1i
    