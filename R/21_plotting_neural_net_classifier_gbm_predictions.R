library(tidyverse)
library(pheatmap)
library(viridis)

date = Sys.Date()

#load RaceID and custom functions
source(file.path("R", "RaceID3_StemID2_class.R"))
source(file.path("R", "sankowski-et-al-functions.R"))

#load data
load(file.path("data", "order_clusters_gbm_micr.RData"))

# Supplementary Figure 3a
#heatmap of test set predictions
load(file.path("data","train-dataset-predictions-verum-labels.Robj"))
cont_tab <- as.matrix(table(verum_labels$Cluster_raceid, verum_labels$Cluster_predict))
pheat <- pheatmap(prop.table(cont_tab,2), cluster_rows = F, cluster_cols = F,display_numbers= cont_tab, cellwidth = 20, cellheight = 20, number_format = "%.0f", fontsize_number = 11, number_color = "black", angle_col = 0)

# Supplementary Figure 3b
df2 <- list()
index=1
datasets <- c("training_set", "test_set", "darmanis-et-al") 
files <- grep("predictions", list.files("data"), value = T)

for (i in files) {
  load(file.path("data",i))
}

df2[datasets] <- list(predictions_train, predictions_test,predictions_gbm) #, predictions_aut_asd, predictions_oligo, 
df_long <- bind_rows(df2, .id="id")
df_long$id <- factor(df_long$id, levels = datasets)
df_long$Cluster <- factor(df_long$Cluster, levels = order_clusters)

#marimekko plot
mosaicGG2(df_long, "id", "Cluster", colors_fig, rect_col = 'black', line_width = 0.1) 