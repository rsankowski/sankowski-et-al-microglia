#load packages
library(tidyverse)
library(pheatmap)
source(file.path("R", "sankowski-et-al-functions.R"))

# Figure 1g - deep neural network classification plot
load(file.path("data", "deep-learning-predictor-control", "jackel-predictions.Robj"))
load(file.path("data", "deep-learning-predictor-control", "masuda-predictions.Robj"))
load(file.path("data", "deep-learning-predictor-control", "mathys-predictions-healthy.Robj"))
load(file.path("data", "deep-learning-predictor-control", "schirmer-predictions-healthy.Robj"))
load(file.path("data", "deep-learning-predictor-control", "test-dataset-predictions.Robj"))
load(file.path("data", "deep-learning-predictor-control", "velmeshev-predictions-healthy.Robj"))
load(file.path("data", "deep-learning-predictor-control", "train-dataset-predictions.Robj"))

df <- bind_rows(
  data.frame(dataset="Train_set", predictions_train),
  data.frame(dataset="Test_set", predictions_test),
  data.frame(dataset="Masuda_et_al", predictions_masuda),
  data.frame(dataset="Velmeshev_et_al", predictions_aut_healthy),
  data.frame(dataset="Schirmer_et_al", predictions_schirmer),
  data.frame(dataset="Jakel_et_al", predictions_jackel),
  data.frame(dataset="Mathys_et_al", predictions_mathys)
)

df$dataset <- factor(df$dataset, levels = c("Train_set","Test_set","Masuda_et_al","Velmeshev_et_al","Schirmer_et_al","Jakel_et_al","Mathys_et_al"))

mosaicGG2(df, X="dataset", FILL = "Cluster", colors = colors_fig)

# Extended data figure 5b
load(file.path("data", "deep-learning-predictor-control", "train-dataset-predictions-verum-labels.Robj"))

mat <- as.matrix.data.frame(table(verum_labels[,3],verum_labels[,2]))
dimnames(mat) <- list(levels(verum_labels$Cluster_raceid), levels(verum_labels$Cluster_raceid))
mat2 <- apply(mat, 2, function(x) x/sum(x))
pheatmap(mat2, show_colnames = T, show_rownames = T, display_numbers = mat, fontsize = 20, font_col = "black", cluster_rows = F, cluster_cols = F)

