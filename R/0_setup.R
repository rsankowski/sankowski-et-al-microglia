#install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

#the following code was adjusted from url https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
list.of.packages <- c(
          "caTools",
          "som",
          "epiR",
          "TSP",
          "igraph",
          "proxy",
          "network",
          "quadprog",
          "hypervolume",
          "rgl",
          "randomForest",
          "tsne",
          "pheatmap",
          "MASS",
          "cluster",
          "mclust",
          "flexmix",
          "lattice",
          "fpc",
          "amap",
          "RColorBrewer",
          "locfit",
          "vegan",
          "Rtsne",
          "scran",
          "DESeq2",
          "ica",
          "som",
          "caTools",
          "tidyverse",
          "assertthat",
          "viridis",
          "readxl",
          "reticulate",
          "keras",
          "tensorflow",
          "tfestimators",
          "clusterProfiler",
          "org.Hs.eg.db"
  )


new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

#before running tensorflow for the deep neural network predictor, please run:
library(tensorflow)
install_tensorflow() #see url for help: https://tensorflow.rstudio.com/installation/
library(keras)
install_keras() # see url for help: https://cran.r-project.org/web/packages/keras/vignettes/index.html

#these packages should download without issues. In case you get error messages please copy these and look them up at
#stackoverflow or discuss them with a bioinformaticion/bioinformaticist at your department. Also, don't hesitate to contact me directly.


