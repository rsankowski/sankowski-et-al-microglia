#figure 1g neural network classifier
#run keras
#Tutorial from: https://www.datacamp.com/community/tutorials/keras-r-deep-learning
#devtools::install_github("rstudio/keras")
#BiocManager::install("reticulate")
#devtools::install_github("rstudio/tensorflow")
# Load in the keras package
#install.packages("keras")
#devtools::install_github("rstudio/tfestimators")
#library(tfestimators)
library(reticulate)
library(keras)
library(tensorflow)
library(tidyverse)
library(tfestimators)

#in the shell i install python 3.7 and make it the default python. then i install anaconda and create the conda environment "r-tensorflow"
#conda create -n r-tensorflow python=3.7 anaconda

#use_condaenv("r-tensorflow")
#install_keras(method = "conda")
#install_tensorflow(method = "conda")

load(file.path("data", "sc_ctrl_microglia_nn.RData"))

genes2 <- gsub("_.*", "", rownames(sc@fdata))

#prepare data
train_matrix <- as.matrix(sc@fdata)
rownames(train_matrix) <- gsub('_.*', '', rownames(train_matrix))

#normalization steps are from this paper: https://www.biorxiv.org/content/biorxiv/early/2019/01/28/532093.full.pdf
cell_ids <- colnames(train_matrix)

train_matrix <- t(train_matrix) %>% normalize()
train_matrix <- log2(train_matrix*10000+1)
train_labels <- sc@cpart
train_data <- train_labels %>%
  cbind(train_matrix) 

dimnames(train_data) <- NULL

#split dataset url: https://www.r-bloggers.com/how-to-implement-random-forests-in-r/
set.seed(100)
train <- sample(nrow(train_data), 0.7*nrow(train_data), replace = FALSE)
x_train <- train_data[train,-1]
x_test <- train_data[-train,-1]

#split dataset url: https://www.r-bloggers.com/how-to-implement-random-forests-in-r/
y_train <- to_categorical(train_data[train,1])[,-1] #-1 see this url to understand why i remove the first column: https://github.com/rstudio/keras/issues/53
y_test <- to_categorical(train_data[-train,1])[,-1] #-1

model <- keras_model_sequential() 

#dimensions are taken from this paper: https://www.biorxiv.org/content/biorxiv/early/2019/01/28/532093.full.pdf
model <- model %>%
  layer_dense(units = 100, activation = 'relu', input_shape = c(792),kernel_regularizer = regularizer_l2(l = 0.005)) %>% 
  layer_dense(units = 50, activation = 'relu',kernel_regularizer = regularizer_l2(l = 0.005)) %>%
  layer_dense(units = 25, activation = 'relu',kernel_regularizer = regularizer_l2(l = 0.005)) %>% 
  layer_dense(units = 8, activation = 'softmax')

model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = 'adam',
  metrics = 'accuracy'
)  

model %>% fit(
  x_train, 
  y_train, 
  epochs = 500, 
  batch_size = 128, 
  validation_split = 0.2
)

plot(history)

model %>% save_model_hdf5("data/deep-learning-model-controls.h5")

model <- load_model_hdf5("data/deep-learning-model-controls.h5")

#evaluate model
#accuracy of the training set
predictions <- model %>% predict_classes(x_train)
predictions

table(train_data[train,1], predictions+1)
mean(train_data[train,1] == predictions+1)

#save predictions on the train dataset
predictions_train <- data.frame("cell_ID" = cell_ids[train], "Cluster" = predictions)
head(predictions_train)
predictions_train$Cluster <- as.factor(predictions_train$Cluster)
levels(predictions_train$Cluster) <- levels(train_labels)
save(predictions_train, file = "data/deep-learning-predictor-control/train-dataset-predictions.Robj")

verum_labels <- data.frame("cell_ID"=cell_ids[train], "Cluster_raceid"=train_data[train,1], "Cluster_predict" = predictions_train$Cluster)
verum_labels$Cluster_raceid <- as.factor(verum_labels$Cluster_raceid)
levels(verum_labels$Cluster_raceid) <- levels(train_labels)
save(verum_labels, file = "data/deep-learning-predictor-control/train-dataset-predictions-verum-labels.Robj")

#random cluster assignment, i.e. guessing
set.seed(79106)
predictions <- model %>% predict_classes(x_train[sample(1:nrow(x_train), nrow(x_train), replace = F),])
predictions

table(train_data[train,1], predictions+1)
mean(train_data[train,1] == predictions+1)


#accuracy of the test set
predictions <- model %>% predict_classes(x_test)
predictions

table(train_data[-train,1], predictions+1)
mean(train_data[-train,1] == predictions+1)

predictions_test <- data.frame("cell_ID" = cell_ids[-train], "Cluster" = predictions)
head(predictions_test)
predictions_test$Cluster <- as.factor(predictions_test$Cluster)
levels(predictions_test$Cluster) <- levels(train_labels)
save(predictions_test, file = "data/deep-learning-predictor-control/test-dataset-predictions.Robj")

##guessing on test set
predictions <- model %>% predict_classes(x_test[sample(1:nrow(x_test), nrow(x_test), replace = F),])
a <- rep({predict_classes(model, x_test[sample(1:nrow(x_test), nrow(x_test), replace = F),])+1}, 100)
predictions

table(train_data[-train,1], predictions+1)
mean(train_data[-train,1] == predictions+1)



#The script below tests the classification of independent human datasets, 
#as some of them have restricted access I only provide the output of the predictions in the data/deep-learning-predictor-control folder


if (F) {
  #autism data- healthy
  load(dataset)
  
  #load metadata
  meta <- read_tsv(metadata)
  micr_healthy <- which(meta$cluster == "Microglia" & meta$diagnosis == "Control") 
  
  #mean age
  meta %>% distinct(individual, age, diagnosis) %>% group_by(diagnosis) %>% summarise(mean_age=mean(age))
  
  ids <- meta$cell[micr_healthy]
  aut_healthy_data <- as.data.frame(as.matrix(micr_counts[, ids]))
  
  #extract conmon genes genes
  aut_healthy_data2 <- aut_healthy_data[genes2,]  
  sum(genes2 %in% rownames(aut_healthy_data2))
  
  #preprocess data
  #missing value imputation
  genes_na_ind <- which(!rownames(aut_healthy_data2) %in% genes2)
  
  aut_healthy_data2[is.na(aut_healthy_data2)] <- 0
  aut_healthy_data3 <- t(aut_healthy_data2) %>% normalize()
  aut_healthy_data3 <- log2(aut_healthy_data3*10000+1)
  genes_na_ind_exp <- matrix(rowMeans(train_matrix[genes_na_ind,]), nrow=(nrow(aut_healthy_data3)), ncol = length(genes_na_ind), byrow = TRUE) 
  
  aut_healthy_data3[,genes_na_ind] <- genes_na_ind_exp
  
  predictions_aut_healthy <- model %>% predict_classes(aut_healthy_data3)
  predictions_aut_healthy
  table(predictions_aut_healthy)
  levels(train_labels)
  predictions_aut_healthy <- as.factor(predictions_aut_healthy)
  levels(predictions_aut_healthy) <- levels(train_labels)
  
  #build data frame of predictions
  predictions_aut_healthy <- data.frame("cell_ID" = colnames(aut_healthy_data), "Cluster" = predictions_aut_healthy)
  head(predictions_aut_healthy)
  
  save(predictions_aut_healthy, file = "data/velmeshev-predictions-healthy.Robj")
  
  
  #schirmer et al data
  #autism data- healthy
  load(dataset)
  
  schirmer_data <- as.data.frame(as.matrix(micr_counts))
  
  #extract conmon genes genes
  schirmer_data2 <- schirmer_data[genes2,]  
  sum(genes2 %in% rownames(schirmer_data2))
  
  #preprocess data
  #missing value imputation
  genes_na_ind <- which(!rownames(schirmer_data2) %in% genes2)
  
  schirmer_data2[is.na(schirmer_data2)] <- 0
  schirmer_data3 <- t(schirmer_data2) %>% normalize()
  schirmer_data3 <- log2(schirmer_data3*10000+1)
  genes_na_ind_exp <- matrix(rowMeans(train_matrix[genes_na_ind,]), nrow=(nrow(schirmer_data3)), ncol = length(genes_na_ind), byrow = TRUE) 
  
  schirmer_data3[,genes_na_ind] <- genes_na_ind_exp
  
  predictions_schirmer <- model %>% predict_classes(schirmer_data3)
  predictions_schirmer
  table(predictions_schirmer)
  levels(train_labels)
  predictions_schirmer <- as.factor(predictions_schirmer)
  levels(predictions_schirmer) <- levels(train_labels)[c(2,4,6,7)]
  
  #build data frame of predictions
  predictions_schirmer <- data.frame("cell_ID" = colnames(schirmer_data), "Cluster" = predictions_schirmer)
  head(predictions_schirmer)
  
  save(predictions_schirmer, file = "data/schirmer-predictions-healthy.Robj")
  
  #mathys et al data
  load(dataset)
  micr_counts <- micr_counts[, colnames(micr_counts) %in% read_csv(metadata)[[1]]]
  
  mathys_data <- as.data.frame(as.matrix(micr_counts))
  
  #extract conmon genes genes
  mathys_data2 <- mathys_data[genes2,]  
  sum(genes2 %in% rownames(mathys_data2))
  
  #preprocess data
  #missing value imputation
  genes_na_ind <- which(!rownames(mathys_data2) %in% genes2)
  
  mathys_data2[is.na(mathys_data2)] <- 0
  mathys_data3 <- t(mathys_data2) %>% normalize()
  mathys_data3 <- log2(mathys_data3*10000+1)
  genes_na_ind_exp <- matrix(rowMeans(train_matrix[genes_na_ind,]), nrow=(nrow(mathys_data3)), ncol = length(genes_na_ind), byrow = TRUE) 
  
  mathys_data3[,genes_na_ind] <- genes_na_ind_exp
  
  predictions_mathys <- model %>% predict_classes(mathys_data3)
  predictions_mathys
  table(predictions_mathys)
  levels(train_labels)
  predictions_mathys <- as.factor(predictions_mathys)
  levels(predictions_mathys) <- levels(train_labels)
  
  #build data frame of predictions
  predictions_mathys <- data.frame("cell_ID" = colnames(mathys_data), "Cluster" = predictions_mathys)
  head(predictions_mathys)
  table(predictions_mathys$Cluster)
  
  save(predictions_mathys, file = "data/mathys-predictions-healthy.Robj")
  
  #masuda et al dataset
  load(dataset)
  
  meta_masuda <- read_csv(metadata)
  
  masuda <- Matrix(as.matrix(sc@expdata[,colnames(sc@expdata) %in% meta_masuda$cell_ID]), sparse=TRUE)
  rownames(masuda) <- gsub("_.*", "", rownames(masuda))
  masuda_data2 <- as.data.frame(as.matrix(masuda))[genes2,]  
  sum(genes2 %in% rownames(masuda_data2))
  
  #preprocess data
  #missing value imputation
  genes_na_ind <- which(!rownames(masuda_data2) %in% genes2)
  
  masuda_data2[is.na(masuda_data2)] <- 0
  masuda_data3 <- t(masuda_data2) %>% normalize()
  masuda_data3 <- log2(masuda_data3*10000+1)
  genes_na_ind_exp <- matrix(rowMeans(train_matrix[genes_na_ind,]), nrow=(nrow(masuda_data3)), ncol = length(genes_na_ind), byrow = TRUE) 
  
  masuda_data3[,genes_na_ind] <- genes_na_ind_exp
  
  predictions_masuda <- model %>% predict_classes(masuda_data3)
  predictions_masuda
  table(predictions_masuda)
  levels(train_labels)
  predictions_masuda <- as.factor(predictions_masuda)
  levels(predictions_masuda) <- levels(train_labels)
  
  #build data frame of predictions
  predictions_masuda <- data.frame("cell_ID" = colnames(masuda), "Cluster" = predictions_masuda)
  head(predictions_masuda)
  table(predictions_masuda$Cluster)
  
  save(predictions_masuda, file = "data/masuda-predictions.Robj")
  
  labels_masuda <-predictions_masuda %>%  bind_cols(data.frame("cluster_paper" = sc@cpart[meta_masuda$cell_ID]))
  
  table(labels_masuda$Cluster, labels_masuda$cluster_paper)
  
  save(labels_masuda, file = "data/masuda-labels-with-deep-learning-predicted-labels.Robj")
  
  #jackel et al
  url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118257&format=file&file=GSE118257%5FMSCtr%5FsnRNA%5FExpressionMatrix%5FR%2Etxt%2Egz'
  tmp <- tempfile()
  ##
  download.file(url,tmp)
  prdata <- read.delim(gzfile(tmp), sep = '\t', stringsAsFactors = F)
  
  #load metadata
  url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118257&format=file&file=GSE118257%5FMSCtr%5FsnRNA%5FFinalAnnotationTable%2Etxt%2Egz'
  tmp <- tempfile()
  ##
  download.file(url,tmp)
  metadata <- read.delim(gzfile(tmp), sep = '\t', stringsAsFactors = F)
  colnames(metadata)[1] <- 'ID'
  micr <- metadata$ID[metadata$Lesion %in% c("Ctrl","NAWM") & metadata$Celltypes == "Microglia_Macrophages"]
  micr <- gsub(":", "\\.", micr)
  
  
  jackel_data2 <- as.data.frame(as.matrix(prdata))[genes2,colnames(prdata) %in% micr]  
  sum(genes2 %in% rownames(jackel_data2))
  
  #preprocess data
  #missing value imputation
  genes_na_ind <- which(!rownames(jackel_data2) %in% genes2)
  
  jackel_data2[is.na(jackel_data2)] <- 0
  jackel_data3 <- t(jackel_data2) %>% normalize()
  jackel_data3 <- log2(jackel_data3*10000+1)
  genes_na_ind_exp <- matrix(rowMeans(train_matrix[genes_na_ind,]), nrow=(nrow(jackel_data3)), ncol = length(genes_na_ind), byrow = TRUE) 
  
  jackel_data3[,genes_na_ind] <- genes_na_ind_exp
  
  predictions_jackel <- model %>% predict_classes(jackel_data3)
  predictions_jackel
  table(predictions_jackel)
  levels(train_labels)
  predictions_jackel <- as.factor(predictions_jackel)
  levels(predictions_jackel) <- c("1","5","3", "2")
  
  #build data frame of predictions
  predictions_jackel <- data.frame("cell_ID" = colnames(prdata)[colnames(prdata) %in% micr], "Cluster" = predictions_jackel)
  head(predictions_jackel)
  
  save(predictions_jackel, file = "data/jackel-predictions.Robj")
  
}

