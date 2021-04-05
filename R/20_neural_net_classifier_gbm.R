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

#load RaceID and custom functions
source(file.path("R", "RaceID3_StemID2_class.R"))
source(file.path("R", "sankowski-et-al-functions.R"))

#load data
load(file.path("data", "sc_gbm_microglia_nn.RData"))
load(file.path("data", "order_clusters_gbm_micr.RData"))
load(file.path("data", "retain_cl_gbm_micr.RData"))
load(file.path("data", "metadata_ctrl_gbm.RData"))

metadata <- df
genes <- sc@fdata %>% rownames
genes2 <- gsub("_.*", "", genes)

#prepare data
train_matrix <- as.matrix(sc@ndata[genes,metadata$cell_ID])

#normalization steps are from this paper: https://www.biorxiv.org/content/biorxiv/early/2019/01/28/532093.full.pdf
cell_ids <- colnames(sc@ndata[genes,metadata$cell_ID])

train_matrix <- t(train_matrix) %>% normalize()
train_matrix <- log2(train_matrix*10000+1)
train_labels <- metadata[["Cluster"]]
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
  layer_dense(units = 100, activation = 'relu', input_shape = c(463),kernel_regularizer = regularizer_l2(l = 0.005)) %>% 
  layer_dense(units = 50, activation = 'relu',kernel_regularizer = regularizer_l2(l = 0.005)) %>%
  layer_dense(units = 25, activation = 'relu',kernel_regularizer = regularizer_l2(l = 0.005)) %>% 
  layer_dense(units = 14, activation = 'softmax')

model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = 'adam',
  metrics = 'accuracy'
)  

model %>% fit(
  x_train, 
  y_train, 
  epochs = 300, 
  batch_size = 128, 
  validation_split = 0.2
)

plot(history)

# save model
model %>% save_model_hdf5(file.path("data", "deep-learning-model-gbm-data.h5"))

# load model
model <- load_model_hdf5(file.path("data", "deep-learning-model-gbm-data.h5"))

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
save(predictions_train, file = file.path("data", "train-dataset-predictions.Robj"))

verum_labels <- data.frame("cell_ID"=cell_ids[train], "Cluster_raceid"=train_data[train,1], "Cluster_predict" = predictions_train$Cluster)
verum_labels$Cluster_raceid <- as.factor(verum_labels$Cluster_raceid)
levels(verum_labels$Cluster_raceid) <- levels(train_labels)
save(verum_labels, file = file.path("data","train-dataset-predictions-verum-labels.Robj"))

#random guessing
set.seed(1234)
predictions <- model %>% predict_classes(x_train[sample(1:nrow(x_train), nrow(x_train), replace = F),])

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
save(predictions_test, file = file.path("data", "test-dataset-predictions.Robj"))

#get auc curve from url: https://blogs.rstudio.com/tensorflow/posts/2018-01-11-keras-customer-churn/
yhat_keras_class_vec <- predict_classes(object = model, x = x_test) %>%
  as.vector()

#test set gbm only
df2 <- df[-train,]
gbms <- which(df2$Diagnosis =="GBM")

#accuracy of the test set
x_test_gbm <- x_test[gbms,]
y_test_gbm <- y_test[gbms]

predictions <- model %>% predict_classes(x_test_gbm)
predictions

table(train_data[-train,1][gbms], predictions+1)
mean(train_data[-train,1][gbms] == predictions+1)

predictions_test_gbm <- data.frame("cell_ID" = cell_ids[-train][gbms], "Cluster" = predictions)
head(predictions_test_gbm)
predictions_test_gbm$Cluster <- as.factor(predictions_test_gbm$Cluster)
levels(predictions_test_gbm$Cluster) <- levels(train_labels)
save(predictions_test_gbm, file = "data/test-dataset-predictions-gbms.Robj")

#get auc curve from url: https://blogs.rstudio.com/tensorflow/posts/2018-01-11-keras-customer-churn/
yhat_keras_class_vec <- predict_classes(object = model, x = x_test) %>%
  as.vector()


# Predicted Class Probability
yhat_keras_prob_vec  <- predict_proba(object = model, x = x_test)

estimates_keras_tbl <- tibble(
  truth      = train_data[-train,1],
  estimate   = predictions_test$Cluster,
  class_prob = yhat_keras_prob_vec
)

estimates_keras_tbl$truth <- as.factor(estimates_keras_tbl$truth)
levels(estimates_keras_tbl$truth) <- levels(predictions_test$Cluster)
estimates_keras_tbl

#accuracy metrics
library(yardstick)
options(yardstick.event_first = FALSE)

# Confusion Table
estimates_keras_tbl %>% conf_mat(truth, estimate)

# Accuracy
estimates_keras_tbl %>% metrics(truth, estimate)

# AUC
estimates_keras_tbl %>% roc_auc(truth, class_prob)

# F1-Statistic
estimates_keras_tbl %>% f_meas(truth, estimate, beta = 1)


#precision recall
tibble(
  precision = estimates_keras_tbl %>% precision(truth, estimate),
  recall    = estimates_keras_tbl %>% recall(truth, estimate)
)

#Classify the darmanis dataset
gbm_data <- x[,grepl("^X", colnames(x))]

metadata <- load("/home/roman/Documents/Single cell analysis/Darmanis et. Barres dataset/RaceID4/darmanis-metadata.Robj")

#exclude peripheral cells 
periph_cells <- paste0("X",df$ID[df$Cluster %in% c(9,3,4,11)])
gbm_data <- gbm_data[,!colnames(gbm_data) %in% periph_cells]
#extract conmon genes genes
gbm_data2 <- gbm_data[genes2,]  #
sum(genes2 %in% rownames(gbm_data2))

#preprocess data
#missing value impution
genes_na_ind <- which(!rownames(gbm_data2) %in% genes2)

gbm_data2[is.na(gbm_data2)] <- 0
gbm_data3 <- t(gbm_data2) %>% normalize()
gbm_data3 <- log2(gbm_data3*10000+1)
genes_na_ind_exp <- matrix(rowMeans(train_matrix[genes_na_ind,]), nrow=(nrow(gbm_data3)), ncol = length(genes_na_ind), byrow = TRUE) 

gbm_data3[,genes_na_ind] <- genes_na_ind_exp

predictions_gbm <- model %>% predict_classes(gbm_data3)
predictions_gbm
table(predictions_gbm)
vals <- as.numeric(names(table(predictions_gbm))) +1
levels(train_labels)
predictions_gbm <- as.factor(predictions_gbm)
levels(predictions_gbm) <- levels(train_labels)[vals]

#build data frame of predictions
predictions_gbm <- data.frame("cell_ID" = colnames(gbm_data), "Cluster" = predictions_gbm)
head(predictions_gbm)
table(predictions_gbm$Cluster)

save(predictions_gbm, file = "data/darmanis-predictions.Robj")
