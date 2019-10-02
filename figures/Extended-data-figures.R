#Extended data figures
#Extended Data Figure 3 
                #plot gene signatures
                cell_signatures <- data.frame("monocytes"=c('CCR2', 'CLEC12A', 'PLAC8', 'FCN1', 'S100A9'),
                                              "microglia"= c('P2RY12', 'CX3CR1', 'CSF1R', 'TMEM119', 'SLC2A5'),
                                              "tcells"=c('TRAC', 'TRBC2', 'CD52', 'IL32', NA),
                                              "myeloid_cells"=c('ITGAM',  'MS4A6A', 'TYROBP', 'CD14', NA),
                                              "oligodendrocytes"=c('MBP',  'MOG', 'MAG', 'PLP1', NA))
                
                for (i in colnames(cell_signatures)) {
                  
                  plt <- plotexptsne2(gene=na.omit(cell_signatures[[i]]) ,.sc=counts_all, point_size=3,line_width=0) + labs(title=i) 
                  print(plt)
                  
                }                 
                

#Extended Data Figure 5b
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
                
                #heatmap of test set predictions
                load("data/deep-learning-predictor-control/train-dataset-predictions-verum-labels.Robj")
                cont_tab <- as.matrix(table(verum_labels$Cluster_raceid, verum_labels$Cluster_predict))
                
                pheat <- pheatmap(cont_tab, display_numbers=T,cluster_rows = F, cluster_cols = F, cellwidth = 40, cellheight = 40, number_format = "%.0f", fontsize_number = 15, number_color = "black", angle_col = 0)
                
                pheat <- pheatmap(prop.table(cont_tab,1), cluster_rows = F, cluster_cols = F,display_numbers=T, cellwidth = 40, cellheight = 40, number_format = "%.0f", fontsize_number = 11, number_color = "black", angle_col = 0)

