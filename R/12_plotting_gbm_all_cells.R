library(tidyverse)

#load RaceID and custom functions
source(file.path("R", "RaceID3_StemID2_class.R"))
source(file.path("R", "sankowski-et-al-functions.R"))

#load RaceID object
load(file.path("data", "sc_gbm_all_cells.RData"))

#please run this code if you would like to have the same t-SNE layout as the original paper
all_cells_nn <- read.csv(file.path("data", "gbm_all_cells_clusters+embeddings_nn.csv"))
sc@tsne <- all_cells_nn[, c("V1", "V2") ]

#define cluster order
if (!file.exists(file.path("data", "order_clusters_gbm_all_cells.RData")))  {
  order_clusters <-  clustheatmap(sc)
  save(order_clusters, file = file.path("data", "order_clusters_gbm_all_cells.RData"))
} else {
    load(file.path("data", "order_clusters_gbm_all_cells.RData"))
  }

# Figure 4a
df <- data.frame(ID=names(sc@cpart), Cluster= factor(sc@cpart, levels=order_clusters), sc@tsne)

#exclude small clusters
retain_cl <- as.numeric(names(table(df$Cluster))[table(df$Cluster) > dim(sc@ndata)[2]/100])
df <- df[df$Cluster %in% retain_cl,]

tsne_plot_no_outline(FILL=df$Cluster, fill_colors = colors_fig) +
  theme(legend.position = "none") 

# Extended data figure 9a - examination of cell types
signature_genes <- data.frame("microglia"= c('P2RY12', 'CX3CR1', 'CSF1R', 'TMEM119', 'SLC2A5'),
                              "monocytes"=c('CCR2', 'CLEC12A', 'PLAC8', 'FCN1', 'S100A9'),
                              "tcell"=c('TRAC', 'TRBC2', 'CD52', 'IL32', NA),
                              "oligodendrocyte"=c('MBP',  'MOG', 'MAG', 'PLP1', NA),
                              "macrophages"=c("MRC1", "MS4A7", "CD163", "LYVE1", "STAB1"),
                              "myeloid"=c('ITGAM',  'MS4A6A', 'TYROBP', 'CD14', NA),
                              "bcells"=c('CD79A', 'IGHG4', 'IGLL5', NA, NA),
                              "astrocyte"=c("GFAP", "HEPACAM","SOX9","AQP4",NA),
                              "apc"=c("CD74", "CD80", "CD86", "HLA-DRA", "CD40"), stringsAsFactors = F)

df <- data.frame(ID=names(sc@cpart), Cluster=sc@cpart, sc@tsne)

for (i in colnames(signature_genes)[1:4]) {
  tryCatch({
    pl <- plotexptsne2(name2id(na.omit(signature_genes[[i]]), rownames(sc@ndata)),.sc=sc, point_size = 4) +
      labs(subtitle = i)
    print(pl)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  on.exit(dev.off())
}     