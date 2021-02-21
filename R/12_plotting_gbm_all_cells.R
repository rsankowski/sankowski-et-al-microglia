library(tidyverse)

#load RaceID and custom functions
source(file.path("R", "RaceID3_StemID2_class.R"))
source(file.path("R", "sankowski-et-al-functions.R"))

#load RaceID object
load(file.path("data", "sc_gbm_all_cells.RData"))

#please run this code if you would like to have the same t-SNE layout as the original paper
all_cells_nn <- read.csv(file.path("data", "gbm_all_cells_clusters+embeddings_nn.csv"))
sc@tsne <- all_cells_nn[, c("V1", "V2") ]

# Figure 9a
df <- data.frame(ID=names(sc@cpart), Cluster=as.factor(sc@cpart), sc@tsne)
tsne_plot_no_outline(FILL=df$Cluster, fill_colors = colors_fig)

# Extended data figure 3a - examination of cell types
signature_genes <- data.frame("monocytes"=c('CCR2', 'CLEC12A', 'PLAC8', 'FCN1', 'S100A9'),
                              "macrophages"=c("MRC1", "MS4A7", "CD163", "LYVE1", "STAB1"),
                              "microglia"= c('P2RY12', 'CX3CR1', 'CSF1R', 'TMEM119', 'SLC2A5'),
                              "tcell"=c('TRAC', 'TRBC2', 'CD52', 'IL32', NA),
                              "myeloid"=c('ITGAM',  'MS4A6A', 'TYROBP', 'CD14', NA),
                              "oligodendrocyte"=c('MBP',  'MOG', 'MAG', 'PLP1', NA),
                              "bcells"=c('CD79A', 'IGHG4', 'IGLL5', NA, NA),
                              "astrocyte"=c("GFAP", "HEPACAM","SOX9","AQP4",NA),
                              "apc"=c("CD74", "CD80", "CD86", "HLA-DRA", "CD40"), stringsAsFactors = F)

df <- data.frame(ID=names(sc@cpart), Cluster=sc@cpart, sc@tsne)

for (i in colnames(signature_genes)) {
  tryCatch({
    pl <- plotexptsne2(name2id(na.omit(signature_genes[[i]]), rownames(sc@ndata)),.sc=sc)
    print(pl)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  on.exit(dev.off())
}     

