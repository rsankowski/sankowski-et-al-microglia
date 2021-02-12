#open needed packages
library(readxl)
library(httr)

#load RaceID and custom functions
source(file.path("R", "RaceID3_StemID2_class.R"))
source(file.path("R", "sankowski-et-al-functions.R"))

#load counts
load(file.path("data", "prdata_gbm.RData"))

ids <- read_csv(file.path("data","gbm_microglia_clusters+embeddings_nn.csv"))[[2]]
prdata <- prdata[, ids]

## RaceID3
# initialize SCseq object with transcript counts
sc <- SCseq(prdata)
# filtering of expression data
sc <- filterdata(sc, mintotal=1500, 
                 minexpr=5, 
                 minnumber=1, 
                 maxexpr=Inf, 
                 downsample=FALSE, 
                 sfn=FALSE, 
                 hkn=FALSE, 
                 dsn=1, 
                 rseed=17000, 
                 CGenes=c("HSP90AA1__chr14", "HERPUD1__chr16", "HSPA1A__chr6", "MTRNR2L1__chrX", "MTRNR2L12__chr3", "MTRNR2L8__chr11", "FOS__chr14", "MALAT1__chr11", "JUN__chr1", "DUSP1__chr5"),  
                 FGenes=c("HSPB1__chr7", 
                          "HSPA6__chr1", 
                          "HSPH1__chr13", 
                          "HSPA1B__chr6", 
                          "HSP90AB1__chr6", 
                          "FOSB__chr19", 
                          "HSP90B1__chr12", 
                          "RPS16__chr19",
                          "DNAJB1__chr19", 
                          "H3F3B__chr17", 
                          "NEAT1__chr11", 
                          "RP11-212I21.4__chr16", 
                          "IVNS1ABP__chr1", 
                          "HSPA8__chr11", 
                          "HIST1H2BG__chr6", 
                          "HSPA5__chr9", 
                          "JUNB__chr19",
                          "ZFP36L1__chr14"), # "SAT1__chrX", "BTG2__chr1", "CDKN1A__chr6"
                 ccor=.4)

#save fdata with 3 counts minimum
write_csv(as.data.frame(rownames(sc@fdata)), "genes-with-minexpr>3-in-3-cells.csv")


# regress out the batch effect
# optional:
#build a table with cell Ids, cluster and conditions
data_t <- as.data.frame(t(sc@ndata))

#Add Cluster numbers
data_t$ID <- rownames(data_t)
data_t$cell_ID <- rownames(data_t)
data_t$Region_gm_wm <- ifelse(grepl('WM', data_t$ID), 'WM', 
                              ifelse(grepl('GM', data_t$ID), 'GM', 'Both')) 

data_t$ID <- gsub('_.*', '', data_t$ID)
data_t$ID <- gsub('-.*', '', data_t$ID)

data_t$ID <- gsub('GM|WM|all|micr|pos|17Pl1|17Pl2', '', data_t$ID)

#load batch info from the supplementary information on the paper website https://www.nature.com/articles/s41593-019-0532-y#Sec30
url1 <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-019-0532-y/MediaObjects/41593_2019_532_MOESM3_ESM.xlsx'
GET(url1, write_disk(tf <- tempfile(fileext = ".xlsx")))
batch_info <- read_excel(tf,sheet = 8, range = "B1:L9")
str(batch_info)

data_t <- left_join(data_t, batch_info)
table(data_t$ID, data_t$Batch)
vars <- as.data.frame(data_t$Batch)
sc@fdata <- varRegression(sc@fdata,vars)

# correct for cell cycle, proliferation, and expression of degradation markers by PCA
x <- CCcorrect(sc@fdata,
               vset=NULL,
               CGenes=c("HSP90AA1__chr14", "HSPA1A__chr6", "MTRNR2L1__chrX", "MTRNR2L12__chr3", "MTRNR2L8__chr11", "FOS__chr14", "MALAT1__chr11", "JUN__chr1", "DUSP1__chr5","HERPUD1__chr16"),
               ccor=.4,
               nComp=25,
               pvalue=.05,
               quant=.01,
               mode="pca")
# number of principal components that have been removed
x$n
# loadings of the first principal component that has been removed
y <- x$pca$rotation[,x$n[1]]
# genes from vset are either enriched in the head or the tail of this list
tail(y[order(y,decreasing=TRUE)],10)
# reassign the corrected expression matrix to sc@fdata
sc@fdata <- x$xcor

# k-medoids clustering
sc <- clustexp(sc,
               clustnr=30,
               bootnr=50,
               metric="pearson",
               do.gap=FALSE,
               sat=TRUE,
               SE.method="Tibs2001SEmax",
               SE.factor=.25,
               B.gap=50,
               cln=14, #new cln - it was based on the graphical output of the plotsaturation function
               rseed=17000,
               FUNcluster="kmedoids",
               FSelect=TRUE)
# compute t-SNE map
sc <- comptsne(sc,
               rseed=15555,
               sammonmap=FALSE,
               initial_cmd=TRUE,
               fast=TRUE,
               perplexity=30)
# detect outliers and redefine clusters
sc <- findoutliers(sc, 
                   outminc=5,
                   outlg=2,
                   probthr=1e-8,
                   thr=2**-(1:40),
                   outdistquant=.95)
# reassign clusters based on random forest
sc <- rfcorrect(sc,rfseed=12345,
                final=TRUE,
                nbfactor=5)

## diagnostic plots
# gap statistics: only if do.gap == TRUE
##plotgap(sc)
# plot within-cluster dispersion as a function of the cluster number: only if sat == TRUE
plotsaturation(sc,disp=TRUE)
# plot change of the within-cluster dispersion as a function of the cluster number: only if sat == TRUE
plotsaturation(sc)
# silhouette of k-medoids clusters
plotsilhouette(sc)
# Jaccard's similarity of k-medoids clusters
plotjaccard(sc)
# barchart of outlier probabilities
plotoutlierprobs(sc)
# regression of background model
plotbackground(sc)
# dependence of outlier number on probability threshold (probthr)
plotsensitivity(sc)
# heatmap of k-medoids cluster
clustheatmap(sc,final=FALSE,hmethod="single")
# heatmap of final cluster
clustheatmap(sc,final=TRUE,hmethod="single")
# highlight k-medoids clusters in t-SNE map
plottsne(sc,final=FALSE)
# highlight final clusters in t-SNE map
plottsne(sc,final=TRUE)
# highlight cell labels in t-SNE map
plotlabelstsne(sc,labels=sub("(\\_\\d+)","",names(sc@ndata)))
# highlight groups of cells by symbols in t-SNE map
plotsymbolstsne(sc,types=sub("(\\_\\d+)$","", names(sc@ndata)))

sc@fcol <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", rep('white', 3), "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",'#984EA3')
save(sc, file = file.path("data", "sc_gbm_myeloid_cells.RData"))

#please run this code if you would like to have the same t-SNE layout and clustering as the original paper
gbm_micr_nn <- read.csv(file.path("data", "gbm_microglia_clusters+embeddings_nn.csv"), row.names = 1)
sc@tsne <- gbm_micr_nn[, c("V1", "V2") ]

table(sc@cpart, gbm_micr_nn$Cluster)

clust_nn <- gbm_micr_nn$Cluster
names(clust_nn) <- rownames(gbm_micr_nn)

assertthat::assert_that(identical(names(sc@cpart), names(clust_nn)))
sc@cpart <- clust_nn
sc@tsne <- gbm_micr_nn[,c("V1", "V2")]

save(sc, file = file.path("data", "sc_gbm_microglia_nn.RData"))
