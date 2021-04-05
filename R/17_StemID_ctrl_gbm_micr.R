#load RaceID and custom functions
source(file.path("R", "RaceID3_StemID2_class.R"))
source(file.path("R", "sankowski-et-al-functions.R"))

#load sc data
load(file.path("data", "sc_gbm_microglia_nn.RData"))
load(file.path("data", "order_clusters_gbm_micr.RData"))
load(file.path("data", "retain_cl_gbm_micr.RData"))
order_clusters <- order_clusters[order_clusters %in% retain_cl]

if (!file.exists(file.path("data", "ltr_stemid_ctrl_gbm.RData"))) {
  ltr <- Ltree(sc)
  # computation of the entropy
  ltr <- compentropy(ltr)
  # computation of the projections for all cells
  ltr <- projcells(ltr,cthr=25,nmode=TRUE)
  # computation of the projections for all cells after randomization
  ltr <- projback(ltr,pdishuf=2000,nmode=TRUE,fast=TRUE,rseed=17000)
  # assembly of the lineage tree
  ltr <- lineagetree(ltr,pthr=0.05,nmode=TRUE,fast=TRUE)
  # determination of significant differentiation trajectories
  ltr <- comppvalue(ltr,pethr=0.05,nmode=TRUE,fast=TRUE)
  
  save(ltr, file = file.path("data", "ltr_stemid_ctrl_gbm.RData"))
} else {
  load(file.path("data", "ltr_stemid_ctrl_gbm.RData"))
}
  ## diagnostic plots
  # histogram of ratio between cell-to-cell distances in the embedded and the input space
  plotdistanceratio(ltr)
  # t-SNE map of the clusters with more than cthr cells including a minimum spanning tree for the cluster medoids
  plotmap(ltr)
  # visualization of the projections in t-SNE space overlayed with a minimum spanning tree connecting the cluster medoids
  plotmapprojections(ltr)
  # lineage tree showing the projections of all cells in t-SNE space
  plottree(ltr,showCells=TRUE,nmode=FALSE,scthr=.3)

# Supplementary Figure 5 
# lineage tree without showing the projections of all cells
plottree(ltr,showCells=FALSE,nmode=TRUE,scthr=.75)

# Figure 6a
#trajectory 11-13
n <- cellsfromtree(ltr,c(11,12,3,6,5,13))

if (!file.exists(file.path("data" , "s1d-StemID-trajectory-11,12,3,6,5,13.RData"))) {
  # filter out lowly expressed genes
  fs  <- filterset(ltr@sc@ndata,n=n$f,minexpr=4,minnumber=1)
  # compute self organizing map (SOM) of co-expressed genes
  s1d <- getsom(fs,nb=1000,k=5,locreg=TRUE,alpha=.5)
  save(s1d, file = file.path("data" , "s1d-StemID-trajectory-11,12,3,6,5,13.RData"))
} else {
    load(file.path("data" , "s1d-StemID-trajectory-11,12,3,6,5,13.RData"))
  }

ps  <- procsom(s1d,corthr=.85,minsom=3)

# coloring scheme for clusters ( vector with colours)
fcol <- colors_fig[retain_cl]
y    <- ltr@sc@cpart[n$f]
# plot average z-score for all modules derived from the SOM
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=TRUE)

#export gene modules for heatmap annotation and GO term analysis
df_ps <- data.frame('node' = ps$nodes, "GENEID" = names(ps$nodes)) 
write.csv(df_ps , '20181214-trajectory-cl11,12,3,6,5,13.csv')
