#load functions
source(file.path("R","RaceID2_class_ext.R"))

##### params #####
do.init    <- T
load.data  <- T

csamp <- c('GBM')

############################################ functions end ############################################

if ( load.data ){
  if (do.init){
    data   <- list()
    data.add <- list()
    cdiff <- list()
    sco   <- list()
    pca   <- list()
    entr  <- list()
    ltree <- list()
    gcl   <- list()
    lgres <- list()
    soPat2   <- list()
    marker <- list()
    net    <- list()
    regentr <- list()
    hyperv  <- list()
    traj    <- list()
  }
  
  gene2iso <- read.csv(file.path("data", "wgEncodeGencodeBasicV24_clean_genes2groups.tsv"),sep="\t",header=FALSE)
  ercc     <- read.csv(file.path("data", "ERCC_Controls_Analysis_length_ext.txt"),sep="\t",header=TRUE)
  for ( n in grep("nb.cell",names(ercc)) ){ercc[,n] <- 5/6 * ercc[,n]} # original CEL-seq
  
  n1 <- 1:97
  n2 <- c(1,98:193)
  n4 <- 1:193
  
  for ( s in csamp ){  
    if (s %in% c('GBM')){
      data.add[[s]] <- list()
      for ( m in c("t","b","c") ){
        z <- list()
        
        if ( s == "GBM" )  { nl <- c(Pat15_GM_1_1="H4_GM_1_1", # corresponding to Pat15
                                     Pat15_GM_2_2="H4_GM_2_2",
                                     Pat15_WM_1_3="H4_WM_1_3",
                                     Pat15_WM_2_4="H4_WM_2_4",
                                     Pat11_Pl1_1_5="H1_micr_1_5", #corresponding to Pat11
                                     Pat11_Pl1_2_6="H1_micr_2_6",
                                     Pat11_Pl2_1_7="H1_all_1_7",
                                     Pat11_Pl2_1_8="H1_all_1_8",
                                     Pat14_Plate1_1_1 = 'H3_1', #corresponding to Pat14
                                     Pat14_Plate1_2_2 = 'H3_2',
                                     Pat14_Plate2_1_3 = 'H3_3',
                                     Pat14_Plate2_2_4 = 'H3_4',
                                     Pat13_GM_1_1 = 'H2_GM_1', #corresponding to Pat13
                                     Pat13_GM_2_2 = 'H2_GM_2',
                                     GBM3_Pl1_1_5 = "GBM3_Pl1_1",
                                     GBM3_Pl1_2_6 = "GBM3_Pl1_2",
                                     GBM3_Pl2_1_7 = "GBM3_Pl2_1",
                                     GBM3_Pl2_2_8 = "GBM3_Pl2_2",
                                     GBM1_1_5 = "GBM1_1",
                                     GBM1_2_6 = "GBM1_2",
                                     GBM4_1_7 = "GBM4_1",
                                     GBM4_2_8 = "GBM4_2",
                                     GBM2_1_5 = "GBM2_1",
                                     GBM2_2_6 = "GBM2_2"
                                      
                                      
        );
        fl <- list()
        for ( i in names(nl) ) fl[[i]] <- n4;  di <- "" }
        
        
        for ( sl in names(nl) ){
          if ( length(di) > 0 ){
            x <- read.csv(file.path("data","counts",paste(sl,".cout",m,".csv",sep="")),sep="\t",header=TRUE)
          }else{
            x <- read.csv(file.path("data","counts",paste(sl,".cout",m,".csv",sep="")),sep="\t",header=TRUE)
          }
          x <- x[,fl[[sl]]]
          x <- merge(data.frame(GENEID=c(as.vector(gene2iso[,1]),as.vector(ercc[,1])),GROUP=c(as.vector(gene2iso[,2]),as.vector(ercc[,1]))),x,by="GENEID",all=TRUE)[,-1]
          names(x)[1] <- "GENEID"
          x[is.na(x[,2]),-1] <- 0
          x <- x[order(x$GENEID),]
          z[[sl]]  <- x
          names(z[[sl]])  <- c("GENEID",paste(nl[[sl]],sub("X","",names(z[[sl]])[-1]),sep="_"))
        }
        for ( i in 1:length(z) ) y <- if ( i == 1 ) z[[i]] else merge(y,z[[i]],by="GENEID")
        row.names(y) <- y$GENEID
        y <- y[,-1]
        y <- y[,apply(y,2,sum)>0]
        if ( m == "t" ){
          data[[s]] <- y
        }else{
          data.add[[s]][[m]] <- y
        }
      }
    }
  }
}


#Create prdata file
prdata <- data[[s]][grep("^(ERCC|MT-)",row.names(data[[s]]),invert=TRUE),] 
cs <- apply(prdata,2,sum)
prdata <- prdata[,cs > 500]
cs <- cs[cs > 500]
f <- t(prdata[name2id("KCNQ1OT1",rownames(prdata)),])/cs < .02

h <- rep(TRUE,nrow(prdata))
for ( g in c("KCNQ1OT1")){
  z <- apply(prdata,1,function(x,y) sigcor(x,y),y=t(prdata[name2id(g,rownames(prdata)),]))
  h <- h & ( is.na(z) | z < .65 )
}
prdata <- prdata[h,]

# save prdata file
save(prdata, file = file.path("data","prdata_gbm.RData"))

