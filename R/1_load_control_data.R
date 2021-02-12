#load functions
source(file.path("R","RaceID2_class_ext.R"))

##### params #####
do.init    <- T
load.data  <- T

csamp <- c('Ctrl')

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
    if (s %in% c('Ctrl')){
      data.add[[s]] <- list()
      for ( m in c("t","b","c") ){
        z <- list()
        
        if ( s == "Ctrl" )  { nl <- c(Pat12_1_1="Pat12_1_1",
                                         Pat12_2_2="Pat12_2_2",
                                         Pat12_1_5="Pat12_1_5",
                                         Pat12_2_6="Pat12_2_6",
                                         Pat15_GM_1_1="Pat15_GM_1_1",
                                         Pat15_GM_2_2="Pat15_GM_2_2",
                                         Pat15_WM_1_3="Pat15_WM_1_3",
                                         Pat15_WM_2_4="Pat15_WM_2_4",
                                         Pat11_Pl1_1_5="Pat11_micr_1_5", 
                                         Pat11_Pl1_2_6="Pat11_micr_2_6",
                                         Pat11_Pl2_1_7="Pat11_all_1_7",
                                         Pat11_Pl2_1_8="Pat11_all_1_8",
                                         Pat10_GM_1_5="Pat10_GM_1", 
                                         Pat10_GM_2_6="Pat10_GM_2", 
                                         Pat10_WM_2_8= 'Pat10_WM_1', 
                                         Pat7_GM_2_2 = 'Pat7_GM_2',
                                         Pat7_WM2_MS_Ctrl_1_7 = 'Pat7_WM_MSCtrl_1',
                                         Pat7_WM2_MS_Ctrl_2_8 = 'Pat7_WM_MSCtrl_2',
                                         Pat7_WM_1_3 = 'Pat7_WM_1',
                                         Pat7_WM_2_4 = 'Pat7_WM_2',
                                         Pat14_Plate1_1_1 = 'Pat14_1',
                                         Pat14_Plate1_2_2 = 'Pat14_2',
                                         Pat14_Plate2_1_3 = 'Pat14_3',
                                         Pat14_Plate2_2_4 = 'Pat14_4',
                                         Pat16_Plate1_1_5 = 'Pat16_1', #this sample was originally loaded with the data in the paper, but not included in the final manuscript. I'm letting it run here to ensure best possible reproducibility in the correlations below
                                         Pat16_Plate1_2_6 = 'Pat16_2', #this sample was originally loaded with the data in the paper, but not included in the final manuscript. I'm letting it run here to ensure best possible reproducibility in the correlations below
                                         Pat16_Plate2_1_7 = 'Pat16_3', #this sample was originally loaded with the data in the paper, but not included in the final manuscript. I'm letting it run here to ensure best possible reproducibility in the correlations below
                                         Pat16_Plate2_2_8 = 'Pat16_4', #this sample was originally loaded with the data in the paper, but not included in the final manuscript. I'm letting it run here to ensure best possible reproducibility in the correlations below
                                         Pat4_GM_1_5 = 'Pat4_GM_1',
                                         Pat4_GM_2_6 = 'Pat4_GM_2',
                                         Pat4_WM_1_11 = 'Pat4_WM_1',
                                         Pat4_WM_1_12 = 'Pat4_WM_2',
                                         Pat6_Plate1_1_1 = 'Pat6_Plate1_1',
                                         Pat6_Plate1_2_2 = 'Pat6_Plate1_2',
                                         Pat6_Plate2_1_6 = 'Pat6_Plate2_1',
                                         Pat6_Plate2_1_6 = 'Pat6_Plate2_1',
                                         Pat2_GM_1_1 = 'Pat2_GM_1',
                                         Pat2_GM_2_2 = 'Pat2_GM_2',
                                         Pat2_WM_1_3 = 'Pat2_WM_1',
                                         Pat2_WM_2_4 = 'Pat2_WM_2',
                                         Pat3_GM_1_1 = 'Pat3_GM_1',
                                         Pat3_GM_2_2 = 'Pat3_GM_2',
                                         Pat3_WM_1_3 = 'Pat3_WM_1',
                                         Pat3_WM_2_4 = 'Pat3_WM_2',
                                         Pat9_GM_1_1 = 'Pat9_GM_1',
                                         Pat9_GM_2_2 = 'Pat9_GM_2',
                                         Pat9_WM_1_5 = 'Pat9_WM_1',
                                         Pat9_WM_2_6 = 'Pat9_WM_2',
                                         Pat8_GM_1_1 = 'Pat8_GM_1',
                                         Pat8_GM_2_2 = 'Pat8_GM_2',
                                         Pat8_WM_1_3 = 'Pat8_WM_1',
                                         Pat8_WM_2_4 = 'Pat8_WM_2_4',
                                         Pat1_GM_1_5 = 'Pat1_GM_1',
                                         Pat1_GM_2_6 = 'Pat1_GM_2',
                                         Pat1_WM_1_7 = 'Pat1_WM_1',
                                         Pat1_WM_2_8 = 'Pat1_WM_2',
                                         Pat5_GM_1_1 = 'Pat5_GM_1',
                                         Pat5_GM_2_2 = 'Pat5_GM_2',
                                         Pat5_WM_1_3 = 'Pat5_WM_1',
                                         Pat5_WM_2_4 = 'Pat5_WM_2',
                                         Pat13_GM_1_1 = 'Pat13_GM_1',
                                         Pat13_GM_2_2 = 'Pat13_GM_2'
                                         
                                         
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
save(prdata, file = file.path("data","prdata.RData"))

