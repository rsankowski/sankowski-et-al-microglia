#Sankowski et al -functions and plots
require(tidyverse)
require(viridis)
source("functions/RaceID3_StemID2_class.R")

#colors
colors_many <- toupper(read_csv('data/20180313_trubetskoy_colors.csv')[[3]])
colors <- toupper(c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999'))
colors_pat <- toupper(c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5'))
colors_fig <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",'#984EA3', "light grey", "grey", "dark grey", "#696969")


#run RaceID
run_raceid <- function(.prdata = prdata,
                       .mintotal = 1500, 
                       .CGenes = c("HSP90AA1__chr14", "HSPA1A__chr6", "MTRNR2L1__chrX", "MTRNR2L12__chr3", "MTRNR2L8__chr11", "FOS__chr14", "MALAT1__chr11", "JUN__chr1", "DUSP1__chr5"),
                       .FGenes = c("HSPB1__chr7", "HSPA6__chr1", "HSPH1__chr13", "HSPA1B__chr6", "HSP90AB1__chr6", "FOSB__chr19", "HSP90B1__chr12", "RPS16__chr19","DNAJB1__chr19", "H3F3B__chr17", "HERPUD1__chr16", "NEAT1__chr11", "RP11-212I21.4__chr16", "IVNS1ABP__chr1", "HSPA8__chr11", "HIST1H2BG__chr6", "HSPA5__chr9", "JUNB__chr19","ZFP36L1__chr14"),
                       #location_batch_info='data/batch_info.csv',
                       PCs_for_clustering = NULL,
                       .cln=cln
) {
  require(tidyverse)
  
  
  data <- SCseq(.prdata)
  # filtering of expression data
  data <- filterdata(data, mintotal=.mintotal, 
                     minexpr=5, 
                     minnumber=1, 
                     maxexpr=Inf, 
                     downsample=FALSE, 
                     sfn=FALSE, 
                     hkn=FALSE, 
                     dsn=1, 
                     rseed=17000, 
                     CGenes=.CGenes,  
                     FGenes=.FGenes,
                     ccor=.4)
  
  # regress out the batch effect
  # optional:
  #build a table with cell Ids, cluster and conditions
  data_t <- as.data.frame(t(data@fdata))
  
  data_t$ID <- rownames(data_t)
  data_t$cell_ID <- rownames(data_t)
  data_t$Region_gm_wm <- ifelse(grepl('WM', data_t$ID), 'WM', 
                                ifelse(grepl('GM', data_t$ID), 'GM', 'Both')) 
  
  data_t$ID <- gsub('_.*', '', data_t$ID)
  data_t$ID <- gsub('-.*', '', data_t$ID)
  
  data_t$ID <- gsub('GM|WM|all|micr|pos|17Pl1|17Pl2', '', data_t$ID)
  
  #regress out batch effects 
  if (F) {
    batch_info <- read_csv(location_batch_info)[,-1]
    data_t <- merge(data_t, batch_info)
    table(data_t$ID, data_t$Batch)
    vars <- as.data.frame(data_t$Batch[data_t$cell_ID %in% colnames(data@fdata)])
    data@fdata <- varRegression(data@fdata,vars)}
  
  # correct for cell cycle, proliferation, and expression of degradation markers by PCA
  # optional:
  x <- CCcorrect(data@fdata,
                 vset=NULL,
                 CGenes=.CGenes,
                 ccor=.4,
                 nComp=,
                 pvalue=.05,
                 quant=.01,
                 mode="pca")
  # number of principal components that have been removed
  x$n
  # loadings of the first principal component that has been removed
  y <- x$pca$rotation[,x$n[1]]
  # genes from vset are either enriched in the head or the tail of this list
  tail(y[order(y,decreasing=TRUE)],10)
  # reassign the corrected expression matrix to data@fdata
  data@fdata <- x$xcor
  
  # k-medoids clustering
  data <- clustexp(data,
                   clustnr=30,
                   bootnr=50,
                   metric="pearson",
                   do.gap=FALSE,
                   sat=TRUE,
                   SE.method="Tibs2001SEmax",
                   SE.factor=.25,
                   B.gap=50,
                   cln=.cln,
                   rseed=17000,
                   FUNcluster="kmedoids",
                   FSelect=TRUE)
  # compute t-SNE map
  data <- comptsne(data,
                   rseed=15555,
                   sammonmap=FALSE,
                   initial_cmd=TRUE,
                   fast=TRUE,
                   perplexity=30)
  # detect outliers and redefine clusters
  data <- findoutliers(data, 
                       outminc=5,
                       outlg=2,
                       probthr=1e-8,
                       thr=2**-(1:40),
                       outdistquant=.95)
  # reassign clusters based on random forest
  data <- rfcorrect(data,rfseed=12345,
                    final=TRUE,
                    nbfactor=5)
  
  data
}

#differential gene expression
differential_test <- function(data, pval.cutoff = 0.05, homeDirectory = home, subDirectory = "Cluster specific genes") {
  mainDir <- homeDirectory
  subDir <- subDirectory
  subsubDir_up <- 'Up'
  subsubDir_down <- 'Down'
  
  if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
  } else {
    dir.create(file.path(mainDir, subDir))
    dir.create(file.path(mainDir, subDir, subsubDir_up))
    dir.create(file.path(mainDir, subDir, subsubDir_down))
    setwd(file.path(mainDir, subDir))
    
  }
  for (i in c(1: max(data@cpart))) {
    cl <- names(data@cpart[data@cpart %in% c(i)])
    rest <- names(data@cpart[data@cpart %in% data@cpart[data@cpart != i]])
    diffexp <- diffexpnb(data@ndata,rest,cl,norm=F,vfit=data@background$vfit) 
    diffexpgenes <- diffexp[["res"]]
    diffexpgenes <- subset(diffexpgenes, diffexpgenes$pval < pval.cutoff)
    diffexpgenes <- subset(diffexpgenes, abs(diffexpgenes$log2FoldChange) > 0)
    diffexpgenes$GENEID <- rownames(diffexpgenes)
    diffexpgenes <- diffexpgenes %>% dplyr::arrange(padj)
    #diffexpgenes$GENEID <- gsub('__.*', '', diffexpgenes$GENEID)
    diffexpgenes$Cluster <- rep(i, nrow(diffexpgenes))
    diffgene_up <- diffexpgenes[diffexpgenes$log2FoldChange > 0, ]
    diffgene_down <- diffexpgenes[diffexpgenes$log2FoldChange < 0, ]
    write.csv(diffgene_up, file = paste0('Up/', "diffgenes_cl", as.character(i), "_up_rest.csv"))
    write.csv(diffgene_down, file = paste0('Down/', "diffgenes_cl", as.character(i), "_down_rest.csv"))
    
    
    tryCatch({
      png(paste0('Up/', date,'-MA-plot-Cl', as.character(i) ,'_rest.png'), res = 300, width =7, height = 7, units = 'in')
      plotdiffgenesnb(diffexp,show_names = T, pthr=.01, lthr = .0001, mthr = .0001, padj=T)
      dev.off()
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    for (j in c(1: max(data@cpart))) {
      cl <- names(data@cpart[data@cpart %in% i])
      rest <- names(data@cpart[data@cpart %in% j])
      diffexp <- diffexpnb(data@ndata,rest,cl,norm=F,vfit=data@background$vfit) 
      diffexpgenes <- diffexp[["res"]]
      diffexpgenes <- subset(diffexpgenes, diffexpgenes$pval < pval.cutoff)
      diffexpgenes <- subset(diffexpgenes, abs(diffexpgenes$log2FoldChange) > 0)
      diffexpgenes$GENEID <- rownames(diffexpgenes)
      diffexpgenes <- diffexpgenes %>% dplyr::arrange(padj)
      #diffexpgenes$GENEID <- gsub('__.*', '', diffexpgenes$GENEID)
      diffexpgenes$Cluster <- rep(i, nrow(diffexpgenes))
      diffgene_up <- diffexpgenes[diffexpgenes$log2FoldChange > 0, ]
      diffgene_down <- diffexpgenes[diffexpgenes$log2FoldChange < 0, ]
      write.csv(diffgene_up, file = paste0('Up/', "diffgenes_cl", as.character(i), "_up_vs_cl", as.character(j), ".csv"))
      write.csv(diffgene_down, file = paste0('Down/', "diffgenes_cl", as.character(i), "_down_vs_cl", as.character(j), ".csv"))
      
      tryCatch({
        png(paste0('Up/', date,'-MA-plot-Cl', as.character(i) , "_vs_cl", as.character(j), ".png"), res = 300, width =7, height = 7, units = 'in')
        plotdiffgenesnb(diffexp,show_names = T, pthr=.01, lthr = .0001, mthr = .0001, padj=T)
        dev.off()
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      
    }
    
  }
  
}

#build a dataset
build_tidy_dataset <- function(data = sc, percent_cluster = 1, location_batch_info='/home/roman/Documents/Single cell analysis/20170921 Healthy microglia/20171210-batch-information.csv') {
  
  cell_numbers <-as.numeric()
  for (i in 1:length(unique(data@cpart)))
  {
    cell_numbers[i] <- length(data@cpart[data@cpart==i])
  }
  names(cell_numbers) <- c(1:length(unique(data@cpart)))
  retain_cl <- as.numeric(names(cell_numbers[cell_numbers > dim(data@ndata)[2]/(100/percent_cluster)]))
  order_clusters <- clustheatmap(data,final=TRUE)
  order_clusters <- order_clusters[order_clusters %in% retain_cl]
  df <- data.frame('Cluster' = data@cpart, data@tsne) 
  df$ID <- rownames(df)
  df$cell_ID <- rownames(df)
  df$Condition <- gsub('_.*', '', rownames(df))
  df <- df[df$Cluster %in% retain_cl,] 
  df$Cluster <- factor(df$Cluster, levels = order_clusters)
  df$Cluster <- droplevels(df$Cluster)
  
  #add batch info
  batch_info <- read_csv(location_batch_info)[,-1]
  colnames(batch_info)[c(1,4)] <- c('Condition', 'Diagnosis')
  df <- left_join(df, batch_info)
  df$Condition <- as.factor(df$Condition)
  df$Condition <- reorder(df$Condition, df$Age)
  df$anon_ID <- df$Condition
  levels(df$anon_ID) <- paste0("Pat", 1:length(unique(df$Condition)))
  df$Region <- ifelse(grepl('WM',df$ID), 'WM', 
                      ifelse(grepl('GM',df$ID), 'GM', 'Mixed'))
  df$Region <- as.factor(df$Region)
  
  df
}

#enrichment test
hyper_test <- function(data1 = sc, data2 = df, var1 = "Cluster", var2 = "Region") {
  require(tidyverse)            
  clusters <- as_tibble(table(data1@cpart))
  colnames(clusters) <- c(var1, 'cluster_size')
  vars <- as_tibble(table(data2[,var1], data2[,var2]))
  colnames(vars) <- c(var1, var2, "freq_var2")
  vars_wide <- spread(vars, var2, freq_var2)
  
  vars_df <- vars_wide %>%
    left_join(clusters)
  
  
  #hypergeometric test
  #Male
  test_df<- data.frame(q=vars_df[,3], 
                       m=sum(vars_df[,3]), 
                       n=sum(vars_df[,2]),
                       k=vars_df[,4])
  
  p_hyper <- apply(test_df, MARGIN = 1, function(x) 1-phyper(x[[1]]-1, x[[2]], x[[3]], x[[4]])) #probability to get q or more successes in populaton
  padj <- p.adjust(p_hyper, method="BH")
  
  
  test_df$p_hyper <- p_hyper
  test_df$padj <- padj
  test_df$Cluster <- vars_df$Cluster
  test_df$Significance <- ifelse(test_df$padj<0.05, '*',
                                 ifelse(test_df$padj<0.01, '**',
                                        ifelse(test_df$padj<0.001, '***','n.s.')))
  test_df$enrichment_var <- colnames(test_df)[1]
  
  #Female
  test_df2 <- data.frame(q=vars_df[,2], 
                         m=sum(vars_df[,2]), 
                         n=sum(vars_df[,3]),
                         k=vars_df[,4])
  p_hyper <- apply(test_df2, MARGIN = 1, function(x) 1-phyper(x[[1]]-1, x[[2]], x[[3]], x[[4]])) #probability to get q or more successes in populaton
  padj <- p.adjust(p_hyper, method="BH")
  
  test_df2$p_hyper <- p_hyper
  test_df2$padj <- padj
  test_df2$Cluster <- vars_df$Cluster
  test_df2$Significance <- ifelse(test_df2$padj<0.05, '*',
                                  ifelse(test_df2$padj<0.01, '**',
                                         ifelse(test_df2$padj<0.001, '***','n.s.')))
  
  test_df2$enrichment_var <- colnames(test_df2)[1]
  
  return(bind_rows(test_df, test_df2))
}

#plots
plotexptsne2 <-   function(gene, .df=df, .sc=counts, point_size=3, logsc=FALSE, line_width=0.25) {
  gene = name2id(gene, rownames(.sc))
  l <- colSums(.sc[rownames(.sc)[rownames(.sc) %in% gene],])
  mi <- min(l)
  ma <- max(l)
  ColorRamp <- colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100)
  ColorLevels <- seq(mi, ma, length = length(ColorRamp))
  v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
  
  kk <- bind_cols(data.frame('l'=l), .df) %>% arrange(l)
  
  if(logsc) {
    plot <- ggplot(kk, aes(V1, V2, fill = log(l))) +
      geom_point(size = point_size, pch = 21, stroke=line_width) +
      scale_fill_gradientn('', colors = ColorRamp) +
      theme_void() +
      labs(title = paste(gene, collapse = ',')) 
    return(plot)
  }
  else {
    plot <- ggplot(kk, aes(V1, V2, fill = l)) +
      geom_point(size = point_size, pch = 21, stroke=line_width) +
      scale_fill_gradientn('', colors = ColorRamp) +
      theme_void() +
      labs(title = paste(gene, collapse = ','))
    return(plot)
  }
  
}

#

cluster_stacked_barplot <- function(data=df,FILL = Condition, fill_colors = colors_pat, fill_variable_name = "Condition") {
  #Stacked cluster plot
  cluster_stack_plot <- ggplot(data, aes(Cluster, fill = FILL, group = FILL)) + 
    geom_bar(position = 'fill', color = 'black', lwd = 0.25) +
    #geom_line(size = 2) +
    theme_minimal() +
    #scale_x_continuous(limits = c(0.5,19.5), breaks = c(1:19)) +
    theme(panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          text=element_text(size=17)) + 
    scale_fill_manual(fill_variable_name, values = fill_colors) +
    labs(title = 'Conditions in Cluster', y='#cells in Cluster / total #cells', x='Cluster')
  
  cluster_stack_plot
}


normalized_stack_plot <- function(data=df, Var1 = 'Cluster', Var2 = "Diagnosis", fill_colors = toupper(c('#f1a340','#998ec3')), line_width=0.25) {
  #normalized stacked map
  pats <- table(data[,Var1], data[,Var2])
  pats_prop <- as.data.frame(prop.table(pats, 2))
  
  ggplot(pats_prop, aes(Var1, Freq, fill = Var2)) +
    geom_bar(stat = 'identity' , position = "fill", color = 'black', lwd=line_width) +
    theme_minimal() +
    theme(panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          text=element_text(size=17)) +
    scale_fill_manual(Var2, values = fill_colors) +
    labs(title = paste0(Var2,' in ', Var1), y='normalized % of cells', x=Var1)
}


#tsne plots

tsne_plot <- function(data = df, FILL = Condition, fill_colors = colors_pat, point_outline = "black", point_size = 3, line_width = 0.25, point_shape=21, alpha_value = 1) {
  tsne_plot <- ggplot(data, aes(V1, V2, fill = FILL, color = FILL)) +
    geom_point(pch = point_shape, size = point_size, stroke = line_width, color = point_outline, alpha = alpha_value) +
    theme(panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          text=element_text(size=17),
          axis.text.x=element_blank(), 
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.title = element_blank(),
          legend.key = element_blank()) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = fill_colors, guide = FALSE) 
  
  tsne_plot 
}

#tsne plot without outline
tsne_plot_no_outline <- function(data = df, FILL = Condition, fill_colors = colors_pat, point_outline = "black", point_size = 3, line_width = 0.25, point_shape=21, alpha_value = 1) {
  tsne_plot <- ggplot(data, aes(V1, V2, fill = FILL, color = FILL)) +
    geom_point(pch = point_shape, size = point_size, stroke = line_width, alpha = alpha_value) +
    theme(panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          text=element_text(size=17),
          axis.text.x=element_blank(), 
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.title = element_blank(),
          legend.key = element_blank()) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = fill_colors, guide = FALSE) 
  
  tsne_plot 
}

#Marimekko plot with stats
mosaicGG <- function(data, X, FILL, rect_col = 'white', line_width = 0.25) {
  require(dplyr)
  require(reshape2)
  #require(ggthemes)
  # Proportions in raw data
  DF <- as.data.frame.matrix(table(data[[X]], data[[FILL]]))
  DF$groupSum <- rowSums(DF)
  DF$xmax <- cumsum(DF$groupSum)
  DF$xmin <- DF$xmax - DF$groupSum
  DF$X <- row.names(DF)
  DF$groupSum <- NULL
  DF_melted <- melt(DF, id = c("X", "xmin", "xmax"), variable.name = "FILL")
  DF_melted <- DF_melted %>%
    group_by(X) %>%
    mutate(ymax = cumsum(value/sum(value)),
           ymin = ymax - value/sum(value))
  
  # Chi-sq test
  results <- chisq.test(table(data[[FILL]], data[[X]])) # fill and then x
  resid <- melt(results$residuals)
  names(resid) <- c("FILL", "X", "residual")
  
  # Merge data
  DF_all <- merge(DF_melted, resid)
  
  # Positions for labels
  DF_all$xposn <- DF_all$xmin + (DF_all$xmax - DF_all$xmin)/2
  index <- DF_all$xmax == max(DF_all$xmax)
  DF_all$yposn <- DF_all$ymin[index] + (DF_all$ymax[index] - DF_all$ymin[index])/2
  
  # Plot
  g <- ggplot(DF_all, aes(ymin = ymin,  ymax = ymax, xmin = xmin,
                          xmax = xmax, fill = residual)) +
    geom_rect(col = rect_col, lwd = line_width) +
    geom_text(aes(x = xposn, label = X),
              y = 1, size = 3, angle = 90, hjust = 1, show.legend = FALSE, check_overlap = TRUE) +
    geom_text(aes(x = max(xmax),  y = yposn, label = FILL),
              size = 3, hjust = 1, show.legend = FALSE, check_overlap = TRUE) +
    scale_fill_gradient2("Residuals") +
    scale_x_continuous(X, expand = c(0,0)) +
    scale_y_continuous("Proportion", expand = c(0,0)) +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(g)
}

#Marimekko plot without stats
mosaicGG2 <- function(data, X, FILL, colors = c(colors_pat,colors_many), rect_col = 'white', line_width = 0.25) {
  require(dplyr)
  require(reshape2)
  #require(ggthemes)
  # Proportions in raw data
  DF <- as.data.frame.matrix(table(data[[X]], data[[FILL]]))
  DF$groupSum <- rowSums(DF)
  DF$xmax <- cumsum(DF$groupSum)
  DF$xmin <- DF$xmax - DF$groupSum
  DF$X <- row.names(DF)
  DF$groupSum <- NULL
  DF_melted <- melt(DF, id = c("X", "xmin", "xmax"), variable.name = "FILL")
  DF_melted <- DF_melted %>%
    group_by(X) %>%
    mutate(ymax = cumsum(value/sum(value)),
           ymin = ymax - value/sum(value))
  
  # Chi-sq test
  results <- suppressWarnings(chisq.test(table(data[[FILL]], data[[X]]))) # fill and then x
  resid <- melt(results$residuals)
  names(resid) <- c("FILL", "X", "residual")
  
  # Merge data
  DF_all <- merge(DF_melted, resid)
  #DF_all  <- DF_melted
  # Positions for labels
  DF_all$xposn <- DF_all$xmin + (DF_all$xmax - DF_all$xmin)/2
  index <- DF_all$xmax == max(DF_all$xmax)
  DF_all$yposn <- DF_all$ymin[index] + (DF_all$ymax[index] - DF_all$ymin[index])/2
  
  # Plot
  g <- ggplot(DF_all, aes(ymin = ymin,  ymax = ymax, xmin = xmin,
                          xmax = xmax, fill = FILL)) +
    geom_rect(col = rect_col, lwd = line_width) +
    geom_text(aes(x = xposn, label = X),
              y = 1, size = 3, angle = 90, hjust = 1, show.legend = FALSE, check_overlap = TRUE) +
    geom_text(aes(x = max(xmax),  y = yposn, label = FILL),
              size = 3, hjust = 1, show.legend = FALSE, check_overlap = TRUE) +
    scale_fill_manual(FILL, values = colors) +
    scale_x_continuous(X, expand = c(0,0)) +
    scale_y_continuous("Proportion", expand = c(0,0)) +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(g)
}


#GO terms enrichment analysis

enrichment_analysis <- function(data1 = sc, 
                                data2 = df_upgenes, 
                                ontology = "BP",
                                stress_genes = c("HSPA1A","MTRNR2L8","MTRNR2L12","HSP90AA1","MALAT1","ZFP36L1","ZFP36","MTRNR2L1","FOS","MALAT1","HSPB*","DUSP1","HSPH1","HSPA*","JUN","HSP90B1","RPS16","DNAJB1","H3F3B","HERPUD1","NEAT1","IVNS1ABP","HIST1H2BG","RP*","XIST")) 
{
  #source("https://bioconductor.org/biocLite.R")
  #biocLite('org.Hs.eg.db')
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(tidyverse)
  keytypes(org.Hs.eg.db)
  require(pheatmap)
  
  data2 <- data2[!grepl(paste0("^(", paste(stress_genes, collapse='|'), ")"), data2$GENEID),  ]
  
  #load background genes
  back_genes <- rownames(data1@ndata)[which(apply(data1@ndata > .1, 1, sum)>0)]
  background <- bitr(sub('_.*', '',back_genes), fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                     OrgDb = org.Hs.eg.db)
  
  #define empty data frame to collect data
  enrich_up <- data.frame(matrix(ncol = 10))
  colnames(enrich_up) <- c('ID','Description', 'GeneRatio', 'BgRatio' ,'pvalue', 'p.adjust', 'qvalue', 'geneID','Count' , 'Cluster')
  
  for (i in seq_along(unique(data2$Cluster)))  {
    
    tryCatch({
      gene <- gsub('_.*', '',data2$GENEID[data2$Cluster == i])
      gene.data2 <- bitr(gene, fromType = "SYMBOL",
                         toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                         OrgDb = org.Hs.eg.db)
      
      ggo <- groupGO(gene     = gene.data2[,3],
                     OrgDb    = org.Hs.eg.db,
                     ont      = ontology,
                     level    = 3,
                     readable = TRUE)
      
      
      ego <- enrichGO(gene          = gene.data2[,3],
                      universe      = background[,3],
                      OrgDb         = org.Hs.eg.db,
                      minGSSize     = 1,
                      ont           = ontology,
                      pool          = TRUE,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
      
      ego_simpl <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = "Wang")
      ego_simpl2 <- ego_simpl[!duplicated(ego_simpl@result$geneID)]
      ego_simpl2$Cluster <- rep(as.character(i, nrow(ego_simpl2)))
      enrich_up <- rbind(enrich_up, ego_simpl2)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  na.omit(enrich_up)
  
  
}

run_gsea <- function(data = df_upgenes, 
                     seed = 79106,
                     term2gene = "GO",
                     stress_genes = c("HSPA1A","MTRNR2L8","MTRNR2L12","HSP90AA1","MALAT1","ZFP36L1","ZFP36","MTRNR2L1","FOS","MALAT1","HSPB*","DUSP1","HSPH1","HSPA*","JUN","HSP90B1","RPS16","DNAJB1","H3F3B","HERPUD1","NEAT1","IVNS1ABP","HIST1H2BG","RP*","XIST")) {
  TERM2GENE <- read.gmt("/home/roman/Documents/Single cell analysis/msigdb.v6.2.symbols.gmt")
  
  if (term2gene == "GO") {
    TERM2GENE <- read.gmt("/home/roman/Documents/Single cell analysis/msigdb.v6.2.symbols.gmt")}
  if (term2gene == "immuno") {
    TERM2GENE <- read.gmt("/home/roman/Documents/Single cell analysis/c7.all.v6.2.symbols.gmt")
  }
  data <- data[!grepl(paste0("^(", paste(stress_genes, collapse='|'), ")"), data$GENEID),  ]
  df_gsea <- list()
  for (i in seq_along(unique(data$Cluster))) {
    tryCatch({
      data1 <- data %>%
        filter(Cluster == i) %>%
        group_by(GENEID) %>%
        filter(foldChange == max(foldChange))
      
      colnames(gene_symbols)[1] <- "GENEID"
      data1 <- left_join(data1, gene_symbols)
      
      gene_names <- data1$foldChange
      names(gene_names) <- data1$GENEID
      geneList <- sort(gene_names, decreasing = T)
      
      df_gsea1 <- GSEA(geneList=geneList, TERM2GENE = TERM2GENE)
      df_gsea[[i]] <- df_gsea1 
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  df_gsea
}

#GO dot plot
go_dot_plot <- function(data = enrich_up, var1 = "Cluster", var2 = "Description", point_size = "GeneCount", FILL = "qvalue", line_width = 0.25, fill_name = "q-value",...) {
  require(tidyverse)
  require(viridis)
  dot_plot <- ggplot(enrich_up, aes_string(x=var1, y=var2, size = point_size, fill= -log10(data[[FILL]]))) + 
    geom_point(pch=21, stroke=line_width) +
    scale_fill_viridis(fill_name) +
    theme_light() +
    theme(text=element_text(size=10),
          axis.title.y=element_blank())
  print(dot_plot)
}

#Histo cell counting
#count images in the target directory
count_images <- function(source_path) {
  
  #counting existing images in dataset
  #count all the available images
  
  distance = data.frame()
  distance_summary_images = data.frame()
  nndistance = data.frame()
  delaunay = data.frame()
  derichelet = data.frame()
  all_coordinates = data.frame()
  delaunay_dist = data.frame()
  
  
  for (i in list.files(source_path)) {
    path=list.files(paste0(source_path, '/', as.character(i)))
    for (file in path) {
      for (tables in list.files(paste0(source_path, '/', as.character(i), '/', as.character(file)))) {
        table_name = paste0(source_path, '/', as.character(i), '/', as.character(file), '/', as.character(tables))
        
        #generate table with summary stats
        distance_summary2 = data.frame(ID = i, 
                                       Condition = file, 
                                       Image = tables
        )
        distance_summary_images = rbind(distance_summary_images, distance_summary2)
        
        
      }
      
      
    }
    
    
    
    
    
    
  }
  
  distance_summary_images$Image <- gsub('image', '', distance_summary_images$Image)
  distance_summary_images$Image <- gsub('.tif', '.csv', distance_summary_images$Image)
  
  return(distance_summary_images)
}

#pattern analysis
pattern_analysis <- function(source_path, spatstat = FALSE) {
  #First version :20181015
  #last update :20181031
  #only usable for 20x magnification images
  
  #list.files()
  distance = data.frame()
  distance_summary = data.frame()
  delaunay = data.frame()
  derichelet = data.frame()
  all_coordinates = data.frame()
  delaunay_dist = data.frame()
  
  
  for (i in list.files(source_path)) {
    path=list.files(paste0(source_path, '/', as.character(i)))
    for (file in path) {
      tryCatch({
        for (tables in list.files(paste0(source_path, '/', as.character(i), '/', as.character(file)))) {
          if (grepl(".csv$", tables)) {table_name = paste0(source_path, '/', as.character(i), '/', as.character(file), '/', as.character(tables))
          b = read.csv(table_name, 1, sep = ',', colClasses='numeric')[,6:7]
          #calibrate the imagej output to microns
          b$X = ((2560/8.51)*0.17) * b$X
          b$Y = ((1920/6.39)*0.17) * b$Y
          
          #Generate data frame with all coordinates
          coord = data.frame(ID = rep(i, nrow(b)), 
                             Region = rep(file, nrow(b)), 
                             Image = rep(tables, nrow(b)), 
                             X = b$X,
                             Y = b$Y)
          all_coordinates = rbind(all_coordinates, coord)
          
          #Compute distances
          distances = as.vector(dist(as.matrix(b)))
          distances2 = data.frame(ID = rep(i, length(distances)), 
                                  Region = rep(file, length(distances)), 
                                  Image = rep(tables, length(distances)), 
                                  Distance = distances)
          distance = rbind(distance, distances2) 
          
          #generate table with summary stats
          distance_summary2 = data.frame(ID = i, 
                                         Condition = file, 
                                         Cell_count = nrow(b), 
                                         Cell_density = nrow(b)/((2560*1920*(0.17)^2/1000000)), 
                                         Mean_cell_dist = mean(distances), 
                                         Median_cell_dist = median(distances), 
                                         Range_cell_dist = max(distances)-min(distances),
                                         Image = tables,
                                         Antigen = strsplit(source_path, "/")[[1]][6])
          distance_summary = rbind(distance_summary, distance_summary2)
          
          if (spatstat) { 
            require(spatstat)
            #Spatstat package anaysis
            if (nrow(b) > 2) {
              #Voronoi tessallation
              derich = ppp(b$X,b$Y, xrange = range(0,(2560 * 0.17)), yrange = range(0,(1920 * 0.17)))
              derich2  = dirichletAreas(derich)
              #build data frame
              derich3= data.frame(ID = rep(i, length(derich2)), 
                                  Region = rep(file, length(derich2)), 
                                  Image = rep(tables, length(derich2)), 
                                  Derichelet_tiles = derich2)
              derichelet = rbind(derichelet, derich3) 
              
              if (nrow(b) > 3) {
                #delaunay
                delaun = delaunay(derich)
                delaun2 = tiles(delaun)
                delaun3 = unlist(lapply(delaun2, area.owin))
                
                #Build data frame
                delaun4 = data.frame(ID = rep(i, length(delaun3)), 
                                     Region = rep(file, length(delaun3)), 
                                     Image = rep(tables, length(delaun3)), 
                                     Delaunay_tiles = delaun3)
                delaunay = rbind(delaunay, delaun4)
                
                #delaunay distance
                delaun10 = data.frame(ID = i,
                                      Region = file, 
                                      Image = tables,
                                      Modus = 'Patient', 
                                      Del_distance = mean(delaunayDistance(derich)))
                delaunay_dist = rbind(delaunay_dist, delaun10)
              }
              
              
            }
            
          }
          
          }
          
        }
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
    }
    
  }
  
  return(distance_summary)
  
}


#make long dataset 
make_data_long <- function(data1 = df, data2 = sc, order_clusters = levels(df$Cluster)) {
  data_t <- data.frame(t(data2@ndata), 'ID' = colnames(data2@ndata))
  
  data_long <- data_t %>% gather('Gene', 'Expression', -ID) %>%
    left_join(df[,c(1,4)])
  
  #order clusters based on hierarchical clustering
  data_long$Cluster <- factor(data_long$Cluster, levels = order_clusters)
  
  return(data_long)
}


#line plot
gene_line_plot <- function(data = data_long, gene) {
  line_plot <- ggplot(data[data$Gene %in% gene,], aes(x=cell_ID, y=Expression, color = Cluster, fill = Cluster)) +
    geom_bar(stat = 'identity') + #width = 0.1, 
    facet_grid(facets = ~Cluster, 
               drop = TRUE, 
               #space = "free", 
               scales = "free", 
               switch = "x",
               space = "free_x") +
    labs(title = gene, y = 'Gene Expression', x = 'Cluster') +
    theme_minimal() +
    theme(axis.line = element_blank(), 
          #axis.title.y = element_blank(),
          #axis.title.x = element_blank(),
          axis.ticks.y = element_blank(), 
          strip.text.x = element_text(),
          #axis.text.y = element_text(size = 10), 
          axis.text.x = element_blank(),#element_text(size = cex.col), 
          strip.background = element_blank(), 
          panel.grid = element_blank(),
          panel.spacing.x = unit(0.2, units = 'line'),
          legend.position = 'None') +
    scale_fill_manual(values =c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",'#984EA3')) +
    scale_color_manual(values = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",'#984EA3'))
  print(line_plot)
  
}

gene_line_plot_go <- function(data = go_term_exp, gene) {
  line_plot <- ggplot(data, aes(x=ID, y=Expression, color = Cluster, fill = Cluster)) +
    geom_bar(stat = 'identity') + #width = 0.1, 
    facet_grid(facets = ~Cluster, 
               drop = TRUE, 
               #space = "free", 
               scales = "free", 
               switch = "x",
               space = "free_x") +
    labs(title = gene, y = 'Gene Expression', x = 'Cluster') +
    theme_minimal() +
    theme(axis.line = element_blank(), 
          #axis.title.y = element_blank(),
          #axis.title.x = element_blank(),
          axis.ticks.y = element_blank(), 
          strip.text.x = element_text(),
          #axis.text.y = element_text(size = 10), 
          axis.text.x = element_blank(),#element_text(size = cex.col), 
          strip.background = element_blank(), 
          panel.grid = element_blank(),
          panel.spacing.x = unit(0.2, units = 'line'),
          legend.position = 'None') +
    scale_fill_manual(values =c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",'#984EA3')) +
    scale_color_manual(values = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",'#984EA3'))
  print(line_plot)
  
}
#generate cumulative gene expression values from enrich_up file

go_term_gene_exp <- function(data1 = enrich_up, data2 = sc, data3 = df) {
  data_t <- data.frame(t(data2@ndata), 'ID' = colnames(data2@ndata))
  data <- data.frame(row.names = colnames(data2@ndata))
  data1 <- na.omit(data1[!duplicated(data1$Description),])
  
  for (i in 1:nrow(data1)) {
    tryCatch({
      n=data1[i,'Description']
      data_ <- dplyr::select(data_t, matches(paste0("^(", c(gsub('/', '|', data1[i,'geneID'])), ")")))
      data_ <- data.frame(rowSums(data_))
      colnames(data_) <- as.character(n)
      data <- cbind(data, data_)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  #add cluster and other data
  data$Cluster <- sc@cpart
  data$Cluster <- factor(data$Cluster, levels = levels(df$Cluster))
  data$ID <- rownames(data)
  data <- na.omit(data)
  data <- data[order(data$Cluster),]
  
  colnames(data) <- gsub("/", "_", colnames(data))
  return(data)
}



#other functions
id2name <- function(x) sub("\\_\\_chr\\w+","",x)

name2id <- function(x,id) {
  ##  id[sub("\\_\\_chr\\w+","",id) %in% x]
  n <- c()
  for ( j in x ){ n <- append(n,id[grep(paste(j,"(\\_\\_chr\\w+|$|\\|)",sep=""),id)])
  }
  n
}

#loading genes
load_data <- function(path) { 
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)
  tables <- lapply(files, read.csv)
  do.call(rbind, tables)
}

#mean expression heatmaps
plot_heatmap <- function(.sc = counts_control, .df=df_control, .up_genes= up_genes_control, .colors = viridis(100)) {
  gene_names <- .up_genes %>%
    filter(padj<0.05, log2FoldChange >1) %>%
    group_by(Cluster) %>%
    dplyr::arrange(Cluster, padj) %>%
    dplyr::distinct(Cluster, GENEID, .keep_all = T) %>% #from url: https://dplyr.tidyverse.org/reference/distinct.html
    top_n(n = 20, wt=log2FoldChange) 
  
  gene_names <- unique(as.character(gene_names$GENEID))
  
  nd <- .sc[gene_names,]
  nd[is.na(nd)] <- 0
  nd  <- t(nd[complete.cases(nd),])
  clust_n <- as.numeric(as.character(.df$Cluster))
  mnd <- as.data.frame(cbind(clust_n,nd[rownames(nd) %in% .df$ID,]))
  mean_mnd <- aggregate(mnd[, 2:dim(mnd)[2]], list(mnd$clust_n), mean)
  row.names(mean_mnd) <- paste("C",mean_mnd$Group.1,sep = "")
  mean_mnd$Group.1 <- NULL
  gene <- as.data.frame(t(log2(mean_mnd)))
  gene <- gene[complete.cases(gene),]
  row.names(gene) <- id2name(rownames(gene))
  pheatmap(gene, cluster_cols = T, color=.colors, cluster_rows=T,fontsize_row = 8, border_color = F, show_rownames = T,show_colnames = T, scale = "row")
  
  
}

#pie_chart
pie_chart <- function(.df=df, .colors=colors_fig) {
  df2 <- .df %>% group_by(Cluster) %>%
    summarise(Freq = n()) %>%
    dplyr::mutate(Ratio = Freq/sum(Freq))
  
  pie_plot <- ggplot(df2, aes(x=factor(1), y = Ratio,fill = Cluster)) +
    geom_bar(position = 'fill', stat = 'identity', width = 1, color='black', lwd=0.1) +
    coord_polar(theta='y') +
    theme(panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          text=element_text(size=17),
          axis.text.x=element_blank(), 
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.title = element_blank(),
          legend.key = element_blank()) +
    scale_fill_manual(values = .colors)
  pie_plot 
  
}
