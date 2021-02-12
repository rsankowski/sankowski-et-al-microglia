library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
keytypes(org.Hs.eg.db)
library(viridis)
library(pheatmap)

#load RaceID and custom functions
source(file.path("R", "RaceID3_StemID2_class.R"))
source(file.path("R", "sankowski-et-al-functions.R"))

#load sc data
load(file.path("data", "sc_ctrl_microglia_nn.RData"))
load(file.path("data", "retain_cl_ctrl.RData"))
load(file.path("data", "order_clustersers_ctrl.RData"))

df <- read.csv(file.path('data','Human_ctrl_micr_unique_up_genes_wo_stress-genes_padj<05.csv'), stringsAsFactors = F, row.names = 1)
df$GENEID <- gsub('_.*', '', df$GENEID)

#load background genes
back_genes <- rownames(sc@ndata)[which(apply(sc@ndata > .1, 1, sum)>0)]
background <- bitr(sub('_.*', '',back_genes), fromType = "SYMBOL",
                   toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                   OrgDb = org.Hs.eg.db)

#define empty data frame to collect data
enrich_up <- data.frame(matrix(ncol = 10))
colnames(enrich_up) <- c('ID','Description', 'GeneRatio', 'BgRatio' ,'pvalue', 'p.adjust', 'qvalue', 'geneID','Count' , 'Cluster')

for (i in unique(df$Cluster))  {
  
  tryCatch({
    gene <- df$GENEID[df$Cluster == i]
    gene.df <- bitr(gene, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                    OrgDb = org.Hs.eg.db)
    head(gene.df)
    
    
    
    
    ggo <- groupGO(gene     = gene.df[,3],
                   OrgDb    = org.Hs.eg.db,
                   ont      = "BP",
                   level    = 3,
                   readable = TRUE)
    
    head(ggo)
    
    ego <- enrichGO(gene          = gene.df[,3],
                    universe      = background[,3],
                    OrgDb         = org.Hs.eg.db,
                    minGSSize     = 1,
                    ont           = "BP",
                    pool          = TRUE,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    head(ego@result, 50)
    
    ego_simpl <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = "Wang")
    head(ego_simpl@result, 50)
    ego_simpl2 <- ego_simpl[!duplicated(ego_simpl@result$geneID)]
    head(ego_simpl2, 50)
    ego_simpl2$Cluster <- rep(as.character(i, nrow(ego_simpl2)))
    enrich_up <- rbind(enrich_up, ego_simpl2)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.csv(na.omit(enrich_up), file.path("data", 'GO_terms_bp_ctrl_micr.csv'))

# Figure 2a
#edit the full list of padj<.05 GO terms
enrich_up <- read.csv(file.path("data", 'GO_terms_bp_ctrl_micr.csv'), stringsAsFactors = F)[,-1]
table(enrich_up$Description)

#load selected GO terms and collapse similar terms
relevant <- read.csv(file.path("data", "20181031-relevant-bp-GO-terms.csv"), stringsAsFactors = F)[[1]]
#Extract relevant GO terms
enrich_up <- enrich_up[enrich_up$Description %in% relevant,]

#collapse duplicate GO terms
enrich_up$Description <- ifelse(grepl('^(chemotaxis|cell chemotaxis|regulation of leukocyte chemotaxis|regulation of chemotaxis)', enrich_up$Description), '(regulation of leukocyte) chemotaxis', 
                                ifelse(grepl('interferon-gamma', enrich_up$Description), '(cellular) response to interferon-gamma', 
                                       ifelse(grepl('T cell activation', enrich_up$Description), '(positive regulation of) T cell activation', 
                                              ifelse(grepl('^(positive regulation of leukocyte migration|mononuclear cell migration|leukocyte migration|regulation of leukocyte migration|macrophage migration)', enrich_up$Description), 'Cell migration', 
                                                     ifelse(grepl('cell-cell adhesion', enrich_up$Description), '(regulation of) leukocyte cell-cell adhesion', 
                                                            ifelse(grepl('regulation of inflammatory response', enrich_up$Description), '(positive) regulation of inflammatory response', 
                                                                   ifelse(grepl('regulation of vasculature development|vascular endothelial growth factor production|angiogenesis', enrich_up$Description), '(positive) regulation of vasculature development', 
                                                                          ifelse(grepl('cytokine secretion', enrich_up$Description), '(positive) regulation of cytokine secretion', 
                                                                                 ifelse(grepl('positive regulation of response to external stimulus|cellular response to external stimulus', enrich_up$Description), '(positive) regulation of response to external stimulus', 
                                                                                        ifelse(grepl('MHC class II', enrich_up$Description), 'antigen processing and presentation of \npeptide/polysaccharide antigen via MHC class II', 
                                                                                               ifelse(grepl('phagocytosis', enrich_up$Description), 'phagocytosis', 
                                                                                                      ifelse(grepl('response to oxidative stress', enrich_up$Description), 'response to oxidative stress', 
                                                                                                             ifelse(grepl('chemokine production', enrich_up$Description), 'chemokine production', 
                                                                                                                    ifelse(grepl('humoral immune response', enrich_up$Description), '(positive regulation of) humoral immune response', 
                                                                                                                           ifelse(grepl('complement activation', enrich_up$Description), '(regulation of) complement activation', 
                                                                                                                                  ifelse(grepl('negative regulation of collateral sprouting|collateral sprouting of intact axon in response to injury', enrich_up$Description), 'response to axon injury',  enrich_up$Description))))))))))))))))





#introduce line breaks
enrich_up$Description2 <- enrich_up$Description
enrich_up$Description <- gsub('exogenous peptide antigen via MHC class I,', 'exogenous \npeptide antigen via MHC class I', enrich_up$Description)

#filter out duplicated terms
enrich_up <- enrich_up[!duplicated(enrich_up[,c('Cluster', 'Description')]),]

#sort the terms
enrich_up <- enrich_up[enrich_up$Cluster %in% retain_cl,]
order_clusters <- order_clusters[order_clusters %in% retain_cl]
enrich_up$Cluster <- factor(enrich_up$Cluster, levels = order_clusters)
enrich_up$Description <- reorder(enrich_up$Description,enrich_up$Description,FUN=length)

#re-order enrich_up based on the levels of Cluster
enrich_up <- enrich_up[with(enrich_up, order(Cluster)),] #from url: https://stackoverflow.com/questions/1296646/how-to-sort-a-dataframe-by-columns
enrich_up$Description <- factor(enrich_up$Description, levels = rev(enrich_up$Description[!duplicated(enrich_up$Description)]))
colnames(enrich_up)[9] <- 'GeneCount'


go_dot_plot()

