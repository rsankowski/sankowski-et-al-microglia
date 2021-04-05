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
load(file.path("data", "sc_gbm_microglia_nn.RData"))
load(file.path("data", "order_clusters_gbm_micr.RData"))
load(file.path("data", "retain_cl_gbm_micr.RData"))

df <- read.csv(file.path('data', 'Human_gbm_microglia_up_genes_padj<05_logfc>1.csv'), row.names = 1)
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
    
    ggo <- groupGO(gene     = gene.df[,3],
                   OrgDb    = org.Hs.eg.db,
                   ont      = "BP",
                   level    = 3,
                   readable = TRUE)
    
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
    
    ego_simpl <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = "Wang")
    head(ego_simpl@result, 50)
    ego_simpl2 <- ego_simpl[!duplicated(ego_simpl@result$geneID)]
    head(ego_simpl2, 50)
    ego_simpl2$Cluster <- rep(as.character(i, nrow(ego_simpl2)))
    enrich_up <- rbind(enrich_up, ego_simpl2)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.csv(enrich_up, file.path("data", "GO_terms_bp_ctrl_vs_gbm_micr.csv"))

# Supplementary Figure 2a
# load GO term enrichment output from the paper supplementary
url1 <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-019-0532-y/MediaObjects/41593_2019_532_MOESM3_ESM.xlsx'
GET(url1, write_disk(tf <- tempfile(fileext = ".xlsx")))
enrich_up <- read_excel(tf,sheet = 10, range = "B1:K453")

#load selected GO terms 
relevant <- read.csv(file.path("data", "relevant-GO-terms-gbm.csv"), stringsAsFactors = F)[[1]]
#Extract selected GO terms
enrich_up <- enrich_up[enrich_up$Description %in% relevant,]

#collapse similar GO terms
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
enrich_up$GeneCount <- as.numeric(enrich_up$GeneCount)
enrich_up$qvalue <- as.numeric(enrich_up$qvalue)

#plot
go_dot_plot()

# Figure 4h
#load the GO enrichment terms from the paper supplementary 
url1 <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-019-0532-y/MediaObjects/41593_2019_532_MOESM3_ESM.xlsx'
GET(url1, write_disk(tf <- tempfile(fileext = ".xlsx")))
enrich_up <- read_excel(tf,sheet = 10, range = "B1:K453")

#load metadata gbm
load(file.path("data", "metadata_ctrl_gbm.RData"))

genes_vegf <- c("C5AR1", "IL6ST", "IL1B")
plotexptsne2(.gene=genes_vegf, point_size = 6) +
  labs(subtitle = "GO:0010573")

genes_antigen_pr <- c("HLA-A", "HLA-B", "HLA-C")
plotexptsne2(.gene=genes_antigen_pr, point_size = 6) +
  labs(subtitle = "GO:0002480")

# Supplementary Figure 2b - line plots
data_long <- make_data_long_cumulative(.genes = genes_antigen_pr, .order_clusters = order_clusters) %>%
  na.omit()
gene_line_plot(.gene="cumulative") +
  labs(title = "antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-independent")

data_long <- make_data_long_cumulative(.genes = genes_vegf, .order_clusters = order_clusters) %>%
  na.omit()
gene_line_plot(.gene="cumulative") +
  labs(title = "vascular endothelial growth factor production")
