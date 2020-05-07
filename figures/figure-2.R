#figure 2 plotting functions

#load packages
require(tidyverse)
require(clusterProfiler)
library(org.Hs.eg.db)

#load data
load("data/counts_control.RData")
load("data/df_control.RData")

#go term analysis
              load("data/up_genes_control.RData") 
              up_genes <- up_genes_control %>%
                filter(padj<.05, log2FoldChange>1) %>%
                distinct(GENEID, Cluster) %>% 
                mutate(GENEID = gsub('_.*', '', as.character(GENEID)))
              
              #run cluster profiler
              
              #load background genes
              back_genes <- rownames(counts_control)[which(apply(counts_control > .1, 1, sum)>0)]
              background <- bitr(sub('_.*', '',back_genes), fromType = "SYMBOL",
                                 toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                                 OrgDb = org.Hs.eg.db)
              
              #define empty data frame to collect data
              enrich_up <- data.frame(matrix(ncol = 10))
              colnames(enrich_up) <- c('ID','Description', 'GeneRatio', 'BgRatio' ,'pvalue', 'p.adjust', 'qvalue', 'geneID','Count' , 'Cluster')
              
              #this will take a while
              for (i in unique(up_genes$Cluster))  {
                
                tryCatch({
                  gene <- up_genes$GENEID[up_genes$Cluster == i]
                  gene.up_genes <- bitr(gene, fromType = "SYMBOL",
                                  toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                                  OrgDb = org.Hs.eg.db)
                  head(gene.up_genes)
                  
                  
                  
                  
                  ggo <- groupGO(gene     = gene.up_genes[,3],
                                 OrgDb    = org.Hs.eg.db,
                                 ont      = "BP",
                                 level    = 3,
                                 readable = TRUE)
                  
                  head(ggo)
                  
                  ego <- enrichGO(gene          = gene.up_genes[,3],
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

              head(enrich_up)

#plot go terms
              enrich_up <- read.csv('data/GO-terms-bp-control.csv', stringsAsFactors = F) %>%
                na.omit() %>%
                #extract biologically relevant go terms
                filter(Description %in% read_csv("data/biologically-relevant-bp-GO-terms-control.csv")[[1]])
              
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
              
              enrich_up <- enrich_up[enrich_up$Cluster %in% retain_cl,]
              ord_clust <- clustheatmap(sc, final = T, hmethod = "single")
              ord_clust <- ord_clust[ord_clust %in% retain_cl]
              enrich_up$Cluster <- factor(enrich_up$Cluster, levels = ord_clust)
              enrich_up$Description <- reorder(enrich_up$Description,enrich_up$Description,FUN=length)
              
              #re-order enrich_up based on the levels of Cluster
              enrich_up <- enrich_up[with(enrich_up, order(Cluster)),] #from url: https://stackoverflow.com/questions/1296646/how-to-sort-a-dataframe-by-columns
              enrich_up$Description <- factor(enrich_up$Description, levels = rev(enrich_up$Description[!duplicated(enrich_up$Description)]))
              colnames(enrich_up)[9] <- 'GeneCount'
              
              dot_plot <- ggplot(enrich_up, aes(Cluster, Description, size = GeneCount, fill= -log2(qvalue))) + #[enrich_up$GeneCount>4,]
                geom_point(pch=21, stroke=0.25) +
                scale_fill_viridis() +
                theme_light() +
                theme(text=element_text(size=10),
                      axis.title.y=element_blank())
              dot_plot
              
              ggsave('bp/20180616-top40-PC_all_GoTerm_dot_plot.pdf', height = 7.5, width = 6, units = 'in')
              
              svg('bp/20180616-top40-PC_all_GoTerm_dot_plot.svg', height = 9, width = 6.5) #, units = 'in', res = 300
              dot_plot #https://stackoverflow.com/questions/15678261/r-ggplot-does-not-work-if-it-is-inside-a-for-loop-although-it-works-outside-of
              dev.off()