#monocle analysis of healthy vs gbm data
#source("https://bioconductor.org/biocLite.R")
#biocLite('monocle')

#load RaceID and custom functions
source(file.path("R", "RaceID3_StemID2_class.R"))
source(file.path("R", "sankowski-et-al-functions.R"))

library(tidyverse)
library(monocle)
library(viridis)

date <- Sys.Date()

#load sc file
load(file.path("data", "sc_gbm_microglia_nn.RData"))

#load metadata
load(file.path("data", "metadata_ctrl_gbm.RData"))

#set up the pdata dataset
colnames(df)[1] <- 'RaceID_Cl'
rownames(df) <- colnames(sc@expdata)[colnames(sc@expdata) %in% df$ID]
feat_df <- data.frame('gene_short_name' = rownames(sc@expdata))
rownames(feat_df) <- rownames(sc@expdata)

cds <- newCellDataSet(as.matrix(sc@expdata[colnames(sc@expdata) %in% df$ID]),
                      phenoData = new("AnnotatedDataFrame", data = df),
                      featureData = new("AnnotatedDataFrame", data = feat_df),
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

#estimate dispersion
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

#filter cells
cds <- detectGenes(cds, min_expr = 0.5)
print(head(fData(cds)))

expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 50))
expressed_genes  <- expressed_genes[!grepl("^(IER|EGR|CTD-|HSPA1A|MTRNR2L8|MTRNR2L12|HSP90AA1|MALAT1|ZFP36L1|ZFP36|FOS|MALAT1|HSPB*|DUSP|HSPH1|HSPA*|JUN|HSP90B1|RPS16|DNAJB1|H3F3B|HERPUD1|NEAT1|IVNS1ABP|HIST1H2BG|RP*|XIST)", expressed_genes)]


#cluster cells
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.5)
HSMM <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds)


#variance explained
cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 11)

#plots
plot_cell_clusters(cds, 1, 2, color = "Cluster")
plot_cell_clusters(cds, 1, 2, color = "Diagnosis")
plot_cell_clusters(cds, 1, 2, color = "RaceID_Cl")

save(cds, file = file.path("data", "monogle-gbm-data.RData"))

#construct ptrajectories
diff_test_res <- differentialGeneTest(cds[expressed_genes,],
                                      fullModelFormulaStr = "~Diagnosis")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds<- orderCells(cds)

#Suplementary Figure 6a
plot_cell_trajectory(cds, color_by = "Diagnosis") +
  scale_color_manual(values = toupper(c('#f1a340','#998ec3')))
plot_cell_trajectory(cds, color_by = "RaceID_Cl")

save(cds, file = file.path("data","-trajectories-monogle-gbm-data.RData" ))
write.csv(diff_test_res, file.path("data", "-monocle-differentially-expr-genes.csv"))

#plot marker genes
gene <- unlist(sc@expdata['TMEM119__chr12', df$ID])
gene <- (gene-min(gene))/(max(gene)-min(gene))
gene <- unlist(sc@expdata['CCL3L3__chr17|CCL3__chr17', df$ID])
plot_cell_trajectory(cds, color_by = gene) +
  scale_color_viridis()
gene <- unlist(sc@expdata['P2RY13__chr3', df$ID])
plot_cell_trajectory(cds, color_by = gene) +
  scale_color_viridis()

# examine the branches of the trajectories
diffgenes2 <- read.csv(file.path('data', 'Human_gbm_microglia_up_genes_padj<05_logfc>1.csv'), row.names = 1)
diffgenes2  <- diffgenes2 [!grepl("^(IER|EGR|CTD-|HSPA1A|MTRNR2L8|MTRNR2L12|HSP90AA1|MALAT1|ZFP36L1|ZFP36|FOS|MALAT1|HSPB*|DUSP|HSPH1|HSPA*|JUN|HSP90B1|RPS16|DNAJB1|H3F3B|HERPUD1|NEAT1|IVNS1ABP|HIST1H2BG|RP*|XIST)", diffgenes2 )]

homeo <- unique(diffgenes2$GENEID[diffgenes2$Cluster %in% c(11,12)])
proinf <- unique(diffgenes2$GENEID[diffgenes2$Cluster %in% c(1:4)])
gbm <- unique(diffgenes2$GENEID[diffgenes2$Cluster %in% c(13,14)])
aging<- unique(diffgenes2$GENEID[diffgenes2$Cluster %in% c(5)])

gene <- colSums(na.omit(sc@ndata[homeo, df$ID]))

plot_cell_trajectory(cds, color_by = gene) + 
  scale_color_gradientn(colors = colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100)) +
  labs(title='homeostatic')

gene <- colSums(na.omit(sc@expdata[proinf, df$ID]))

plot_cell_trajectory(cds, color_by = gene) + 
  scale_color_gradientn(colors = colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100)) +
  labs(title='proinflammatory')

gene <- colSums(na.omit(sc@ndata[aging, df$ID]))

plot_cell_trajectory(cds, color_by = gene) + 
  scale_color_gradientn(colors = colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100))+
  labs(title='aging')

gene <- colSums(na.omit(sc@ndata[gbm, df$ID]))

plot_cell_trajectory(cds, color_by = gene) + 
  scale_color_gradientn(colors = colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100))+
  labs(title='gbm')

#examine the expression of individual genes
for (i in unique(diffgenes2$GENEID)) {
  gene <- colSums(na.omit(sc@ndata[i, df$ID]))
  plt <- plot_cell_trajectory(cds, color_by = gene) + 
    scale_color_gradientn(colors = colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100))+
    labs(title='gbm')
  print(plt)
  
}

#Suplementary Figure 6b
diffgenes <- read_csv(file.path("data", "Human_gbm_micr_unique_up_genes_wo_stress-genes_log2fc>1.csv"))$GENEID
diffgenes  <- diffgenes[!grepl("^(IER|EGR|CTA-|CTD-|MTRNR2L8|MTRNR2L12|HSP90AA1|MALAT1|ZFP36L1|ZFP36|FOS|MALAT1|HSPB*|DUSP|HSPH1|HSPA*|JUN|HSP90B1|RPS16|DNAJB1|H3F3B|HERPUD1|NEAT1|IVNS1ABP|HIST1H2BG|RP*|XIST)", diffgenes)]
marker_genes <- row.names(subset(fData(cds),
                                 gene_short_name %in% diffgenes))

diff_test_res <- differentialGeneTest(cds[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))

plot_pseudotime_heatmap(cds[sig_gene_names,],
                        num_clusters = 5,
                        cores = 1,
                        show_rownames = T)
