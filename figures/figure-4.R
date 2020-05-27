#figure 4 plotting functions

#load functions
source("functions/sankowski-et-al-functions.R")

#load data
load("data/metadata-ctrl-gbm.RData")
load("data/unnormalized_counts_ctrl_gbm.RData")
load("data/counts-ctrl-gbm.RData")
#load go terms
load("data/GO-terms-bp-gbm-description-unchanged.RData")

#make df$Diagnosis a factor
df$Diagnosis <- factor(df$Diagnosis, levels = c("GBM","Control"))
df$Patient_ID <- factor(df$Patient_ID, levels = paste0(c("Ctrl", "GBM"), rep(1:4, each=2)))

#Fig 4b
#cluster
tsne_plot(FILL=df$Cluster) +
  scale_fill_manual(values=colors_fig)

#Fig 4c
#diagnoses
tsne_plot(FILL=df$Diagnosis) +
  scale_fill_manual(values = toupper(c('#f1a340','#998ec3')))

#Fig 4d
#heatmap
mean_heatmap(sc, .genes = name2id(c("HIF1A", "CX3CR1", "TMEM119", "CD163"), rownames(sc@ndata)),.retain_cl = unique(sc@cpart))

#Fig 4e
#patient id
mosaicGG2(df, "Cluster", "Patient_ID") +
  scale_fill_brewer(palette = "Paired")

#Fig 4f
#marimekko plots
#Diagnosis
mosaicGG2(df, X="Cluster",FILL = "Diagnosis") +
  scale_fill_manual(values = toupper(c('#f1a340','#998ec3')))

#Fig 4g
#tsne plots with gene expression
#cx3cr1
plotexptsne2("CX3CR1")

#vegf
plotexptsne2("VEGFA")

#hif1a
plotexptsne2("IFI44")

#line plots
counts2 <- counts[ c("CX3CR1", "VEGFA", "IFI44"),]
data_long <- t(counts2) %>%
  as.data.frame() %>%
  rownames_to_column(var = "cell_ID") %>%
  pivot_longer(CX3CR1: IFI44, names_to = "Gene", values_to = "Expression") %>%
  right_join(df)

gene_line_plot(data_long, gene = "CX3CR1")
gene_line_plot(data_long, gene = "VEGFA")
gene_line_plot(data_long, gene = "IFI44")

#Fig 4h
#plot go terms
genes <- c(strsplit(enrich_up[enrich_up$ID == "GO:0010573",'geneID'], split ='/')[[1]])
counts2 <- counts[genes,]
plotexptsne2(gene = genes, .sc=counts2, line_width = 0)

genes <- c(strsplit(enrich_up[enrich_up$ID == "GO:0002480",'geneID'], split ='/')[[1]])
counts2 <- counts[genes,]
plotexptsne2(gene = genes, .sc=counts2, line_width = 0)


#suppl figure 2a
load("data/GO-terms-bp-gbm.RData")
dot_plot <- ggplot(enrich_up, aes(Cluster, Description, size = GeneCount, fill= -log2(qvalue))) +
  geom_point(pch=21, stroke=0.25) +
  scale_fill_viridis() +
  theme_light() +
  theme(text=element_text(size=10),
        axis.title.y=element_blank())
dot_plot

#suppl figure 2b
