#open needed packages
library(readxl)
library(httr)
library(tidyverse)
library(pheatmap)
library(MASS)

#load RaceID and custom functions
source(file.path("R", "RaceID3_StemID2_class.R"))
source(file.path("R", "sankowski-et-al-functions.R"))

#load sc data
load(file.path("data", "sc_gbm_microglia_nn.RData"))

#set cluster order
if (!file.exists(file.path("data","order_clusters_gbm_micr.RData"))) {
  order_clusters <- clustheatmap(sc,final=TRUE)
  save(order_clusters, file=file.path("data", "order_clusters_gbm_micr.RData"))
} else {
  load(file.path("data", "order_clusters_gbm_micr.RData"))
}

#exclude small clusters
if (!file.exists(file.path("data", "retain_cl_gbm_micr.RData"))) {
  cell_numbers <-as.numeric()
  for (i in 1:length(unique(sc@cpart)))
  {
    cell_numbers[i] <- length(sc@cpart[sc@cpart==i])
  }
  names(cell_numbers) <- c(1:length(unique(sc@cpart)))
  retain_cl <- as.numeric(names(cell_numbers[cell_numbers > dim(sc@ndata)[2]/100]))
  
  save(retain_cl, file = file.path("data", "retain_cl_gbm_micr.RData"))
  
} else {
  load(file.path("data", "retain_cl_gbm_micr.RData"))
} 

order_clusters <- order_clusters[order_clusters %in% retain_cl]

#build a table with cell Ids, cluster and Patient_IDs
df <- data.frame('Cluster' = sc@cpart, sc@tsne) 
df$ID <- rownames(df)
df$Patient_ID <- gsub('_.*', '', rownames(df))
df$Patient_ID2 <- factor(gsub("H", "Ctrl", df$Patient_ID), levels = c(paste0(c("Ctrl", "GBM"), rep(1:4, each =2))))
df <- df[df$Cluster %in% retain_cl,] 
df$Cluster <- factor(df$Cluster, levels = order_clusters)
df$Diagnosis <- factor(ifelse(grepl('GBM', df$Patient_ID), 'GBM', 'Control'), levels = c("GBM", "Control"))

#load batch info from the supplementary information on the paper website https://www.nature.com/articles/s41593-019-0532-y#Sec30
url1 <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-019-0532-y/MediaObjects/41593_2019_532_MOESM3_ESM.xlsx'
GET(url1, write_disk(tf <- tempfile(fileext = ".xlsx")))
batch_info <- read_excel(tf,sheet = 8, range = "B1:L9")
str(batch_info)

colnames(batch_info)[1] <- "Patient_ID"
df <- left_join(df, batch_info)

#export metadata
save(df, file = file.path("data", "metadata_ctrl_gbm.RData"))

# Figure 4b - cluster t-SNE
tsne_plot_no_outline(FILL=df$Cluster, fill_colors = colors_fig, point_size = 4) 

# Figure 4c - diagnoses t-SNE
tsne_plot_no_outline(FILL=df$Diagnosis, fill_colors = toupper(c('#f1a340','#998ec3')), point_size = 4) 

# Figure 4d - heatmap of cluster marker genes
up_genes <- read.csv(file.path('data', 'Human_gbm_microglia_up_genes_padj<05_logfc>1.csv'), row.names = 1)
#filter unannotated genes and genes associated with tissue dissection
up_genes2 <- up_genes[!grepl("^(HSPA1A|MTRNR2L8|MTRNR2L12|HSP90AA1|MALAT1|ZFP36L1|ZFP36|FOS|MALAT1|HSPB*|DUSP|HSPH1|HSPA*|JUN|HSP90B1|RPS16|DNAJB1|H3F3B|HERPUD1|NEAT1|IVNS1ABP|HIST1H2BG|RP*|XIST)", up_genes$GENEID),  ]

gene_names <- up_genes2 %>%
  filter(padj<0.05, log2FoldChange >1) %>%
  group_by(Cluster) %>%
  dplyr::arrange(Cluster, padj) %>%
  dplyr::distinct(Cluster, GENEID, .keep_all = T) %>% #from url: https://dplyr.tidyverse.org/reference/distinct.html
  top_n(n = 20, wt=log2FoldChange) %>%
  ungroup

#filter out other unannotated/non-protein coding genes
gene_names <- gene_names %>%
  filter(!grepl("^(CTD-|CTC-|CTA-)", .$GENEID))
plot_heatmap(sc, .gene_names = gene_names)

# Figure 4e
mosaicGG2(df, "Cluster", "Patient_ID2") +
  scale_fill_brewer(palette = "Paired")

# Figure 4f
mosaicGG2(df, "Cluster", "Diagnosis", colors = toupper(c('#f1a340','#998ec3'))) ; ggsave("plots/cluster_diagnosis_marimekko.pdf", useDingbats=F)

# Figure 4g
genes <- c("CX3CR1__chr3", "VEGFA__chr6", "IFI44__chr1")

#gene expression tsne plots
plotexptsne2(.gene=genes[1], point_size = 6)
plotexptsne2(.gene=genes[2], point_size = 6)
plotexptsne2(.gene=genes[3], point_size = 6)

#gene expression line plots
#for this plot one needs a long data plot
data_long <- make_data_long(.genes = genes) %>% na.omit
#generate selected line plots
gene_line_plot(.gene=genes[1])
gene_line_plot(.gene=genes[2])
gene_line_plot(.gene=genes[3])

# Figure 5d
genes <- c("^AIF1", "P2RY12", "TMEM119", "HLA-DRA", "CD68", "CD74", "^SPP1")
plotexptsne2(.gene=genes[1], point_size = 6)
plotexptsne2(.gene=genes[2], point_size = 6)
plotexptsne2(.gene=genes[3], point_size = 6)
plotexptsne2(.gene=genes[4], point_size = 6)
plotexptsne2(.gene=genes[5], point_size = 6)
plotexptsne2(.gene=genes[6], point_size = 6)

# Figure 5f
spp1 <- read_csv2(file.path("data","counting_GBM_SPP1.csv")) %>%
  mutate(Localization = gsub("WM", "Ctrl", .$Localization)) %>%
  filter(Condition == "spp1_pos") %>%
  dplyr::select(ID, Localization, Percentage) 
spp1 %>%
  ggplot(aes(Localization, Percentage, fill=Localization)) +
  geom_boxplot(width=.4) +
  geom_jitter(width = .2, height = 0.01, pch=21, size = 9, stroke = 0.25) +
  scale_fill_manual(guide=F, values = toupper(c('#998ec3','#f1a340'))) +
  theme_minimal() +
  theme(text = element_text(size =30)) +
  labs(y= "SPP1+ cells (%)", x=element_blank()) 

# Figure 5g
plotexptsne2(.gene=genes[7], point_size = 6)

# Extended data figure 9b
clustheatmap(sc, final = T)

#Extended data figure 10
genes <- c("P2RY12__chr3", "SELPLG__chr12", "CD163__chr12", "APOE__chr19", "HIF1A__chr14")

#gene expression tsne plots
plotexptsne2(.gene=genes[1], point_size = 6)
plotexptsne2(.gene=genes[2], point_size = 6)
plotexptsne2(.gene=genes[3], point_size = 6)
plotexptsne2(.gene=genes[4], point_size = 6)
plotexptsne2(.gene=genes[5], point_size = 6)

#gene expression line plots
#for this plot one needs a long data plot
data_long <- make_data_long(.genes = genes) %>% na.omit
#generate selected line plots
gene_line_plot(.gene=genes[1])
gene_line_plot(.gene=genes[2])
gene_line_plot(.gene=genes[3])
gene_line_plot(.gene=genes[4])
gene_line_plot(.gene=genes[5])

# Supplementary Figure 1
data_t <- as.data.frame(t(as.matrix(sc@ndata)))
data_t$ID <- rownames(data_t)
data_t <- data_t %>% left_join(df)
data_t$Diagnosis <- factor(data_t$Diagnosis, levels = c('Control', 'GBM'))

genes <- c("CSF1R__chr5", "P2RY12__chr3",  "HLA-DRA__chr6", "VEGFA__chr6", "HIF1A__chr14", "CD163__chr12", "SPP1__chr4", "APOE__chr19", "IFI44__chr1")

for (gene in genes) {
  data <- data.frame(Cluster = data_t[["Cluster"]], Expression=data_t[[gene]], Diagnosis=data_t[["Diagnosis"]])
  plt <- ggplot(na.omit(data), aes(Cluster, Expression, fill = Diagnosis)) +
  geom_violin(scale = 'width', lwd=0.25) +
  geom_boxplot(width=0.5, outlier.shape = NA, position=position_dodge(1))+
  scale_fill_manual(values = toupper(c('#998ec3','#f1a340'))) +
  theme_minimal() +
  labs(y='Gene expression', title=gene) 

  print(plt)
  ggsave(file.path("plots", paste0(gene,"boxplots.pdf")), useDingbats=F)
}

