#open needed packages
library(readxl)
library(httr)
library(tidyverse)
library(pheatmap)
library(viridis)

#load RaceID and custom functions
source(file.path("R", "RaceID3_StemID2_class.R"))
source(file.path("R", "sankowski-et-al-functions.R"))

#load sc data
load(file.path("data", "sc_ctrl_microglia_nn.RData"))

#set cluster order
if (!file.exists(file.path("data","order_clusters_ctrl.RData"))) {
  order_clusters <- clustheatmap(sc,final=TRUE)
  save(order_clusters, file=file.path("data", "order_clusters_ctrl.RData"))
} else {
  load(file.path("data", "order_clusters_ctrl.RData"))
}

#exclude small clusters
if (!file.exists(file.path("data", "retain_cl_ctrl.RData"))) {
cell_numbers <-as.numeric()
for (i in 1:length(unique(sc@cpart)))
{
  cell_numbers[i] <- length(sc@cpart[sc@cpart==i])
}
names(cell_numbers) <- c(1:length(unique(sc@cpart)))
retain_cl <- as.numeric(names(cell_numbers[cell_numbers > dim(sc@ndata)[2]/100]))
retain_cl <- retain_cl[retain_cl != 4]

save(retain_cl, file = file.path("data", "retain_cl_ctrl.RData"))

} else {
  load(file.path("data", "retain_cl_ctrl.RData"))
} 
# exclude cluster with high expression of ADAM10 and PIWIL1 #this cluster may be different if the analysis is rerun
order_clusters <- order_clusters[order_clusters %in% retain_cl]

#build a table with cell Ids, cluster and Patient_IDs
df <- data.frame('Cluster' = sc@cpart, sc@tsne) 
df$ID <- rownames(df)
df$cell_ID <- rownames(df)
df$Patient_ID <- gsub('_.*', '', rownames(df))
df <- df[df$Cluster %in% retain_cl,] 
order_clusters <- order_clusters[order_clusters %in% retain_cl]
df$Cluster <- factor(df$Cluster, levels = order_clusters)
#define region
df$Region <- ifelse(grepl('WM',df$ID), 'WM', 
                    ifelse(grepl('GM',df$ID), 'GM', 'Mixed'))

df$Region <- factor(df$Region, levels = c('GM', 'WM', 'Mixed'))

#load batch info from the supplementary information on the paper website https://www.nature.com/articles/s41593-019-0532-y#Sec30
url1 <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-019-0532-y/MediaObjects/41593_2019_532_MOESM3_ESM.xlsx'
GET(url1, write_disk(tf <- tempfile(fileext = ".xlsx")))
batch_info <- read_excel(tf,sheet = 1 )
str(batch_info)

colnames(batch_info)[c(1,7, 12)] <- c("Patient_ID", 'Diagnosis', "Age_bin")
df <- left_join(df, batch_info)

#number of genes per cell
df$Gene_number <- apply(sc@ndata[, df$ID], 2, function(x) sum(x>0.1))

#Plot clusters
df_clusters <- df %>%
  group_by(Cluster) %>%
  summarize(mean_V1 = mean(V1),
            mean_V2 = mean(V2))
df <- df %>% 
  left_join(df_clusters)
save(df, file = file.path("data", "metadata_ctrl.RData"))
# Figure 1c - heatmap of cluster marker genes
up_genes <- read.csv(file.path('data', 'Human_ctrl_microglia_up_genes_padj<05_logfc>1.csv'), row.names = 1)
#filter unannotated genes and genes associated with tissue dissection
up_genes2 <- up_genes[!grepl("^(HSPA1A|MTRNR2L8|MTRNR2L12|HSP90AA1|MALAT1|ZFP36L1|ZFP36|FOS|MALAT1|HSPB*|DUSP1|HSPH1|HSPA*|JUN|HSP90B1|RPS16|DNAJB1|H3F3B|HERPUD1|NEAT1|IVNS1ABP|HIST1H2BG|RP*|XIST)", up_genes$GENEID),  ]

gene_names <- up_genes2 %>%
  filter(padj<0.05, log2FoldChange >1) %>%
  group_by(Cluster) %>%
  dplyr::arrange(Cluster, padj) %>%
  dplyr::distinct(Cluster, GENEID, .keep_all = T) %>% #from url: https://dplyr.tidyverse.org/reference/distinct.html
  top_n(n = 20, wt=log2FoldChange) %>%
  ungroup

#filter out other unannotated/non-protein coding genes
gene_names <- gene_names %>%
  filter(!grepl("^(AC005614.3|CTD-|CTC-|AC074117.10|AC011747.4|LACTB2-AS1)", .$GENEID))
plot_heatmap(sc, .gene_names = gene_names)

# Figure 1d - marimekko plot of patient IDs
df$Patient_ID <- factor(df$Patient_ID, levels = paste0("Pat", 1:15))
mosaicGG2(df, "Cluster", "Patient_ID", c(colors_pat, colors_many), rect_col = 'black', line_width = 0.1)

#Figure 1e - cluster plot
tsne_plot_no_outline(df, FILL = df$Cluster, fill_colors = colors_fig)

#Figure 1f - cluster pie chart
pie_chart()

#Figure 1h - cluster pie chart
genes <- c("CX3CR1__chr3", "HLA-DRA__chr6", "SPP1__chr4", "CCL2__chr17")

#gene expression tsne plots
plotexptsne2(.gene=genes[1], point_size = 4)
plotexptsne2(.gene=genes[2], point_size = 4)
plotexptsne2(.gene=genes[3], point_size = 4)
plotexptsne2(.gene=genes[4], point_size = 4)

#gene expression line plots
#for this plot one needs a long data plot
data_long <- make_data_long(.genes = genes)
#generate selected line plots
genes <- c("CX3CR1__chr3", "HLA.DRA__chr6", "SPP1__chr4", "CCL2__chr17")
gene_line_plot(.gene=genes[1])
gene_line_plot(.gene=genes[2])
gene_line_plot(.gene=genes[3])
gene_line_plot(.gene=genes[4])

#figure 3d - Regions t-SNE map
tsne_plot_no_outline(df[df$Region!="Mixed",], FILL=df[df$Region!="Mixed",]$Region, fill_colors = colors_many) +
  scale_fill_manual(values = toupper(c('#ef8a62', '#67a9cf'))) +
  scale_color_manual(values = toupper(c('#ef8a62', '#67a9cf')))

#figure 3e - Regions marimekko chart
df[df$Region!="Mixed",] %>%
  droplevels() %>%
  mosaicGG2(data=.,FILL = "Region", X="Cluster") +
  scale_fill_manual(values = toupper(c('#ef8a62', '#67a9cf')))

#statistical testing for regions enrichment in clusters
#to achieve best results for statistical testing, we conducted it for all clusters that are larger than 1% of all cells, including C4
retain_cl2 <- names(table(sc@cpart)[table(sc@cpart) > dim(sc@ndata)[2]/100])
df2 <- data.frame(ID=names(sc@cpart[sc@cpart %in% retain_cl2]), Cluster = sc@cpart[sc@cpart %in% retain_cl2], Region=ifelse(grepl('WM',names(sc@cpart[sc@cpart %in% retain_cl2])), 'WM', 
                                                                        ifelse(grepl('GM',names(sc@cpart[sc@cpart %in% retain_cl2])), 'GM', 'Mixed')))
clusters <- data.frame(table(sc@cpart))
colnames(clusters) <- c('Cluster', 'freq_clust')
regions <- as.data.frame(table(df2$Cluster, df2$Region))
colnames(regions) <- c("Cluster", "Region", "Freq_region")

regions_wide <- spread(regions, Region, Freq_region)

regions_df <- regions_wide %>%
  left_join(clusters)

regions_df$net_freq_clust <- regions_df$freq_clust-regions_df$Mixed

#hypergeometric test for grey matter
reg_clust_ordered_gm <- data.frame(q=regions_df$GM, 
                                   m=sum(regions_df$GM), 
                                   n=sum(regions_df$WM),
                                   k=regions_df$net_freq_clust)
reg_clust_ordered_gm %>%
  mutate(p_hyper = apply(., MARGIN = 1, function(x) 1-phyper(x[[1]]-1, x[[2]], x[[3]], x[[4]])), #probability to get q or more successes in populaton
         padj = p.adjust(p_hyper, method="BH"),
         Cluster = regions_df$Cluster,
         Significance = case_when(padj<0.05 & padj>0.01 ~ '*',
                                  padj<0.01 & padj>0.001 ~ '**',
                                  padj<0.001 ~ '***',
                                  TRUE ~ 'n.s.'))

#hypergeometric test for white matter
reg_clust_ordered_wm <- data.frame(q=regions_df$WM, 
                                   m=sum(regions_df$WM), 
                                   n=sum(regions_df$GM),
                                   k=regions_df$net_freq_clust)
reg_clust_ordered_wm %>%
  mutate(p_hyper = apply(., MARGIN = 1, function(x) 1-phyper(x[[1]]-1, x[[2]], x[[3]], x[[4]])), #probability to get q or more successes in populaton
         padj = p.adjust(p_hyper, method="BH"),
         Cluster = regions_df$Cluster,
         Significance = case_when(padj<0.05 & padj>0.01 ~ '*',
                                  padj<0.01 & padj>0.001 ~ '**',
                                  padj<0.001 ~ '***',
                                  TRUE ~ 'n.s.'))

#figure 3h
tsne_plot_no_outline(FILL=df$Age) +
  scale_fill_viridis('viridis', direction = -1)+
  scale_color_viridis('viridis', direction = -1)

#figure 3i
mosaicGG2(data=df,X="Cluster", FILL = "Age_bin", colors = viridis(3, direction = -1)) 

#statistical testing - ages
#again, I included all clusters that are >1% in size. To this end I updated the df2 data frame
df2 <- df2 %>%
  mutate(Patient_ID = gsub("_.*", "", .$ID)) %>%
  left_join(unique(df[, c("Patient_ID", "Age", "Age_bin")]))
clusters <- data.frame(table(sc@cpart))
colnames(clusters) <- c('Cluster', 'freq_clust')
ages <- as.data.frame(table(df2$Cluster, df2$Age_bin))
colnames(ages) <- c("Cluster", "Age_bin", "Freq_region")
levels(ages$Age_bin) <- c('young', 'high', 'middle')

ages_wide <- spread(ages, Age_bin, Freq_region)

ages_df <- ages_wide %>%
  left_join(clusters)

#young
age_clust_ordered_young <- data.frame(q=ages_df$young, 
                                      m=sum(ages_df$young), 
                                      n=sum(ages_df$middle, ages_df$high),
                                      k=ages_df$freq_clust)

age_clust_ordered_young %>%
  mutate(p_hyper = apply(., MARGIN = 1, function(x) 1-phyper(x[[1]]-1, x[[2]], x[[3]], x[[4]])), #probability to get q or more successes in populaton
         padj = p.adjust(p_hyper, method="BH"),
         Cluster = ages_df$Cluster,
         Significance = case_when(padj<0.05 & padj>0.01 ~ '*',
                                  padj<0.01 & padj>0.001 ~ '**',
                                  padj<0.001 ~ '***',
                                  TRUE ~ 'n.s.'))

#middle
age_clust_ordered_middle <- data.frame(q=ages_df$middle, 
                                       m=sum(ages_df$middle), 
                                       n=sum(ages_df$young, ages_df$high),
                                       k=ages_df$freq_clust)

age_clust_ordered_middle %>%
  mutate(p_hyper = apply(., MARGIN = 1, function(x) 1-phyper(x[[1]]-1, x[[2]], x[[3]], x[[4]])), #probability to get q or more successes in populaton
         padj = p.adjust(p_hyper, method="BH"),
         Cluster = ages_df$Cluster,
         Significance = case_when(padj<0.05 & padj>0.01 ~ '*',
                                  padj<0.01 & padj>0.001 ~ '**',
                                  padj<0.001 ~ '***',
                                  TRUE ~ 'n.s.'))

#high
age_clust_ordered_high <- data.frame(q=ages_df$high, 
                                     m=sum(ages_df$high), 
                                     n=sum(ages_df$young, ages_df$middle),
                                     k=ages_df$freq_clust)

age_clust_ordered_high %>%
  mutate(p_hyper = apply(., MARGIN = 1, function(x) 1-phyper(x[[1]]-1, x[[2]], x[[3]], x[[4]])), #probability to get q or more successes in populaton
         padj = p.adjust(p_hyper, method="BH"),
         Cluster = ages_df$Cluster,
         Significance = case_when(padj<0.05 & padj>0.01 ~ '*',
                                  padj<0.01 & padj>0.001 ~ '**',
                                  padj<0.001 ~ '***',
                                  TRUE ~ 'n.s.'))

# Figure 3j
spp1_df <- read.csv2(file.path("data", "counting_SPP1.csv"), stringsAsFactors = F)[,-c(1, 4:9)]

#tidy up the df
colnames(spp1_df)[c(1,4,5)] <- c('ID','percent_spp1','total_iba1')
spp1_df$percent_spp1 <- as.numeric(spp1_df$percent_spp1)
spp1_df$spp1_reactivity <- ifelse(grepl('-$', spp1_df$Region), 'neg','pos')
spp1_df$Region <- gsub(".$", '', spp1_df$Region)
spp1_df$Age_bin <- ifelse(spp1_df$Age < 30, '<30 years', 
                          ifelse(spp1_df$Age >50, '>50 years', '30-50 years'))

spp1_df$Age_bin <- factor(spp1_df$Age_bin, levels = c('<30 years', '30-50 years' , '>50 years'))

#plot

spp1_df_bar <- spp1_df[spp1_df$spp1_reactivity == 'pos',] %>% 
  group_by(Age_bin, Region) %>%
  summarise(mean_pos = mean(percent_spp1),
            median_pos = median(percent_spp1),
            sd_pos = sd(percent_spp1),
            se_pos = sd(percent_spp1)/sqrt(n()))

age_spp1_plot <- ggplot(spp1_df_bar, aes(Age_bin, mean_pos, fill = Age_bin)) +
  geom_col() +  
  geom_errorbar(aes(ymin = mean_pos - se_pos, ymax = mean_pos + se_pos), width=0, lwd=1)+
  geom_point(data=spp1_df[spp1_df$spp1_reactivity=='pos',], aes(Age_bin, percent_spp1, fill= Age_bin), pch=21, size=5, color='white') +
  scale_fill_viridis(direction = -1, discrete = T) +
  facet_wrap(~Region, ncol = 1) +
  theme_minimal() +
  theme(text = element_text(size=20),
        legend.key = element_blank())+
  coord_flip() +
  labs(y='% SPP1+ cells', x = 'Age bin')

age_spp1_plot 

#Extended data Figure 3b - cell distance heatmap
clustheatmap(sc,final=TRUE)

#Extended data Figure 4
#Figure 1h - cluster pie chart
genes <- c("TMEM119__chr12", "APOE__chr19", "CD74__chr5", "IFI44L__chr1","LPL__chr8")

#gene expression tsne plots
plotexptsne2(.gene=genes[1], point_size = 4)
plotexptsne2(.gene=genes[2], point_size = 4)
plotexptsne2(.gene=genes[3], point_size = 4)
plotexptsne2(.gene=genes[4], point_size = 4)
plotexptsne2(.gene=genes[5], point_size = 4)

#gene expression line plots
#for this plot one needs a long data plot
data_long <- make_data_long(.genes = genes) %>% na.omit
#generate selected line plots
gene_line_plot(.gene=genes[1])
gene_line_plot(.gene=genes[2])
gene_line_plot(.gene=genes[3])
gene_line_plot(.gene=genes[4])
gene_line_plot(.gene=genes[5])

#Extended data Figure 7b
data_t <- as.data.frame(t(as.matrix(sc@ndata)))
data_t$ID <- rownames(data_t)
data_t <- data_t %>% left_join(df)
data_t$Age_bin <- factor(data_t$Age_bin, levels= c("<30 years", "30-50 years", ">50 years"))

ggplot(na.omit(data_t[data_t$Region != 'Mixed',]), aes(Cluster, SPP1__chr4-0.1, fill = Age_bin)) +
  geom_violin(scale = 'width', lwd=0.25) +
  geom_boxplot(color='dark grey', width=0.5, position=position_dodge(0.9), outlier.shape = NA)+
  scale_fill_viridis(discrete = T, direction = -1) +
  theme_minimal() +
  labs(y='SPP1 expression') +
  facet_wrap(~ Region, nrow = 2)

library(MASS)
mod <- glm.nb(SPP1__chr4 ~ Cluster * Age_bin * Region, data = na.omit(data_t[data_t$Region != 'Mixed',]))
summary(mod)
mod_aov <- aov(mod)
broom::tidy(TukeyHSD(mod_aov))
