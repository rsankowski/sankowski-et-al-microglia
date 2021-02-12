load(file.path("data", "sc_ctrl_microglia_nn.RData"))

mainDir <- getwd()
subDir <- file.path("data", "Cluster_specific_genes")
subsubDir_up <- 'Up'

if (file.exists(subDir)){
  setwd(file.path(mainDir, subDir))
} else {
  dir.create(file.path(mainDir, subDir))
  dir.create(file.path(mainDir, subDir, subsubDir_up))
  
  
}
for (i in c(1: max(sc@cpart))) {
  cl <- names(sc@cpart[sc@cpart %in% c(i)])
  rest <- names(sc@cpart[sc@cpart %in% sc@cpart[sc@cpart != i]])
  diffexp <- diffexpnb(sc@ndata,rest,cl,norm=F,vfit=sc@background$vfit) 
  diffexpgenes <- diffexp[["res"]]
  diffexpgenes <- subset(diffexpgenes, diffexpgenes$pval < 0.05)
  diffexpgenes <- subset(diffexpgenes, abs(diffexpgenes$log2FoldChange) > 0)
  diffexpgenes$GENEID <- rownames(diffexpgenes)
  diffexpgenes <- diffexpgenes %>% dplyr::arrange(padj)
  diffexpgenes$Cluster <- rep(i, nrow(diffexpgenes))
  diffgene_up <- diffexpgenes[diffexpgenes$log2FoldChange > 0, ]
  write.csv(diffgene_up, file = file.path(mainDir, subDir, "Up", paste0("diffgenes_cl", as.character(i), "_up_rest.csv")))
  
  for (j in c(1: max(sc@cpart))) {
    cl <- names(sc@cpart[sc@cpart %in% i])
    rest <- names(sc@cpart[sc@cpart %in% j])
    diffexp <- diffexpnb(sc@ndata,rest,cl,norm=F,vfit=sc@background$vfit) 
    diffexpgenes <- diffexp[["res"]]
    diffexpgenes <- subset(diffexpgenes, diffexpgenes$pval < 0.05)
    diffexpgenes <- subset(diffexpgenes, abs(diffexpgenes$log2FoldChange) > 0)
    diffexpgenes$GENEID <- rownames(diffexpgenes)
    diffexpgenes <- diffexpgenes %>% dplyr::arrange(padj)
    diffexpgenes$Cluster <- rep(i, nrow(diffexpgenes))
    diffgene_up <- diffexpgenes[diffexpgenes$log2FoldChange > 0, ]
    write.csv(diffgene_up, file = file.path(mainDir, subDir, "Up", paste0("diffgenes_cl", as.character(i), "_up_vs_cl", as.character(j), ".csv")))
    
  }
  
}

#export diffgenes
load_data <- function(path) { 
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)
  tables <- lapply(files, read.csv)
  do.call(rbind, tables)
}

up_genes <- load_data(file.path('data','Cluster_specific_genes',"Up"))
up_genes$GENEID <- as.character(up_genes$GENEID)

#remove not annotated or dissociation associated genes
up_genes2 <- up_genes[!grepl("^(HSPA1A|MTRNR2L8|MTRNR2L12|HSP90AA1|MALAT1|ZFP36L1|ZFP36|FOS|MALAT1|HSPB*|DUSP1|HSPH1|HSPA*|JUN|HSP90B1|RPS16|DNAJB1|H3F3B|HERPUD1|NEAT1|IVNS1ABP|HIST1H2BG|RP*|XIST)", up_genes$GENEID),  ]
up_genes2$GENEID <- gsub('_.*', '', up_genes2$GENEID)

up_genes %>% 
  filter(padj<0.05, log2FoldChange > 1) %>% 
  dplyr::arrange(Cluster) %>%
  write.csv(file.path('data', 'Human_ctrl_microglia_up_genes_padj<05_logfc>1.csv'))

up_genes2 %>% 
  filter(padj<0.05, up_genes2$log2FoldChange > 1) %>% 
  dplyr::arrange(Cluster) %>%
  dplyr::distinct(Cluster, GENEID) %>%
  write.csv(file.path('data','Human_ctrl_micr_unique_up_genes_wo_stress-genes_log2fc>1.csv'))

up_genes2 %>% 
  filter(padj<0.05) %>% 
  dplyr::arrange(Cluster) %>%
  dplyr::distinct(Cluster, GENEID) %>%
  write.csv(file.path('data','Human_ctrl_micr_unique_up_genes_wo_stress-genes_padj<05.csv'))

