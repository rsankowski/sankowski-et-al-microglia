#Extended data figures
#Extended Data Figure 3 
#plot gene signatures
cell_signatures <- data.frame("monocytes"=c('CCR2', 'CLEC12A', 'PLAC8', 'FCN1', 'S100A9'),
                              "microglia"= c('P2RY12', 'CX3CR1', 'CSF1R', 'TMEM119', 'SLC2A5'),
                              "tcells"=c('TRAC', 'TRBC2', 'CD52', 'IL32', NA),
                              "myeloid_cells"=c('ITGAM',  'MS4A6A', 'TYROBP', 'CD14', NA),
                              "oligodendrocytes"=c('MBP',  'MOG', 'MAG', 'PLP1', NA))

for (i in colnames(cell_signatures)) {
  
  plt <- plotexptsne2(gene=na.omit(cell_signatures[[i]]) ,.sc=counts_all, point_size=3,line_width=0) + labs(title=i) 
  print(plt)
  
}                 
