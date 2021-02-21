## sankowski-et-al-microglia
This is the code repository for the analyses in the paper by Sankowski, Böttcher et al, Nat. Neurosci, 2019 (https://www.nature.com/articles/s41593-019-0532-y). It should be understood as a step by step instruction on how to reproduce the findings for the paper. Since it is written in R, users need to install the R programming language on their system. The comprehensive R Archive Network (CRAN) is a reliable source for the current R distribution: https://cran.r-project.org/

For more ease in using the code I recommend additionally downloading the Rstudio IDE here: https://rstudio.com/products/rstudio/download/. To ensure reproducibility of the findings I will start with the instructions on which packages to use and then walk through the steps that I had taken when working on this exciting project. Files are numbered in the recommended order that they should be run starting at 0 to set up the system. The R folder contains the code, the data folder all the data associated with this project. To ensure that the code runs across systems the file.path function is consistently used. The current version of my system with the current package versions can be found under "data/system-information.csv".

If you have any questions, please don't hesitate to raise an issue or contact me directly at roman.sankowski@uniklinik-freiburg.de.

In the following I will list some comments for the respective files in the bin/ folder.

To get started please load the counts data from here: 
https://drive.google.com/file/d/15k-vvl3nOvBh2OQeUCAfUOCsTNi-u3cx/view?usp=sharing
Put and unpack the folder inside the data folder and you should be good to go!

### 0_setup.R
This step will check for the packages available on your system and add the ones that are not present yet. Since, at least at the time of the preparation of the associated manuscript some of the packages were only available on bioconductor, I have set up the installation with BiocManager. This way packages available on CRAN and on bioconductor will be installed. For some packages, you will need to make sure that you have certain software installed on your system, e.g. Macs need to have Xcode for some packages. Unfortunately, I may not be competent enought to remotely help in every case when issues arise. Therefore, I ask you to address these problems with your local bioinformaticion/bioinformaticist (https://www.biostars.org/p/1183/).
What I can tell you is that the packages for deep learning were by far the most difficult to install at the time of the analysis of these data. In these cases please refer to the documentation of the respective packages. When it worked for me I was quite happy and proud :)

### 1_load_control_data.R
This is a somewhat convoluted script that I had used to load the data back in 2018. When re-running all analysis I realized that in order to more reliably reproduce the clustering of the paper you need not only the .coutt.csv, but also the other cout*.csv files that you will find in the "data/counts" folder. Another unexpected behavior during the loading process is the necessity for the counts files of Pat16. The counts of the other patients were initially loaded together with this case, and the counts table is different when this sample is not loaded concomitantly. Even though we have not included this sample before running RaceID, it is included for the loading. The lines causing this unexpected behavior are lines 148-155, where cells expressing more than 20% of the low-quality cell marker KCNQ1OT1 and genes with a correlation with KCNQ1OT1 < .65 were excluded. Since the sample was not included in the final manuscript, it is excluded from the following analyses. We have updated the GEO repository associated with this paper with these files. All counts files can be found here: https://drive.google.com/file/d/15k-vvl3nOvBh2OQeUCAfUOCsTNi-u3cx/view?usp=sharing. 

### 2_RaceID_with_all_control_cells.R
At this point your computer will for the first time need quite some time to run RaceID on all control cells contained in this study. This script will return the sc object that you can use to obtain the equivalent of figure 1b containing different cell types. Please note that the layout of your t-SNE map will be different as the underlying algorithm has a random process that leads to different embeddings every time. You can transfer the embeddings from the original paper based on the provided table ("data/all_ctrl_cells_clusters+embeddings_nn.csv"). The results of the RaceID analysis are saved under "data/sc_all_ctrl_cells.RData".

### 3_plotting_all_control_cells.R
This code contains the graphical output that was used for figures 1b and Extended data Figure 3a. Note that the layout of the t-SNE map may be different from the paper due to the random character of the embedding. If you want the same t-SNE layout of the function please use the provided table ("data/all_ctrl_cells_clusters+embeddings_nn.csv") as described in lines 7-12 of the script. To make sure that the clustering output of your RaceID run was consistent with the clustering in the original paper please run the code in line 12. If the clustering is consistent, you should only see numbers >0 on the diagonal and the other numbers should be 0. 
- Lines 14-16 contain the code for figure 1b
- Lines 18-37 contain the code for extended data figure 3a

### 4_RaceID_with_microglia_control_cells.R
This is the code that was the basis for most of the analyses for figures 1-3. The used parameters are copied and adopted from the RaceID manual as it was in 2018 when the analysis was performed. For your reference I added the manual in the data/ folder under "Reference_manual_RaceID3_StemID2.pdf". Genes associated with tissue dissociation were removed from RaceID clustering. Again, please note that since the t-SNE algorithm has a random component your t-SNE map will look different from the one that I have obtained after my analysis. The clustering, however is expected to be stable. One aspect that caused internal discussions was the output of the findoutliers function on the RaceID object. The purpose of this function is to find "rare" cell types by inducing some degree of overclustering, for a reference see the original publication by Grün et al, Nature 2015 (https://doi.org/10.1038/nature14966). Running this function gave rise to a somewhat bizzare cluster 14 that consisted of two distinct clouds of ~30 cells each that were positioned at opposite ends of the t-SNE map as outliers from clusters 2 and 9, respectively. Since, as shown in the paper, clusters 2 and 9 have substantial transcriptional and biological differences, we found the fusion of cells from both clusters counterintuitive. We had several options on how to proceed and decided to keep the outlier cells from clusters 2 and 9 with their respective parent clusters as they were assigned by the original clustering before outlier detection. An alternative would have been to further subcluster cluster 14 and separate both cell clouds. However, this would have led to the exclusion of these cells from analysis as the resulting clusters would have fallen under the 1% cutoff that we had set for the inclustion of clusters into the analysis. We wanted to avoid exclustion of cells as much as possible to maximize the insights gained from our data. 
If you wish to obtain plots and analyses with the embeddings and clustering of the original paper please enforce both from the file "data/ctrl_microglia_clusters+embeddings_nn.csv". This is done in lines 161-172 of the present script. In order to run the analysis on an identical RaceID object as was used for the paper, please run the following analyses with the object saved under "data/sc_ctrl_microglia_nn.RData". If you run the analyses with the object that you obtained from running this script yourself, your downstream results will be identical or consistent with the paper. 

### 5_diffgenes_ctrl.R
Differential gene expression testing was conducted similarly to the DESeq2 algorithm as described by Love et al, 2014 (https://doi.org/10.1186/s13059-014-0550-8). For this analysis comparisons between each cluster and all other pooled clusters and all other individual clusters were conducted. The output was saved in the folder "data/Cluster_specific_genes/Up". It is exported for analysis using lines 45-68. In order to ensure that genes associated with tissue dissociation overlay potential biological effects, a number of genes in question was removed from the list before plotting. The gene list obtained here was used to generate figure panels 1c, 1h, figure 2, figure 3h-j as well as Extended data figures 4-5 and 7d.

### 6_plotting_microglia_cells.R
This script generates the graphical output for figure panels 1c-f, 1h, figure 2, figure 3d-e, h-j as well as Extended data figures 3b, 4-5 and 7d. The code is annotated to ensure better readability. In lines 15-21 the cluster order is set in order for transcriptionally similar clusters to appear next to each other. This is achieved using the output of the clustheatmap function of RaceID. The graphical output of this function, a cell-to-cell distance heatmap, is shown in Extended data figure 3b of the manuscript. 
- Lines 23-40 exclude clusters containing fewer than 1% of all cells and filtering the cluster order accordingly. 
- Lines 42-74 generate a metadata table that is used for most abovementioned plots. 
- Lines 76-92 show the plot to generate a heatmap of the top 20 differentially expressed genes per cluster as shown in figure 1c. 
- Lines 94-96 show the marimekko patient plot from figure 1d. 
- Lines 98-99 show the code for figure 1e. 
- Lines 101-102 show the code for figure 1f. 
- Lines 104-121 show the code for figure 1h. 
- Lines 123-126 show the code for figure 3d.
- Lines 128-173 show the code for figure 3e and the statistical testing for regional enrichment.
- Lines 175-178 show the code for figure 3h.
- Lines 

### 7_neural_net_classifier_ctrl_microglia.R
This script shows the training and prediction using a deep neural network classifier for cell subset prediction. This classifier was developed in order to compare our cells to published human microglia datasets. Luckily me, keras was released at that time so I utilized it. The installation was a bit tricky, but doable with the documentation of the packages. I have added the required packages in the "0_setup.R" script. At the beginning of the present script you will also find commented out sections on the steps that I have taken to get this script up and running. Unfortunately, it still takes quite a few steps to get the packages up and running. To save th work for you, I have saved the outputs of the predictors under "data/deep-learning-predictor-controls"
The model was trained on the fdata matrix of the sc object. The training matrix was randomly sampled to make up 70% of the data. It was normalized, scaled and log2-transformed like shown before (https://www.biorxiv.org/content/biorxiv/early/2019/01/28/532093.full.pdf). The cluster assignments were used as the labels that the model was trained on. A sequential neural network was used. The architecture was as used by others before (https://www.biorxiv.org/content/biorxiv/early/2019/01/28/532093.full.pdf); namely it contained three layers with 100, 50 and 25 neurons with relU activation functions. The output layer contained 8 neurons with a softmax activation function. The 8 neurons corresponded to the 8 possible cluster labels. The resulting model is saved as "deep-learning-model-controls.h5". The predictions on each of the dataset were saved in the folder "data/deep-learning-predictor-controls".

### 8_plot_neural_net_classifier_ctrl_predictions.R
- Lines 6-27 contain the code to generate figure 1g.
- Lines 29-35 show the code to generate extended data figure 5b. Note that the displayed numbers are slightly different, as the data come from a re-run of the classifier.

### 9_GO_term_analysis_ctrl_microglia.R
This script leads to the data for figures 2 and Extended data figure 6.
- Lines 133-140 code for Figure 2b
- Lines 142-144 code for Figure 2c
- Lines 146-160 code for Extended data figure 6

### 10_load_gbm_data.R
This script loads the 4 glioblastoma samples and their age-matched controls. Since we saw age-dependent transcriptional changes, we included control cells from control individuals in a similar age range as the glioblastoma samples. To make the number of cells balanced, we compared 4 controls versus 4 glioblastoma. Like in the control samples, cells expressing more than 20% of the low-quality cell marker gene KCNQ1OT1 were excluded. Also genes with a correlation with KCNQ1OT1 < .65, ERCC and mitochondrial genes were excluded. 

### 11_RaceID_GBM_all_cells.R
This script runs RaceID on the glioblastoma cells and their age-matched controls without excluding non-microglia cells.

### 12_plotting_gbm_all_cells.R
This script generates the plots for Figure 9a  and Extended data figure 3a. 

### 13_RaceID_gbm_micr.R
This code runs RaceID on the glioblastoma microglia cells and their age-matched microglia controls. If you want to make sure that the clustering is consistent with the clustering from the published manuscript, please run lines 164-168. Furthermore, if you wish to run the downstream analysis and plots using the t-SNE layout of the published manuscript, please run the following lines 170-177. The output of these lines will be saven in the /data folder.

### 14_diffgenes_gbm.R
This code runs the differential gene expression analysis in a consistent way with script 5_diffgenes_ctrl.R.

### 15_plotting_gbm_microglia_cells.R
