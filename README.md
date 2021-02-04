### sankowski-et-al-microglia
This is the code repository for the analyses in the paper by Sankowski, Böttcher et al, Nat. Neurosci, 2019 (https://www.nature.com/articles/s41593-019-0532-y). It should be understood as a step by step instruction on how to reproduce the findings for the paper. Since it is written in R, users need to install the R programming language on their system. The comprehensive R Archive Network (CRAN) is a reliable source for the current R distribution: https://cran.r-project.org/

For more ease in using the code I recommend additionally downloading the Rstudio IDE here: https://rstudio.com/products/rstudio/download/. To ensure reproducibility of the findings I will start with the instructions on which packages to use and then walk through the steps that I had taken when working on this exciting project. Files are numbered in the recommended order that they should be run starting at 0 to set up the system. The R folder contains the code, the data folder all the data associated with this project. To ensure that the code runs across systems the file.path function is consistently used. The current version of my system with the current package versions can be found under "data/system-information.csv".

If you have any questions, please don't hesitate to raise an issue or contact me directly at roman.sankowski@uniklinik-freiburg.de.

In the following I will list some comments for the respective files.

####0_setup.R
This step will check for the packages available on your system and add the ones that are not present yet. Since, at least at the time of the preparation of the associated manuscript some of the packages were only available on bioconductor, I have set up the installation with BiocManager. This way packages available on CRAN and on bioconductor will be installed. For some packages, you will need to make sure that you have certain software installed on your system, e.g. Macs need to have Xcode for some packages. Unfortunately, I may not be competent enought to remotely help in every case when issues arise. Therefore, I ask you to address these problems with your local bioinformaticion/bioinformaticist (https://www.biostars.org/p/1183/).
What I can tell you is that the packages for deeplearning were by far the most difficult to install at the time of the analysis of these data. In these cases please refer to the documentation of the respective packages. When it worked for me I was quite happy and proud :)

####1_load_control_data.R
This is a somewhat convoluted script that I had used to load the data back in 2018. When re-running all analysis I realized that in order to more reliably reproduce the clustering of the paper you need not only the .coutt.csv, but also the other cout*.csv files. We have updated the GEO repository associated with this paper with these files. All counts files can be found in the data/counts folder. 

####2_RaceID_with_all_control_cells.R
At this point your computer will for the first time need quite some time to run RaceID on all control cells contained in this study. This script will return the sc object that you can use to obtain the equivalent of figure 1b containing different cell types. Please note that the layout of your t-SNE map will be different as the underlying algorithm has a random process that leads to different embeddings every time. You can transfer the embeddings from the original paper based on the provided table ("data/all_ctrl_cells_clusters+embeddings.csv"2).

####3_RaceID_with_microglia_control_cells.R
This is the code that was the basis for most of the analyses for figures 1-3. The used parameters are copied adopted from the RaceID manual as it was in 2018 when the analysis was performed. 
