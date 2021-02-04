### sankowski-et-al-microglia
This is the code repository for the analyses in the paper https://www.nature.com/articles/s41593-019-0532-y.
It should be understood as a step by step instruction on how to reproduce the findings for the paper. Since it is written in R, users need 
to install the R programming language on their system. The comprehensive R Archive Network is a reliable source for the current R distribution:
https://cran.r-project.org/

For more ease in using the code I recommend additionally downloading the Rstudio IDE here: https://rstudio.com/products/rstudio/download/. 
To ensure reproducibility of the findings I will start with the instructions on which packages to use and then walk through the steps that I had 
taken when working on this exciting project. Files are numbered in the recommended order that they should be run starting at 0 to set up the system.
The bin folder contains the code, the data folder all the data associated with this project. To ensure that the code runs across systems the file.path 
function is consistently used. The current version of my system with the current package versions can be found under "data/system-information.csv".

If you have any questions, please don't hesitate to raise an issue or contact me directly at roman.sankowski@uniklinik-freiburg.de.


