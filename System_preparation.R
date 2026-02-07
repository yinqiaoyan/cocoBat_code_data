##############################################################
################## System preparation ########################
##############################################################
rm(list=ls())

### Install CRAN packages, if necessary
if (!require("scuttle", quietly = TRUE))
  install.packages("scuttle", version="1.20.0")
if (!require("sva", quietly = TRUE))
  install.packages("sva", version="3.58.0")
if (!require("harmony", quietly = TRUE))
  install.packages("harmony", version="1.2.4")
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2", version="4.0.1")
if (!require("umap", quietly = TRUE))
  install.packages("umap", version="0.2.10.0")
if (!require("truncnorm", quietly = TRUE))
  install.packages("truncnorm", version="1.0-9")
if (!require("MCMCpack", quietly = TRUE))
  install.packages("MCMCpack", version="1.7-1")
if (!require("mvtnorm", quietly = TRUE))
  install.packages("mvtnorm", version="1.3-3")
if (!require("RcppArmadillo", quietly = TRUE))
  install.packages("RcppArmadillo", version="15.2.3-1")
if (!require("SingleCellExperiment", quietly = TRUE))
  install.packages("SingleCellExperiment", version="1.32.0")


## For Windows users, please update RTools to the latest version.



## Install R package cocoBat from a local file
## (Please make sure to install all required dependencies manually in advance)
install.packages("./cocoBat_1.0.tar.gz", repos = NULL, type="source")
## Or install from GitHub 
## (Automatically install dependencies)
devtools::install_github("yinqiaoyan/cocoBat", dependencies = TRUE)




