##############################################################
#######################  BaristaSeq  #########################
##############################################################

########################## Note ##############################
## Choose the working directory as the directory where 
## this code file locates.
##############################################################

########################## Note ##############################
## Please first install required R packages.
## R: System_preparation.R
##############################################################


rm(list=ls())
library(harmony)


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())

## Load data
data_mat  <- read.csv("../input_data/Yall_data_mat.csv", header = TRUE, row.names = 1)
meta_data <- read.csv("../input_data/Yall_meta_data.csv", header = TRUE, row.names = 1)

harmony_embeddings <- harmony::RunHarmony(
  data_mat  = data_mat,      
  meta_data = meta_data,    
  vars_use  = "batch",
  do_pca    = FALSE          
)

write.csv(harmony_embeddings,  "../result_data/harmony_res.csv", row.names = TRUE)





