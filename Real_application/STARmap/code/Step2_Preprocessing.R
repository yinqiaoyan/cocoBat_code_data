##############################################################
#######################    STARmap   #########################
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
library(scuttle)


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


for (i in 1:3) {
  geneData_raw = read.csv(paste0("../input_data/gene_count_matrix_raw", i, ".csv"))
  coord = read.csv(paste0("../input_data/coordinates", i, ".csv"))
  
  ## Normalization: Take logarithm
  cell_names = geneData_raw[, 1]
  geneData_raw_noName = geneData_raw[, -1]
  rownames(geneData_raw_noName) = cell_names
  sce <- SingleCellExperiment(assays=list(counts=t(geneData_raw_noName)),
                              colData = coord)
  sce <- logNormCounts(sce)
  write.csv(sce@assays@data@listData[["logcounts"]], file = paste0("../input_data/logcounts", i, ".csv"))
}




