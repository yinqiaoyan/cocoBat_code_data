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
library(scuttle)


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())

for (i in 1:3) {
  geneData_raw = read.csv(paste0("../input_data/gene_count_matrix_raw", i, ".csv"))
  coord = read.csv(paste0("../input_data/coordinates", i, ".csv"))
  truth_labels = read.csv(paste0("../input_data/ground_truth", i, ".csv"))
  
  ## Remove first column and add row names
  gene_names = geneData_raw[, 1]
  geneData_raw_noName = as.matrix(geneData_raw[, -1])
  rownames(geneData_raw_noName) = gene_names
  print(nrow(geneData_raw_noName))
  # 1525, 2042, 1690
  
  ## remove cells with NO gene expression
  rowSum_cells = rowSums(geneData_raw_noName)
  geneData_raw_noName = geneData_raw_noName[rowSum_cells > 0, ]
  print(nrow(geneData_raw_noName))
  # 1525, 2042, 1683 (remove 7 cells with NO gene expression)
  
  sce <- SingleCellExperiment(assays=list(counts=t(geneData_raw_noName)))
  sce <- logNormCounts(sce)
  
  
  write.csv(sce@assays@data@listData[["logcounts"]], file = paste0("../input_data/logcounts", i, ".csv"))
  
  if (sum(rowSum_cells == 0) > 0) {
    coord_new = coord[rowSum_cells > 0, ]
    truth_labels_new = truth_labels[rowSum_cells > 0, ]
    
    write.csv(coord_new, file = paste0("../input_data/coordinates", i, ".csv"), row.names = F)
    write.csv(truth_labels_new, file = paste0("../input_data/ground_truth", i, ".csv"), row.names = F)
  }
}




