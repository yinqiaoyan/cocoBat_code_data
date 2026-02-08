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

## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## 3 batches
# Batch1: 1049 × 166
# Batch2: 1053 × 166
# Batch3: 1088 × 166
p = 3
K = 4
G = 166

Y_mat = vector("list", p)
coord_save = vector("list", p)
true_labels_save = vector("list", p)
num_samples = numeric(p)
dists_save = vector("list", p)
dists_tri_save = vector("list", p)

for (i in 1:3) {
  print(i)
  coord = read.csv(paste0("../input_data/coordinates", i, ".csv"))
  logcounts_data = read.csv(paste0("../input_data/logcounts", i, ".csv"))
  truth_labels = read.csv(paste0("../input_data/ground_truth", i, ".csv"))$ground_truth
  
  ### Preprocessing
  gene_names = logcounts_data[, 1]
  logcounts_data_noName = as.matrix(logcounts_data[, -1])
  rownames(logcounts_data_noName) = gene_names
  colnames(coord) = c("coord_x", "coord_y")
  
  num_samples[i] = ncol(logcounts_data_noName)
  
  ### Compute distance matrix
  dists = as.matrix(dist(coord, diag = T, upper = T, method = "manhattan"))
  dists_tri = dists[lower.tri(dists)]
  
  ### Save data
  dists_save[[i]] = dists
  dists_tri_save[[i]] = dists_tri
  Y_mat[[i]] = logcounts_data_noName
  coord_save[[i]] = coord
  true_labels_save[[i]] = truth_labels
}


## Save
save(Y_mat,
     coord_save,
     true_labels_save,
     file = "../input_data/listed_data.RData")

save(dists_tri_save,
     dists_save, 
     file = "../input_data/dist_mat_STARmap.RData")



## Save data for Harmony
Y_all <- do.call(cbind, Y_mat)    # dim: G × sum(num_samples)
rownames(Y_all) <- paste0("Gene", seq_len(G))

cell_names <- unlist(
  lapply(seq_len(p), function(b) {
    paste0("Batch", b, "_Cell", seq_len(num_samples[b]))
  })
)
colnames(Y_all) <- cell_names

# batch information
batch_vec <- rep(paste0("Batch", seq_len(p)), num_samples)
# cell type information
celltype_vec <- unlist(true_labels_save)
unique(celltype_vec)

meta_df <- data.frame(
  batch    = factor(batch_vec),
  celltype = factor(celltype_vec)
)
rownames(meta_df) <- colnames(Y_all)

# Save data
data_mat  <- t(Y_all)      # dim: sum(num_samples) × G
meta_data <- meta_df 

write.csv(data_mat,  "../input_data/Yall_data_mat.csv", row.names = TRUE)
write.csv(meta_data, "../input_data/Yall_meta_data.csv", row.names = TRUE)




