##############################################################
####################### Simulation 2 #########################
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
library(truncnorm)
library(MCMCpack) # rinvgamma
library(mvtnorm) # rmvnorm
library(RcppArmadillo)


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


Rcpp::sourceCpp("./utils_generation.cpp")


## Data generation
# batch number
p = 3
# spot number in each batch
ni <- 500 
# gene number
G <- 1000
# cell type number
K <- 3

nrows = 10
ncols = ni / nrows 
alpha_g_true = 2


set.seed(426)
### --- true mean --- ###
mu_true_val = c(0, 5, 10)  # dim: K
mu_true_mat = matrix(0, nrow = G, ncol = K)
# dim: G * K, with first column all zeros

for (k in 2:K) {
  mu_true_mat[, k] = rnorm(G, mean=mu_true_val[k], sd=0.1)
}

### --- additive batch effect --- ###
gamma_true_val = c(0, -2, 4)  # dim: p
gamma_true_mat = matrix(0, nrow = G, ncol = p)
# dim: G * p, with first column all zeros

for (i in 2:p) {
  gamma_true_mat[, i] = rnorm(G, mean=gamma_true_val[i], sd=0.1)
}


### --- scale batch effect --- ###
delta2_true_val = c(0.5, 1.5, 2.5)  # dim: p  
delta2_true_mat = matrix(NA, nrow = G, ncol = p)
# dim: G * p

for (i in 1:p) {
  delta2_true_mat[, i] = rnorm(G, mean=delta2_true_val[i], sd=0.1)
}

### --- corr batch effect --- ###
eta_true_val = c(1.5, 3.5, 5.5)  
eta_true_mat = matrix(NA, nrow = G, ncol = p)
# dim: G * p

for (i in 1:p) {
  eta_true_mat[, i] = rnorm(G, mean=eta_true_val[i], sd=0.1)
}


sigma2_g = 0.1

### --- coordinates --- ###
X_coord = rep(1:ncols, times = nrows)
Y_coord = rep(1:nrows, each = ncols)

true_labels_save = vector("list", p)
true_labels_save[[1]] = c(rep(1, 4*ncols), rep(2, 3*ncols), rep(3, 3*ncols))
true_labels_save[[2]] = c(rep(1, 3*ncols), rep(2, 4*ncols), rep(3, 3*ncols))
true_labels_save[[3]] = c(rep(1, 3*ncols), rep(2, 3*ncols), rep(3, 4*ncols))

coords = cbind(X_coord, Y_coord)
dists = as.matrix(dist(coords, diag = T, upper = T, method = "euclidean"))
dists_tri = dists[lower.tri(dists)]

## Generate data
Y_mat = GenerateYmat(p, G, ni, dists, mu_true_mat, gamma_true_mat, delta2_true_mat,
                     eta_true_mat, true_labels_save, alpha_g_true, sigma2_g)

coord_save = list(coords, coords, coords)
dists_save = list(dists, dists, dists)
dists_tri_save = list(dists_tri, dists_tri, dists_tri)

## Save
save(Y_mat, coord_save, true_labels_save,
     mu_true_mat, gamma_true_mat, delta2_true_mat, eta_true_mat,
     file = "../input_data/simulated_data.RData")


save(dists_tri_save,
     dists_save, file = "../input_data/dist_mat_sim2.RData")



## Save data for Harmony
Y_all = Reduce(cbind, Y_mat)
rownames(Y_all) <- paste0("Gene", seq_len(G))
cell_names <- unlist(
  lapply(seq_len(p), function(b) {
    paste0("Batch", b, "_Cell", seq_len(ni))
  })
)
colnames(Y_all) <- cell_names

# batch information
batch_vec <- rep(paste0("Batch", seq_len(p)), each = ni)
# cell type information
celltype_vec <- unlist(true_labels_save) 

meta_df <- data.frame(
  batch    = factor(batch_vec),
  celltype = factor(celltype_vec)
)
rownames(meta_df) <- colnames(Y_all)

coords_mat <- cbind(
  x = rep(X_coord, times = p),
  y = rep(Y_coord, times = p)
)
meta_df$x <- coords_mat[, "x"]
meta_df$y <- coords_mat[, "y"]

# Save data
data_mat  <- t(Y_all)      # dim: (p*ni) × G
meta_data <- meta_df 

write.csv(data_mat,  "../input_data/Yall_data_mat.csv", row.names = TRUE)
write.csv(meta_data, "../input_data/Yall_meta_data.csv", row.names = TRUE)




