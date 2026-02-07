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
library(cocoBat)


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


load("../input_data/listed_data.RData")

cocoBat_res = run_eb_mcmc(Y_mat,
                          true_labels_save,
                          dists_tri_save,
                          dists_save,
                          random_seed = 30,
                          maxit  = 100,
                          num_iters = 100,
                          n_threads = 8,
                          verbose   = TRUE)


save(cocoBat_res, file = "../result_data/cocoBat_res.RData")





