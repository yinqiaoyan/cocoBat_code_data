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
library(sva) # ComBat


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


load("../input_data/simulated_data.RData")

num_samples = sapply(true_labels_save, length)
edata = Reduce(cbind, Y_mat)
batch = c(rep(seq_along(num_samples), num_samples))

df_labels = data.frame(labels = unlist(true_labels_save))
modcombat = model.matrix(~-1 + as.factor(labels), data=df_labels)
modcombat = modcombat[, -1]

combat_edata = ComBat(dat=edata, batch=batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE,
                      mean.only = FALSE, ref.batch = NULL, BPPARAM = bpparam("SerialParam"))

write.csv(t(combat_edata),  "../result_data/ComBat_res.csv", row.names = TRUE)




