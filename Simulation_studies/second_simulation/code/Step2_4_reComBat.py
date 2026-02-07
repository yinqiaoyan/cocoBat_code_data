##############################################################
####################### Simulation 2 #########################
##############################################################

########################## Note ##############################
## Choose the working directory as the directory where 
## this code file locates.
##############################################################

import warnings
warnings.filterwarnings("ignore")
import pandas as pd
from reComBat import reComBat

## Load data
data_df = pd.read_csv("../input_data/Yall_data_mat.csv", index_col=0)   # cells × genes
meta_df = pd.read_csv("../input_data/Yall_meta_data.csv", index_col=0)  # cells × meta

assert (data_df.index.equals(meta_df.index)), "Cell order is inconsistent. Please check the CSV file."

batches = meta_df["batch"]
data_df = data_df.astype(float)

print("data_df shape:", data_df.shape)   
print("unique batches:", batches.unique())

## Run reComBat

combat = reComBat(
    model="elastic_net",   # elastic net regression
    parametric=True,       # parametric EB
    mean_only=False,       # adjust mean + variance (location + scale)
    optimize_params=True,  # optimize EB parameters automatically
    verbose=True
)

## fit + transform
corrected_df = combat.fit_transform(data_df, batches)
print("corrected_df shape:", corrected_df.shape) 

## Save results
corrected_df.to_csv("../result_data/reComBat_res.csv")
