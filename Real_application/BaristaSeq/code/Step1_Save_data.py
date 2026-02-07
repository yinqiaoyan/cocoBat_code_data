##############################################################
#######################  BaristaSeq  #########################
##############################################################

########################## Note ##############################
## Choose the working directory as the directory where 
## this code file locates.
##############################################################

import scanpy as sc
import pandas as pd

FilePath = "../input_data/"
for i in [1, 2, 3]:
    adata = sc.read_h5ad(FilePath + f"Slice_{i}.h5ad")

    df = adata.to_df()
    df.to_csv(FilePath + f"gene_count_matrix_raw{i}.csv", sep=",")

    pd.DataFrame(adata.obsm["spatial"]).to_csv(FilePath + f"coordinates{i}.csv", index=False)

    adata.obs.to_csv(FilePath + f"ground_truth{i}.csv", index=True)
