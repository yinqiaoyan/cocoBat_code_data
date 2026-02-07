##############################################################
####################### Simulation 1 #########################
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
library(ggplot2)
library(umap)


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Load data
load("../input_data/simulated_data.RData")
meta_df <- read.csv("../input_data/Yall_meta_data.csv", header = TRUE, row.names = 1)

load("../result_data/cocoBat_res.RData")
ComBat_res   <- read.csv("../result_data/ComBat_res.csv", header = TRUE, row.names = 1)
Harmony_res  <- read.csv("../result_data/harmony_res.csv", header = TRUE, row.names = 1)
reComBat_res <- read.csv("../result_data/reComBat_res.csv", header = TRUE, row.names = 1)
batch = meta_df$batch
true_labels_vec = meta_df$celltype

correct_data = cocoBat_res[["corrected_data"]]


methods = c("ComBat", "reComBat", "Harmony")

for (method_name in methods) {
  if (method_name == "ComBat") { 
    ComBat_list <- split(seq_len(nrow(ComBat_res)), meta_df$batch)
    ComBat_list <- lapply(ComBat_list, function(idx) t(ComBat_res[idx, , drop = FALSE]))
  } else if (method_name == "reComBat") { 
    reComBat_list <- split(seq_len(nrow(reComBat_res)), meta_df$batch)
    reComBat_list <- lapply(reComBat_list, function(idx) t(reComBat_res[idx, , drop = FALSE]))
  } else if (method_name == "Harmony") { 
    Harmony_list <- split(seq_len(nrow(Harmony_res)), meta_df$batch)
    Harmony_list <- lapply(Harmony_list, function(idx) t(Harmony_res[idx, , drop = FALSE]))
  }
}


##################################
#####    Figure:heatmaps     #####
##################################
# heatmap
methods = c("rawData", "cocoBat", "ComBat", "reComBat", "Harmony")

for (method_name in methods) {
  if (method_name == "rawData") { 
    use_data <- Y_mat
  } else if (method_name == "cocoBat") { 
    use_data <- correct_data
  } else if (method_name == "ComBat") { 
    use_data <- ComBat_list
  } else if (method_name == "reComBat") { 
    use_data <- reComBat_list
  } else if (method_name == "Harmony") { 
    use_data <- Harmony_list
  }
  
  for (batch_num in 1:3) {
    tmpMatrix = use_data[[batch_num]]
    dim(tmpMatrix)
    
    cor_matrix <- cor(tmpMatrix)  
    cor_df <- as.data.frame(cor_matrix)
    cor_df$row <- rownames(cor_df)  
    cor_df <- tidyr::gather(cor_df, col, value, -row) 
    cor_df$row = factor(cor_df$row, levels=rev(unique(cor_df$row)))
    cor_df$col = factor(cor_df$col, levels=unique(cor_df$col))
    
    pheatmap <- ggplot(cor_df, aes(x = col, y = row, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(
        low = "#1E88E5", mid = "#FFFFFF", high = "#D82F2F",
        midpoint = 0,
        limits = c(-0.5, 1),          
        oob = scales::squish,         
        breaks = c(-0.5, 0, 0.5, 1),
        name = NULL
      )+
      labs(
        title = NULL,                         
        x = "Cells",       
        y = "Cells",      
        fill = NULL                          
      ) +
      scale_x_discrete(position = "top") +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32),
      
        legend.key.height = grid::unit(1.2, "cm"),
        legend.key.width  = grid::unit(0.6, "cm"),
        legend.text = element_text(size = 30),

        plot.title = element_blank()
      )
    ggsave(paste0("../figures/heatmap_", method_name, "_batch", batch_num, ".png"), pheatmap, width = 8, height = 7, dpi = 100)
  }
}



################################
#####    Figure:violin     #####
################################
# violin plots

methods = c("rawData", "cocoBat", "ComBat", "reComBat", "Harmony")

for (cell_type_id in 1:3) {
  for (cor_method in c("pearson", "spearman")) {
    for (method_name in methods) {
      
      if (method_name == "rawData") { 
        Y1_R1 = Y_mat[[1]][, true_labels_save[[1]] == cell_type_id]
        Y2_R1 = Y_mat[[2]][, true_labels_save[[2]] == cell_type_id]
        Y3_R1 = Y_mat[[3]][, true_labels_save[[3]] == cell_type_id]
      } else if (method_name == "cocoBat") { 
        Y1_R1 = correct_data[[1]][, true_labels_save[[1]] == cell_type_id]
        Y2_R1 = correct_data[[2]][, true_labels_save[[2]] == cell_type_id]
        Y3_R1 = correct_data[[3]][, true_labels_save[[3]] == cell_type_id]
      } else if (method_name == "ComBat") { 
        Y1_R1 = ComBat_list[[1]][, true_labels_save[[1]] == cell_type_id]
        Y2_R1 = ComBat_list[[2]][, true_labels_save[[2]] == cell_type_id]
        Y3_R1 = ComBat_list[[3]][, true_labels_save[[3]] == cell_type_id]
      } else if (method_name == "reComBat") { 
        Y1_R1 = reComBat_list[[1]][, true_labels_save[[1]] == cell_type_id]
        Y2_R1 = reComBat_list[[2]][, true_labels_save[[2]] == cell_type_id]
        Y3_R1 = reComBat_list[[3]][, true_labels_save[[3]] == cell_type_id]
      } else if (method_name == "Harmony") { 
        Y1_R1 = Harmony_list[[1]][, true_labels_save[[1]] == cell_type_id]
        Y2_R1 = Harmony_list[[2]][, true_labels_save[[2]] == cell_type_id]
        Y3_R1 = Harmony_list[[3]][, true_labels_save[[3]] == cell_type_id]
      }
      
      cor_values_list = list(batch1 = numeric(ncol(Y1_R1)),
                             batch2 = numeric(ncol(Y2_R1)),
                             batch3 = numeric(ncol(Y3_R1)))
      names(cor_values_list) = c("Batch 1", "Batch 2","Batch 3")
      
      unlist(Map(function(x){length(x)}, cor_values_list))
      
      tmpDataList = list(Y1_R1, Y2_R1, Y3_R1)
      
      for (i in 1:3) {
        cor_matrix <- cor(tmpDataList[[i]], method = cor_method)
        cor_values_list[[i]] <- cor_matrix[lower.tri(cor_matrix)]
      }
      
      df <- stack(cor_values_list)
      
      pp <- ggplot(df, aes(x = ind, y = values)) +
        geom_violin(scale = "width", width = 0.9, trim = TRUE, linewidth = 1.2) +
        coord_cartesian(ylim = c(-0.35, 1)) +
        scale_y_continuous(
          breaks = seq(-0.25, 1, by = 0.25),
          minor_breaks = NULL
        ) +
        labs(x = NULL, y = "Correlation coefficients") +
        theme_classic() +
        theme(
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          
          panel.grid.major.y = element_line(linewidth = 0.4, linetype = "dashed"),
          panel.grid.major.x = element_blank(),
          
          axis.text.x  = element_text(size = 38, margin = margin(t = 10)),
          axis.text.y  = element_text(size = 30),
          axis.title.y = element_text(size = 30)
        )
      
      ggsave(paste0("../figures/violin_", method_name, "_", cor_method, "_cellType", cell_type_id, ".png"), 
             pp, width = 8, height = 7, dpi = 100)
    }
  }
}





##############################
#####    Figure:UMAP     #####
##############################
# umap
plot_color=c("#A674DE", "#529D3B", "#F08633")

batch = meta_df$batch
batch <- gsub("Batch1", "Batch 1", batch, fixed = TRUE)
batch <- gsub("Batch2", "Batch 2", batch, fixed = TRUE)
batch <- gsub("Batch3", "Batch 3", batch, fixed = TRUE)
true_labels_vec = meta_df$celltype
methods = c("rawData", "cocoBat", "ComBat", "reComBat", "Harmony")

for (method_name in methods) {
  if (method_name == "rawData") { 
    use_data <- t(Reduce(cbind, Y_mat))
  } else if (method_name == "cocoBat") { 
    use_data <- t(Reduce(cbind, correct_data))
  } else if (method_name == "ComBat") { 
    use_data <- ComBat_res
  } else if (method_name == "reComBat") { 
    use_data <- reComBat_res
  } else if (method_name == "Harmony") { 
    use_data <- Harmony_res
  }
  
  data.umap = umap(use_data, n_components = 2, random_state = 30) 
  layout <- data.umap[["layout"]] 
  layout <- data.frame(layout)
  final <- cbind(layout, batch)
  dim(final)
  
  use_final = final
  xx = use_final[, 1]
  yy = use_final[, 2]
  ppdata = data.frame(x = xx, y = yy, c = use_final[, 3])
  pp <- ggplot(ppdata, aes(x = x, y = y)) +
    geom_point(aes(color = factor(c)), size = 0.2) +
    labs(x = "Dimension 1", y = "Dimension 2", color = NULL) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title.x = element_text(size = 30, margin = margin(t = 10)),
      axis.title.y = element_text(size = 30),
      
      axis.ticks = element_blank(),
      axis.text  = element_blank(),
      
      legend.title = element_text(size = 30),
      legend.text  = element_text(size = 30),
      legend.key.size = unit(2, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 10))) + 
    scale_color_manual(values=plot_color)
  
  ggsave(paste0("../figures/umap_", method_name, ".png"), pp, width = 8, height = 7, dpi = 100)
}




#################################
#####    Figure:scatter     #####
#################################
# scatter plots
plot_color = c("#A674DE", "#529D3B", "#F08633")

for (i in 1:3) {
  plot_dat <- meta_df[meta_df$batch == paste0("Batch", i), ]
  pp <- ggplot(data = plot_dat, aes(x=x, y=y)) +
    geom_point(aes(color = factor(celltype)), size = 12) +
    scale_y_reverse(breaks = 1:10) +
    scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50), limits = c(-1, 52)) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      
      axis.line = element_line(color = "black", linewidth = 2),
      
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text  = element_text(size = 20),
      axis.ticks = element_line(linewidth = 1.5),
      axis.ticks.length = unit(0.3, "cm"),
      
      legend.text  = element_text(size = 30),
      legend.title = element_text(size = 30),
      legend.key.size = unit(2, "cm")
    ) +
    scale_color_manual(values = plot_color, name = "Cell types")
  ggsave(paste0("../figures/cell_loc_batch", i, ".png"), pp, width = 7, height = 7, dpi = 100)
}







