library(Seurat)
library(stringr)
library(glmnet)
library(data.table)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

jinhong_folder <- './jinhong_deviance/'
setwd('./jinhong_deviance/data/')
source('./code_paper/v2high_dimensional_mean_testing_function.R')
source('./code_paper/supporting_function.R')
source('./code_paper/direct_ginv_function.R')

version <- 2
# cell_type <- 'T8'
if(version == 1){
cell_types <- c('T4', 'T8', 'cM', 'B', 'NK')
for(cell_type in cell_types){
  raw_gene <- data.table(read.csv(paste0('resid_', cell_type,'.csv')))
  raw_label <- data.table(read.csv(paste0('case_ctrl_', cell_type,'.csv')))
  raw_merged <- merge(raw_label, raw_gene, by = 'X')
  rownames(raw_merged) <- raw_merged$X
  raw_merged <- raw_merged[, -1]
  colnames(raw_merged)[1] <- 'case'
  # raw_merged$case <- sample(raw_merged$case)
  
  group1_data <- raw_merged[case == 1,]
  group1_data <- group1_data[, -1]
  group0_data <- raw_merged[case == 0,]
  group0_data <- group0_data[, -1]
  group0_data <- data.frame(group0_data)
  group1_data <- data.frame(group1_data)
  
  pooled_data <- rbind(group1_data, group0_data)
  pooled_data <- as.data.frame(apply(pooled_data,2,mynorm <- function(x){return (x/sd(x))}))
  
  group1_data <- pooled_data[1:nrow(group1_data),]
  group0_data <- pooled_data[(nrow(group1_data)+ 1):nrow(pooled_data),]
  
  # pooled_data <- as.matrix(pooled_data)
  # correlation_matrix <- nrow(pooled_data)^(-1) * t(pooled_data) %*% pooled_data
  # heatmap(abs(correlation_matrix)^1.1)
  lasso_result <-  anchored_lasso_testing(group1_data, group0_data, 
                                          pca_method = 'sparse_pca', n_folds = 10)
  simple_result <-  simple_pc_testing(group1_data, group0_data, 
                                          pca_method = 'sparse_pca', n_folds = 10, num_latent_factor = 3)
  debiased_result <-  debiased_pc_testing(group1_data, group0_data, 
                                      pca_method = 'sparse_pca', n_folds = 10, num_latent_factor = 3)
  
  result_to_save <- list(simple_pc_result = simple_result,
                         debiased_pc_result = debiased_result,
                         simple_lasso_test_result = lasso_result)
  
  result_file_name <- paste0(jinhong_folder,'result/', cell_type, '_cell_version_1.rds')
  saveRDS(result_to_save, file = result_file_name)
}
}else if(version == 2){
    cell_types <- c('T4', 'T8', 'cM', 'B', 'NK')
    for(cell_type in cell_types){
      raw_gene <- data.table(read.csv(paste0('resid_', cell_type,'.csv')))
      raw_label <- data.table(read.csv(paste0('case_ctrl_', cell_type,'.csv')))
      raw_merged <- merge(raw_label, raw_gene, by = 'X')
      rownames(raw_merged) <- raw_merged$X
      raw_merged <- raw_merged[, -1]
      colnames(raw_merged)[1] <- 'case'
      # raw_merged$case <- sample(raw_merged$case)
      
      group1_data <- raw_merged[case == 1,]
      group1_data <- group1_data[, -1]
      group0_data <- raw_merged[case == 0,]
      group0_data <- group0_data[, -1]
      group0_data <- data.frame(group0_data)
      group1_data <- data.frame(group1_data)
      
      pooled_data <- rbind(group1_data, group0_data)
      pooled_data <- as.data.frame(apply(pooled_data,2,mynorm <- function(x){return (x/sd(x))}))
      
      group1_data <- pooled_data[1:nrow(group1_data),]
      group0_data <- pooled_data[(nrow(group1_data)+ 1):nrow(pooled_data),]
      
      # pooled_data <- as.matrix(pooled_data)
      # correlation_matrix <- nrow(pooled_data)^(-1) * t(pooled_data) %*% pooled_data
      # heatmap(abs(correlation_matrix)^1.1)
      # lasso_result <-  anchored_lasso_testing(group1_data, group0_data, 
      #                                         pca_method = 'sparse_pca', n_folds = 10)
      simple_result <-  simple_pc_testing(group1_data, group0_data, 
                                          pca_method = 'sparse_pca', n_folds = 10, num_latent_factor = 10)
      debiased_result <-  debiased_pc_testing(group1_data, group0_data, 
                                              pca_method = 'sparse_pca', n_folds = 10, num_latent_factor = 10)
      
      result_to_save <- list(simple_pc_result = simple_result,
                             debiased_pc_result = debiased_result)
      
      result_file_name <- paste0(jinhong_folder,'result/', cell_type, '_cell_version_2.rds')
      saveRDS(result_to_save, file = result_file_name)
    }
}

