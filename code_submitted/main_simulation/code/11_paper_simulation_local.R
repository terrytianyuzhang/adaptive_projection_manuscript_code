
work_directory <- './main_simulation/'

SIMULATION_BATCH <- 'InfAlternative_hard'

source(paste0(work_directory, 'code/v3high_dimensional_mean_testing_function.R'))
source(paste0(work_directory, 'code/supporting_function.R'))
source(paste0(work_directory, 'code/simulation_setting_parameter.R'))
num_features <- c(100)
num_repeat <- 10
print(paste0("running setting ", SIMULATION_BATCH))

library(PMA)
library(glmnet)
library(MASS)
library(data.table)
library(highmean)

stretch_one_repeat <- function(result_list, method = 'simple'){
  one_repeat_result <- matrix(repeat_index)
  one_repeat_result <- cbind(one_repeat_result, matrix(result_list$test_statistics, nrow = 1))
  one_repeat_result <- cbind(one_repeat_result, matrix(result_list$test_statistics_before_studentization, nrow = 1))
  one_repeat_result <- cbind(one_repeat_result, matrix(result_list$standard_error, nrow = 1))
  one_repeat_result <- data.frame(one_repeat_result)
  
  colnames(one_repeat_result) <- c("repeat_index",
                                   paste0('test_statistics_pc', 1:length(result_list$test_statistics)),
                                   paste0('parameter_estimate_pc', 1:length(result_list$test_statistics)),
                                   paste0('standard_error_pc', 1:length(result_list$test_statistics)))
  
  one_repeat_result <- data.frame(one_repeat_result)
  one_repeat_result$method <- method
  one_repeat_result$sample_size_1 <- sample_size_1
  one_repeat_result$sample_size_2 <- sample_size_2
  
  return(one_repeat_result)
}

#hard means hard-thresholding
set.seed(2018)

if(SIMULATION_BATCH %in% c('InfAlternative_hard', 'InfGlobal_hard', 'InfProject_hard')){
  all_result <- data.frame()
  result_file_name <- paste0(work_directory,
                             "result/",
                             "SETTING_", SIMULATION_BATCH,
                             "_my_methods",
                             ".rds" )
  
  for(sample_size_1 in candidate_sample_size_1){
    sample_size_2 <- sample_size_1
    for(number_feature in num_features){
      for(non_zero_number_feature in non_zero_number_features){
        
        true_mean_1 <- matrix(c(rep(non_zero_mean_1, non_zero_number_feature), 
                                rep(0, number_feature - non_zero_number_feature)), 
                              ncol = 1)
        
        true_mean_2 <- matrix(c(rep(non_zero_mean_2_part1, non_zero_number_feature), 
                                rep(non_zero_mean_2_part2, non_zero_number_feature), 
                                rep(0, number_feature - 2 * non_zero_number_feature)), 
                              ncol = 1)
        
        pc1 <- c(rep(1, non_zero_number_feature), 
                 rep(0, number_feature - non_zero_number_feature))
        pc1 <- pc1/norm(pc1, type = '2')
        
        pc2 <- c(rep(0, non_zero_number_feature), 
                 rep(1, non_zero_number_feature), 
                 rep(0, number_feature - 2 * non_zero_number_feature))
        pc2 <- pc2/norm(pc2, type = '2')
        
        
        simulation_covariance <- 100*pc1 %*% t(pc1) + 50*pc2 %*% t(pc2)
        simulation_covariance <- simulation_covariance + diag(noise_variance, nrow(simulation_covariance))
        
        for(repeat_index in 1:num_repeat){
          print(SIMULATION_BATCH)
          print(repeat_index)
          
          set.seed(repeat_index)
          # sample_1 <- data.frame(mvrnorm(sample_size_1,
          #                                mu = true_mean_1,
          #                                Sigma = simulation_covariance))
          # sample_2 <- data.frame(mvrnorm(sample_size_2,
          #                                mu = true_mean_2,
          #                                Sigma = simulation_covariance))
          sample_1 <- generate_zero_inflated_normal(sample_size_1,
                                                    true_mean_1, simulation_covariance, mask_probability)$sample
          sample_2 <- generate_zero_inflated_normal(sample_size_2,
                                                    true_mean_2, simulation_covariance, mask_probability)$sample
          
          
          
          simple_pc_result <- simple_pc_testing(sample_1, sample_2, 
                                                num_latent_factor = 2, 
                                                n_folds = 5, 
                                                pca_method = 'dense_pca')
          
          anchor_pc1_result <- anchored_lasso_testing(sample_1,
                                                       sample_2, 
                                                       pca_method = "hard",
                                                       num_latent_factor = 1)
          anchor_pc2_result <- anchored_lasso_testing(sample_1,
                                                     sample_2, 
                                                     pca_method = "hard",
                                                     num_latent_factor = 2)
          
          debiased_pc_result <- debiased_pc_testing(sample_1, 
                                                    sample_2, 
                                                    pca_method = "hard", 
                                                    num_latent_factor = 2)
          
          
          one_repeat_result <- stretch_one_repeat(simple_pc_result, 'simple')
          all_result <- combine_data_frame_of_diff_column(all_result, one_repeat_result)
          
          one_repeat_result <- stretch_one_repeat(debiased_pc_result, 'debiased')
          all_result <- combine_data_frame_of_diff_column(all_result, one_repeat_result)
          
          one_repeat_result <- stretch_one_repeat(anchor_pc1_result, 'lasso_pc1')
          all_result <- combine_data_frame_of_diff_column(all_result, one_repeat_result)
          
          one_repeat_result <- stretch_one_repeat(anchor_pc2_result, 'lasso_pc2')
          all_result <- combine_data_frame_of_diff_column(all_result, one_repeat_result)
          
          
          saveRDS(all_result, file = result_file_name)
        }
        
        # result_sub_folder <- paste0(work_directory, 
        #                             "data/", SIMULATION_BATCH, "/")
        # dir.create(result_sub_folder)
        # result_file_name <- paste0(result_sub_folder, 
        #                            "SETTING_", SIMULATION_BATCH,
        #                            "_my_methods",
        #                            ".rds" )
        
        
        
      }
    }
  }
}

