library(RSpectra)

direct_ginv_testing <- function(sample_1,
                                sample_2 = NULL,
                                pca_method = "sparse_pca",
                                mean_method = "naive",
                                num_latent_factor = 1,
                                n_folds = 5){
  ##sample_2 is the control sample and the covariance matrix estimation is based on control only, this is verified
  
  ###split the samples into n_folds splits
  set.seed(1)
  sample_1 <- as.data.frame(sample_1)
  split_data <- vector("list", n_folds)
  sample_1_split_index <- index_spliter(1:nrow(sample_1),
                                        n_folds = n_folds)
  
  if(!is.null(sample_2)){
    sample_2 <- as.data.frame(sample_2)
    sample_2_split_index <- index_spliter(1:nrow(sample_2),
                                          n_folds = n_folds)
  }
  
  hyperparameter_shared_between_folds <<- -1
  ###process each split
  for(split_index in 1:n_folds){
    print(paste0("processiong fold ", split_index))
    
    sample_1_cross <- sample_1[sample_1_split_index[[split_index]], ]
    sample_1_nuisance <- sample_1[-sample_1_split_index[[split_index]], ]
    split_data[[split_index]]$sample_1_cross_index <- sample_1_split_index[[split_index]]
    
    if(!is.null(sample_2)){
      split_data[[split_index]]$sample_2_cross_index <- sample_2_split_index[[split_index]]
      sample_2_cross <- sample_2[sample_2_split_index[[split_index]], ]
      sample_2_nuisance <- sample_2[-sample_2_split_index[[split_index]], ]
    }else{
      sample_2_cross <- sample_2_nuisance <- NULL
    }
    
    
    split_data[[split_index]]$nuisance_collection <- direct_estimate_nuisance(sample_1_nuisance, 
                                                                          sample_2_nuisance, 
                                                                          mean_method = mean_method,
                                                                          num_latent_factor = num_latent_factor)
    split_data[[split_index]]$influence_function_value <- direct_evaluate_influence_function(sample_1_cross,
                                                                                             sample_2_cross,
                                                                                             nuisance_collection = split_data[[split_index]]$nuisance_collection)
    # testing data, should remove later
    
    # to_look_new <- split_data[[split_index]]$influence_function_value
    # to_look_new_multi <- split_data[[split_index]]$influence_function_value
    # to_look_old <- evaluate_influence_function(sample_1_cross, sample_2_cross, nuisance_collection = split_data[[split_index]]$nuisance_collection)
    # summary(to_look_new$for_variance_subject_1 - to_look_old$for_variance_subject_1)
    # to_look_new$influence_eigenvector_each_subject_1[1:10]
    # to_look_new_multi$influence_eigenvector_each_subject_1[1:10, 1]
    ##############
    
    if(!is.null(sample_2)){
      
      inner_product_1 <- split_data[[split_index]]$influence_function_value$inner_product_1
      inner_product_2 <- split_data[[split_index]]$influence_function_value$inner_product_2
      split_data[[split_index]]$test_statistic <- apply(X = inner_product_2, MARGIN = 2, FUN = mean) - apply(X = inner_product_1, MARGIN = 2, FUN = mean)
      #debias
      # sample_1_correction <- split_data[[split_index]]$influence_function_value$influence_eigenvector_each_subject_1
      sample_2_correction <- split_data[[split_index]]$influence_function_value$influence_eigenvector_each_subject_2
      
      split_data[[split_index]]$test_statistic <- split_data[[split_index]]$test_statistic + apply(X = sample_2_correction, MARGIN = 2, FUN = mean)
      
      ##variance
      split_data[[split_index]]$variance_sample_1 <- apply(X = split_data[[split_index]]$influence_function_value$for_variance_subject_1, MARGIN = 2, FUN = var)
      split_data[[split_index]]$variance_sample_2 <- apply(X = split_data[[split_index]]$influence_function_value$for_variance_subject_2, MARGIN = 2, FUN = var)
      
    }else{
      sample_1_individual <- split_data[[split_index]]$influence_function_value$influence_each_subject_1
      split_data[[split_index]]$variance_sample_1 <- apply(X = sample_1_individual, MARGIN = 2, FUN = var)
      split_data[[split_index]]$test_statistic <- apply(X = sample_1_individual, MARGIN = 2, FUN = mean)
    }
    
  }
  
  
  ####now combine the folds
  combine_1_factor_over_splits <- function(latent_factor_index){
    test_statistics_before_studentization <- 0
    variance_sample_1 <- variance_sample_2 <- 0
    
    for(split_index in 1:n_folds){
      inner_product_projection_direction <- crossprod(split_data[[1]]$nuisance_collection$estimate_leading_pc[,latent_factor_index], 
                                                      split_data[[split_index]]$nuisance_collection$estimate_leading_pc[,latent_factor_index])
      
      same_sign <- sign(inner_product_projection_direction)
      print(inner_product_projection_direction)
      if(same_sign == 0){
        print('the projection directions are orthogonal')
        same_sign <- 1
      }
      
      test_statistics_before_studentization <- test_statistics_before_studentization + same_sign * split_data[[split_index]]$test_statistic[latent_factor_index]
      
      variance_sample_1 <- variance_sample_1 + split_data[[split_index]]$variance_sample_1[latent_factor_index]
      
      if(!is.null(sample_2)){
        variance_sample_2 <- variance_sample_2 + split_data[[split_index]]$variance_sample_2[latent_factor_index]
      }
    }
    
    test_statistics_before_studentization <- test_statistics_before_studentization/n_folds
    variance_sample_1 <- variance_sample_1/n_folds
    
    if(!is.null(sample_2)) variance_sample_2 <- variance_sample_2/n_folds
    
    if(!is.null(sample_2)){
      standard_error <- sqrt((nrow(sample_1))^(-1) * variance_sample_1 + 
                               (nrow(sample_2))^(-1) *variance_sample_2)
    }else{
      standard_error <- sqrt((nrow(sample_1))^(-1) * variance_sample_1)
    }
    
    test_statistics <- test_statistics_before_studentization/standard_error
    
    tuple_one_latent_factor <- c(test_statistics_before_studentization, standard_error, test_statistics)
    return(tuple_one_latent_factor)
  }
  
  test_statistics_before_studentization <- standard_error <- test_statistics <- rep(NULL,num_latent_factor)
  for(latent_factor_index in 1:num_latent_factor){
    tuple_one_latent_factor <- combine_1_factor_over_splits(latent_factor_index)
    
    test_statistics_before_studentization[latent_factor_index] <- tuple_one_latent_factor[1]
    standard_error[latent_factor_index] <- tuple_one_latent_factor[2]
    test_statistics[latent_factor_index] <- tuple_one_latent_factor[3]
  }
  
  
  return(list(test_statistics = test_statistics,
              standard_error = standard_error,
              test_statistics_before_studentization = test_statistics_before_studentization,
              split_data = split_data))
}


direct_estimate_nuisance <- function(nuisance_sample_1,
                                     nuisance_sample_2 = NULL,
                                     pca_method = "sparse_pca",
                                     mean_method = "hard",
                                     num_latent_factor = 1){
  
  nuisance_sample_1 <<- nuisance_sample_1
  nuisance_sample_2 <<- nuisance_sample_2
  
  feature_number <- ncol(nuisance_sample_1)
  
  nuisance_sample_1 <- as.data.frame(nuisance_sample_1)
  nuisance_sample_size_1 <- nrow(nuisance_sample_1)
  
  ###estimate mean vector 
  estimate_mean_1 <- colMeans(nuisance_sample_1)
  # estimate_mean_1[abs(estimate_mean_1) < mean_threshold_1] <- 0
  if(mean_method == 'hard'){
    variance_feature_wise_centered <- apply(nuisance_sample_1, 2, var)
    mean_threshold <- 2*sqrt(variance_feature_wise_centered)*sqrt(2*log(2*feature_number*nuisance_sample_size_1)/nuisance_sample_size_1)
    estimate_mean_1[abs(estimate_mean_1) < mean_threshold] <- 0
  }

  centered_sample_1 <- sweep(nuisance_sample_1, 2, estimate_mean_1)
  
  #######
  
  if(!is.null(nuisance_sample_2)){
    
    nuisance_sample_2 <- as.data.frame(nuisance_sample_2)
    nuisance_sample_size_2 <- nrow(nuisance_sample_2)
    
    ###estimate mean vector
    estimate_mean_2 <- colMeans(nuisance_sample_2)
    # estimate_mean_2[abs(estimate_mean_2) < mean_threshold_2] <- 0
    
    if(mean_method == 'hard'){
      variance_feature_wise_centered <- apply(nuisance_sample_2, 2, var)
      mean_threshold <- 2*sqrt(variance_feature_wise_centered)*sqrt(2*log(2*feature_number*nuisance_sample_size_2)/nuisance_sample_size_2)
      estimate_mean_2[abs(estimate_mean_2) < mean_threshold] <- 0
    }
    
    centered_sample_2 <- sweep(nuisance_sample_2, 2, estimate_mean_2)
    
    ########
    sample_centered <- as.matrix(rbind(centered_sample_1, centered_sample_2))
    sample_centered <- as.matrix(centered_sample_2)
    
  }else{
    estimate_mean_2 <- NULL
    sample_centered <- as.matrix(centered_sample_1)
  }
  
  sample_cov <- (nrow(sample_centered))^(-1) * (t(sample_centered) %*% sample_centered)
  sample_cov_shrink <- sample_cov
  # sample_cov_shrink[abs(sample_cov_shrink) < 0.3 * max(abs(diag(sample_cov))) ] <- 0
  sample_cov_shrink[abs(sample_cov_shrink) < sqrt(log(nrow(sample_cov_shrink))/nrow(sample_centered))] <- 0
  diag(sample_cov_shrink) <- diag(sample_cov)
  # shrink_eigen <- eigen(sample_cov_shrink)
  
  
  if(pca_method == "sparse_pca"){
    ###estimate sparse principle component
    
    if(hyperparameter_shared_between_folds < 0){
      # cv_result <- SPC.cv(sample_centered,
      #                     sumabsvs = seq(1, sqrt(floor(feature_number)), length = 10))
      # hyperparameter_shared_between_folds <<- cv_result$bestsumabsv
      cv_result <- SPC.cv(sample_centered)
      # hyperparameter_shared_between_folds <<- max(1, 0.7*cv_result$bestsumabsv) ##I USED THIS FOR LUPAS REAL-DATA ANALYSIS
      hyperparameter_shared_between_folds <<- cv_result$bestsumabsv
    }
    print(hyperparameter_shared_between_folds)
    estimate_leading_pc <- SPC(sample_centered,
                               K = num_latent_factor,
                               sumabsv = hyperparameter_shared_between_folds)$v
    estimate_leading_pc <- estimate_leading_pc/apply(X = estimate_leading_pc, MARGIN = 2, FUN = norm, type = '2')
    
    colnames(estimate_leading_pc) <- paste0('pc', 1:ncol(estimate_leading_pc))
    rownames(estimate_leading_pc) <- colnames(nuisance_sample_2)
    
    ##estimate eigenvalue
    data_times_pc <- sample_centered %*% estimate_leading_pc
    estimate_eigenvalue <- (apply(X = data_times_pc, MARGIN = 2, 
                                  FUN = norm, type = '2'))^2/nrow(sample_centered)
    # estimate_leading_pc <- estimate_leading_pc/norm(estimate_leading_pc, type = '2')
  }else if(pca_method == "dense_pca"){
    estimate_leading_pc <- irlba(sample_centered, nv = num_latent_factor)$v
    colnames(estimate_leading_pc) <- paste0('pc', 1:ncol(estimate_leading_pc))
    rownames(estimate_leading_pc) <- colnames(nuisance_sample_2)
    
    ##estimate eigenvalue
    data_times_pc <- sample_centered %*% estimate_leading_pc
    estimate_eigenvalue <- (apply(X = data_times_pc, MARGIN = 2, 
                                  FUN = norm, type = '2'))^2/nrow(sample_centered)
  }else if(pca_method == 'eigen_covariance'){
    shrink_eigen <- eigs_sym(sample_cov_shrink, num_latent_factor)
    
    estimate_leading_pc <- shrink_eigen$vectors[,1:num_latent_factor]
    estimate_leading_pc <- matrix(estimate_leading_pc, ncol = num_latent_factor)
    rownames(estimate_leading_pc) <- colnames(nuisance_sample_2)
    
    
    estimate_leading_pc[abs(estimate_leading_pc) < min(0.02, sqrt(log(nrow(sample_cov_shrink))/nrow(sample_centered)))] <- 0
    for(latent_index in 1:num_latent_factor){
      if(norm(estimate_leading_pc[,latent_index], type = '2') == 0){
        print('opps')
        estimate_leading_pc[,latent_index] <- shrink_eigen$vectors[,latent_index]
      }else{
        estimate_leading_pc[,latent_index] <- estimate_leading_pc[,latent_index] / norm(estimate_leading_pc[,latent_index], type = '2')
      }
      
    }
    
    estimate_eigenvalue <- shrink_eigen$values[1:num_latent_factor]
    
  }
  
  return(list(estimate_leading_pc = estimate_leading_pc,
              estimate_mean_1 = estimate_mean_1,
              estimate_mean_2 = estimate_mean_2,
              estimate_eigenvalue = estimate_eigenvalue,
              estimate_cov = sample_cov_shrink))
}

direct_evaluate_influence_function <- function(cross_fitting_sample_1,
                                               cross_fitting_sample_2 = NULL,
                                               nuisance_collection){
  cross_fitting_sample_1 <<- cross_fitting_sample_1
  cross_fitting_sample_2 <<- cross_fitting_sample_2
  nuisance_collection <<- nuisance_collection
  
  
    
  cross_fitting_sample_1 <- as.data.frame(cross_fitting_sample_1)
  cross_fitting_sample_2 <- as.data.frame(cross_fitting_sample_2)
  estimate_cov <- nuisance_collection$estimate_cov
  estimate_leading_pc_all <- nuisance_collection$estimate_leading_pc
  estimate_eigenvalue_all <-  nuisance_collection$estimate_eigenvalue
  
  ###process each eigenvector
  influence_function_all <- inner_product_1_all <- inner_product_2_all <- for_variance_subject_1_all <- for_variance_subject_2_all <- NULL
  for(latent_index in 1:ncol(estimate_leading_pc_all)){
    estimate_leading_pc <- estimate_leading_pc_all[,latent_index]
    estimate_eigenvalue <- estimate_eigenvalue_all[latent_index]
    
    #####inner product statistics, also used in the simple statistics
    estimate_mean_1 <- nuisance_collection$estimate_mean_1
    inner_product_1 <- as.matrix(cross_fitting_sample_1) %*% estimate_leading_pc
    
    estimate_mean_2 <- nuisance_collection$estimate_mean_2
    inner_product_2 <- as.matrix(cross_fitting_sample_2) %*% estimate_leading_pc
    
    mean_diff <- estimate_mean_1 - estimate_mean_2
    
    ####influence function of the eigenvector
    estimate_mean_2 <- matrix(estimate_mean_2, ncol = 1)
    sudo_inv_part <- ginv(estimate_eigenvalue * diag(1, ncol(cross_fitting_sample_2)) - estimate_cov)
    influence_function <- NULL
    for(sample_2_index in 1:nrow(cross_fitting_sample_2)){
      one_data <- matrix(as.matrix(cross_fitting_sample_2)[sample_2_index,], ncol = 1)
      covariance_deviation <- (one_data - estimate_mean_2) %*% t(one_data - estimate_mean_2) - estimate_cov
      covariance_times_eigen_vector <- covariance_deviation %*% estimate_leading_pc
      before_time_mean_diff <- sudo_inv_part %*% covariance_times_eigen_vector
      # plot(before_time_mean_diff)
      one_number_one_sample <- crossprod(mean_diff, before_time_mean_diff)
      influence_function <- c(influence_function, one_number_one_sample)
    }
    
    for_variance_subject_1 <- inner_product_1 
    for_variance_subject_2 <- inner_product_2 + influence_function
    
    inner_product_1_all <- cbind(inner_product_1_all, inner_product_1)
    inner_product_2_all <- cbind(inner_product_2_all, inner_product_2)
    for_variance_subject_1_all <- cbind(for_variance_subject_1_all, for_variance_subject_1)
    for_variance_subject_2_all <- cbind(for_variance_subject_2_all, for_variance_subject_2)
    influence_function_all <- cbind(influence_function_all, influence_function)
  }
  
  
  
  # centered_sample_1 <- as.matrix(sweep(cross_fitting_sample_1, 2, estimate_mean_1)) 
  # centered_1_times_pc <- inner_product_centered_1 <- centered_sample_1 %*% estimate_leading_pc
  # 
  # centered_sample_2 <- as.matrix(sweep(cross_fitting_sample_2, 2, estimate_mean_2))
  # centered_2_times_pc <- inner_product_centered_2 <- centered_sample_2 %*% estimate_leading_pc
  # 
  # mean_diff <- estimate_mean_1 - estimate_mean_2
  # centered_1_times_meandiff <- centered_sample_1 %*% mean_diff
  # centered_2_times_meandiff <- centered_sample_2 %*% mean_diff
  # 
  # meandiff_times_pc <- as.numeric(crossprod(mean_diff, estimate_leading_pc))
  # two_sample_weight <- nrow(cross_fitting_sample_1)/(nrow(cross_fitting_sample_1) + nrow(cross_fitting_sample_2))
  # 
  # influence_eigenvector_each_subject_1 <- (centered_1_times_meandiff - meandiff_times_pc * centered_1_times_pc) * centered_1_times_pc
  # influence_eigenvector_each_subject_1 <- two_sample_weight * (nuisance_collection$estimate_eigenvalue - nuisance_collection$estimate_noise_variance)^(-1) * influence_eigenvector_each_subject_1
  # 
  # influence_eigenvector_each_subject_2 <- (centered_2_times_meandiff - meandiff_times_pc * centered_2_times_pc) * centered_2_times_pc
  # influence_eigenvector_each_subject_2 <- (1 - two_sample_weight) * (nuisance_collection$estimate_eigenvalue - nuisance_collection$estimate_noise_variance)^(-1) * influence_eigenvector_each_subject_2
  # 
  
  
  
  return(list(inner_product_1 = inner_product_1_all,
              inner_product_2 = inner_product_2_all,
              influence_eigenvector_each_subject_2 = influence_function_all,
              for_variance_subject_1 = for_variance_subject_1_all,
              for_variance_subject_2 = for_variance_subject_2_all))

}

