# Below is the link to Vivek's code
#https://colab.research.google.com/drive/10ADmaLN83mOHdUQ7hMSgsjUW9Fka4G--?usp=sharing#scrollTo=Zfs4BXQ8iBF9
library(glmnet)
library(caret)
library(pracma)
library(MASS)
library(grpreg)
library(PMA)
# ---------------------------
# Utility Functions
# ---------------------------

index_spliter <- function(array, n_folds = 5){
  
  # array <- 1:99
  
  # Calculate the length of each part
  part_length <- length(array) %/% n_folds
  
  # Create an empty list to store the parts
  parts <- vector("list", n_folds)
  
  # Randomly shuffle the array
  shuffled_array <- sample(array)
  
  # Split the shuffled array into parts
  for (fold_index in 1:n_folds) {
    start_index <- (fold_index - 1) * part_length + 1
    
    if(fold_index < n_folds){
      end_index <- fold_index * part_length
    }else{
      end_index <- length(array)
    }
    
    parts[[fold_index]] <- shuffled_array[start_index:end_index]
  }
  
  return(parts)
}

validate_and_convert_data <- function(data, name) {
  if (!inherits(data, c("matrix", "data.frame"))) {
    stop(paste(name, "should be a matrix or a data frame."))
  }
  
  return(as.matrix(data))
}

check_non_null_and_identical_colnames <- function(data_list) {
  # Check for null or empty column names and ensure all column names are identical
  colnames_data <- lapply(data_list, colnames)
  
  # Check for null or empty column names
  for (i in 1:length(colnames_data)) {
    if (any(is.null(colnames_data[[i]])) || any(colnames_data[[i]] == "")) {
      stop(paste("Dataset", i, "contains null or empty column names. Please give the input data proper column names."))
    }
  }
  
  # Check if all column names are identical
  if (!all(sapply(colnames_data, function(x) all(x == colnames_data[[1]])))) {
    stop("The column names across the datasets are not identical. Please make sure all datasets have the same column names.")
  }
}

normalize_and_split <- function(df1, df2) {
  # Combine
  combined <- rbind(df1, df2)
  
  # Compute pooled mean and SD
  pooled_mean <- colMeans(combined)
  pooled_sd <- apply(combined, 2, sd)
  
  # Center and scale
  normalized <- scale(combined, center = pooled_mean, scale = pooled_sd)
  
  # Split back
  df1_norm <- normalized[1:nrow(df1), , drop = FALSE]
  df2_norm <- normalized[(nrow(df1) + 1):nrow(combined), , drop = FALSE]
  
  return(list(df1 = df1_norm, df2 = df2_norm))
}

check_data_for_folds <- function(data, n_folds) {
  if (nrow(data) < n_folds) stop("Not enough rows to create folds.")
}

# ---------------------------
# Statistical and Modeling Helpers
# ---------------------------

fit_lasso <- function(control_train, treat_train,
                      lambda_type = c("lambda.min", "lambda.1se"),
                      classifier_method = c("lasso", "group_lasso"),
                      group = NULL) {
  
  lambda_type <- match.arg(lambda_type)
  classifier_method <- match.arg(classifier_method)
  
  X_train <- rbind(control_train, treat_train)  
  y_train <- c(rep(0, nrow(control_train)), rep(1, nrow(treat_train)))  
  
  if (classifier_method == "group_lasso") {
    if (is.null(group)) stop("Group vector must be provided for group lasso.")
    
    lasso_result <- cv.grpreg(X = X_train, 
                              y = y_train, 
                              group = group, 
                              penalty = "grLasso", 
                              family = 'binomial')
    
  } else {
    lasso_result <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1)
  }
  # Step 3: Extract LASSO coefficients at the optimal lambda
  # Extract coefficients at the optimal lambda
  beta_est <- coef(lasso_result, s = lambda_type)[-1]
  names(beta_est) <- colnames(X_train)
  
  # Calculate threshold for small coefficients
  
  n_effect <- 2 * min(nrow(control_train), nrow(treat_train))
  max_beta_element <- max(abs(beta_est))
  
  # Vectorized operation to threshold small coefficients to zero
  threshold <- max_beta_element * n_effect^(-1/3)
  beta_est[abs(beta_est) < threshold] <- 0
  
  # if(length(beta_est) == 50) beta_est[6:50] <- 0
  
  return(beta_est)
  
}

estimate_leading_pc <- function(control, pca_method = c("dense_pca", "sparse_pca")) {
  # Center the data
  centered_control <- sweep(as.data.frame(control), 2, colMeans(control))
  sample_centered <- as.matrix(centered_control)
  feature_number <- ncol(sample_centered)
  
  # Match PCA method argument
  pca_method <- match.arg(pca_method)
  print(pca_method)
  # Handle edge case: 1-dimensional data
  if (feature_number == 1) {
    warning("There is only one dimension and PCA is requested.")
    pc <- matrix(1, ncol = 1)
    names(pc) <- colnames(control)
    return(pc / sqrt(sum(pc^2)))
  }
  
  # Force dense PCA for low-dimensional data
  if (feature_number <= 30) {
    message("Dimension too small, switching method to dense PCA.")
    pca_method <- "dense_pca"
  }
  
  # Run PCA
  if (pca_method == "dense_pca") {
    print('conducting dense PCA')
    pc <- array(irlba::irlba(sample_centered, nv = 1)$v)
  } else if (pca_method == "sparse_pca") {
    print('conducting sparse PCA')
    cv_result <- PMA::SPC.cv(sample_centered)
    pc <- PMA::SPC(
      sample_centered,
      K = 1,
      sumabsv = cv_result$bestsumabsv
    )$v
  }
  
  # Normalize and return
  names(pc) <- colnames(control)
  return(pc / sqrt(sum(pc^2)))
}

# ---------------------------
# Visualization Functions
# ---------------------------

collect_active_features <- function(test_result, voting_method = c("majority_voting"), 
                                    group = NULL, group_threshold = 1) {
  fold_data <- test_result$fold_data
  voting_method <- match.arg(voting_method)
  n_folds <- length(fold_data)
  active_features_list <- vector("list", n_folds)
  
  # Collect non-zero features for each fold
  for (i in 1:n_folds) {
    if (!is.null(fold_data[[i]]$final_beta)) {
      beta <- fold_data[[i]]$final_beta
    } else {
      beta <- fold_data[[i]]$classifier_coef
    }
    
    non_zero_features <- names(beta[abs(beta) > 1e-10])
    active_features_list[[i]] <- non_zero_features
  }
  
  # Flatten and count
  all_active_features <- unlist(active_features_list)
  feature_counts <- table(all_active_features)
  
  # Apply majority voting
  if (voting_method == 'majority_voting') {
    active_features <- names(feature_counts[feature_counts > n_folds / 2])
  }
  
  # Group handling
  if (!is.null(group)) {
    if (is.null(names(group))) {
      if (!is.null(fold_data[[1]]$final_beta)) {
        names(group) <- names(fold_data[[1]]$final_beta)
      } else {
        names(group) <- names(fold_data[[1]]$classifier_coef)
      }
    }
    
    group_nonzero_counts <- table(group[active_features])
    active_groups <- names(group_nonzero_counts[group_nonzero_counts >= group_threshold])
    
    return(list(
      active_features = active_features,
      active_groups = active_groups
    ))
  }
  
  return(active_features)
}
