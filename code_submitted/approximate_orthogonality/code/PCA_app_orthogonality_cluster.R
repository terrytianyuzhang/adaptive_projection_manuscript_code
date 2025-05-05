library(MASS)  # For mvrnorm
library(Matrix)  # For block matrix construction
work_directory <- './approximate_orthogonality'

generate_zero_inflated_block_gaussian <- function(n, p, mu, block_size, var, cov, zero_prob) {
  # n: Number of samples
  # p: Number of dimensions (features)
  # mu: Mean vector (length p)
  # block_size: Size of each block in correlation structure
  # rho: Vector of within-block correlations (length at least ceiling(p/block_size))
  # zero_prob: Probability of zero inflation (0 to 1)
  
  # Number of blocks
  num_blocks <- ceiling(p / block_size)
  
  # Construct block correlation matrix
  Sigma <- matrix(0, p, p)
  for (i in seq_len(num_blocks)) {
    idx <- ((i - 1) * block_size + 1):min(i * block_size, p)
    Sigma[idx, idx] <- cov[i]  # Set within-block correlation
    diag(Sigma)[idx] <- var[i]  # Set diagonal to 1
  }
  
  # Ensure positive definiteness
  Sigma <- as.matrix(nearPD(Sigma)$mat)
  
  # Generate multivariate Gaussian samples
  Z <- mvrnorm(n, mu = mu, Sigma = Sigma)
  
  # Apply zero inflation
  mask <- matrix(rbinom(n * p, size = 1, prob = 1 - zero_prob), n, p)
  Z <- Z * mask
  Z[Z<0.5] <- 0
  return(Z)
}

# Example usage
set.seed(42)
n <- 50  # Number of samples
p <- 300   # Number of dimensions

block_size <- 10  # Size of correlation blocks
cov <- c(rep(1.6, 1), rep(0.8, 1), rep(0, p/block_size - 2))# Correlation within blocks
var <- c(rep(2, 1), rep(1, 1), rep(1, p/block_size - 2))

mu1 <- rep(1,p)
mu2 <- c(rep(1, block_size), rep(2, block_size-5), rep(1, p-2*block_size+5))
zero_prob <- 0.5
library(HMC)
all_result <- data.frame()
for(repeat_index in 1:1000){
  print(repeat_index)
  X1 <- generate_zero_inflated_block_gaussian (5*n, p, mu1, block_size, var, cov, zero_prob)
  X2 <- generate_zero_inflated_block_gaussian (n, p, mu2, block_size, var, cov, zero_prob)
  
  test_result <- simple_pc_testing(sample_1 = X1,
                                   sample_2 = X2, 
                                   pca_method = 'sparse_pca', num_latent_factor = 2, n_folds = 5)
  test_result$test_statistics
  
  direction_indicator <- test_result$split_data[[1]]$nuisance_collection$estimate_leading_pc[,2][block_size + 1]
  
  
  all_result <- rbind(all_result,
                      data.frame(repeat_index = repeat_index,
                                 test_stat = test_result$test_statistics[1],
                                 test_stat2 = sign(direction_indicator)* test_result$test_statistics[2]))
  if(repeat_index %% 100 == 0)
    saveRDS(all_result, glue::glue(work_directory, '/data/zero_inflated_normal.rds'))
}
