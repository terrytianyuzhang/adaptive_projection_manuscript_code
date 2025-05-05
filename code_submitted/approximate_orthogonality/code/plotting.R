work_directory <- './approximate_orthogonality/'
all_result <- readRDS(glue::glue(work_directory, '/data/zero_inflated_normal.rds'))

library(ggplot2)
library(data.table)
all_result <- data.table(all_result)
all_result_long <- melt(all_result, id.vars = "repeat_index", 
                        variable.name = "stat_type", value.name = "value")

x_vals <- seq(-4, 4, length.out = 100)
y_vals <- dnorm(x_vals, mean = 0, sd = 1)
df <- data.frame(x = x_vals, y = y_vals)
p1 <- ggplot() +
  geom_histogram(aes(x = value, y = ..density.., fill = stat_type), 
                 alpha = 0.6, position = "identity", bins = 30, color = "black", data = all_result_long) +
  geom_line( aes(x = x, y = y), color = "black", size = 1, data = df) +  # Overlay standard normal curve
  scale_fill_manual(
    values = c("test_stat" = "#0073C2", "test_stat2" = "#EFC000"), 
    name = "PC index",  
    labels = c("PC1", "PC2")
  ) +
  theme_minimal(base_size = 14) +
  labs(title = NULL,
       x = "Test Statistics",  # Use LaTeX-style expression in x-axis label
       y = NULL)+
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))  


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

set.seed(42)
n <- 50*10  # Number of samples
p <- 300   # Number of dimensions

block_size <- 10  # Size of correlation blocks
cov <- c(rep(1.6, 1), rep(0.8, 1), rep(0, p/block_size - 2))# Correlation within blocks
var <- c(rep(2, 1), rep(1, 1), rep(1, p/block_size - 2))

mu1 <- rep(1,p)
mu2 <- c(rep(1, block_size), rep(2, block_size-5), rep(1, p-2*block_size+5))
zero_prob <- 0.5

X1 <- generate_zero_inflated_block_gaussian (5*n, p, mu1, block_size, var, cov, zero_prob)
marginal_df <- data.frame(X1 = X1[,1])
sum(X1 == 0)/ (5*n*p)

p2 <- ggplot(marginal_df, aes(x = X1, y = ..density..)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.7, fill = '#73C200') +
  theme_minimal(base_size = 14) +
  labs(x = "First Dimension of X", y = "Density") +
  theme(
    legend.position = "none",  # No legend needed for a single-variable histogram
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

library(ggpubr)
ggarrange(plotlist = list(p2, p1), ncol = 2, labels = c('A', 'B'))
ggsave(glue::glue(work_directory, '/data/marginal_n_test_stat.pdf'),
       width = 6, height =3.5)
