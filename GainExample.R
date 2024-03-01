# Example usage:
library(reticulate)
  n <- 2000
  d <- 3
  X <- matrix(rnorm(n * d, mean = 0, sd = 2), ncol = d, nrow = n)
  
  # M: Pattern matrix
  M <- apply(X, 2, function(x) sample(c(0, 1), size = length(x), replace = TRUE, prob = c(1 - 0.2, 0.2)))
  # X.NA: Matrix with missing values
  X_NA<-X
  X_NA[M == 1] <- NA
  
  # GAIN parameters
  gain_parameters <- list(
    batch_size = 64L,
    hint_rate = 0.9,
    alpha = 10,
    iterations = 10000L
  )
  
  # Impute missing values using GAIN
  X_imputed <- gain(X_NA, gain_parameters)
  
  save(X_imputed, file="X_imputed.Rdata")
  