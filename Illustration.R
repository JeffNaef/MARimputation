library(ggplot2)
library(gridExtra)

set.seed(123)

# Number of samples per pattern
n <- 5000

# ===== EXAMPLE 1: Uniform distributions =====
# Pattern m1: X2* ~ U([0,1])
X2_m1 <- runif(n, 0, 1)
# Pattern m2: X2* ~ U([1,2])
X2_m2 <- runif(n, 1, 2)

# Common noise epsilon ~ U([0,1])
epsilon_m1 <- runif(n, 0, 1)
epsilon_m2 <- runif(n, 0, 1)

# X1* = X2* * epsilon, so X1*|X2* ~ U([0, X2*])
X1_m1 <- X2_m1 * epsilon_m1
X1_m2 <- X2_m2 * epsilon_m2

# Create data frame for Example 1
df1 <- data.frame(
  X1 = c(X1_m1, X1_m2),
  X2 = c(X2_m1, X2_m2),
  Pattern = factor(rep(c("M = m_1", "M = m_2"), each = n))
)

# Plot Example 1
p1 <- ggplot(df1, aes(x = X2, y = X1, color = Pattern)) +
  geom_point(alpha = 0.5, size = 2) +
  # Add the conditional relationship line: X1 ranges from 0 to X2
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "black", size = 1.2) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed", 
              color = "black", size = 1.2) +
  scale_color_manual(values = c("M = m_1" = "#E69F00", "M = m_2" = "#56B4E9")) +
  labs(
    #title = "Example 1: Uniform Distributions",
    #subtitle = "X₁|X₂ ~ U([0, X₂]) is invariant across patterns",
    x = expression(X[2]),
    y = expression(X[1])
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(-0.1, 2.1), ylim = c(-0.1, 2.1))

# ===== EXAMPLE 2: Gaussian mixture =====
library(MASS)

# Pattern m1: (X1, X2) | M=m1 ~ N(mean1, Sigma1)
mean1 <- c(0, 0)
Sigma1 <- matrix(c(2, 1, 1, 1), nrow = 2)

# Pattern m2: (X1, X2) | M=m2 ~ N(mean2, Sigma2)
mean2 <- c(5, 5)
Sigma2 <- matrix(c(2, 1, 1, 1), nrow = 2)

# Generate samples
samples_m1 <- mvrnorm(n, mu = mean1, Sigma = Sigma1)
samples_m2 <- mvrnorm(n, mu = mean2, Sigma = Sigma2)

# Create data frame for Example 2
df2 <- data.frame(
  X1 = c(samples_m1[, 1], samples_m2[, 1]),
  X2 = c(samples_m1[, 2], samples_m2[, 2]),
  Pattern = factor(rep(c("M = m_1", "M = m_2"), each = n))
)

# Calculate the conditional relationship: X1|X2 for both patterns
# For bivariate normal, the conditional mean is: E[X1|X2] = μ1 + Σ12/Σ22 * (X2 - μ2)
# With our covariance matrix, Σ12 = 1, Σ22 = 1, so slope = 1

# For m1: E[X1|X2] = 0 + 1*(X2 - 0) = X2
# For m2: E[X1|X2] = 5 + 1*(X2 - 5) = X2
# So both have the same conditional relationship: E[X1|X2] = X2

# Plot Example 2
p2 <- ggplot(df2, aes(x = X2, y = X1, color = Pattern)) +
  geom_point(alpha = 0.5, size = 2) +
  # Add the conditional relationship line: E[X1|X2] = X2
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "black", size = 1.2) +
  scale_color_manual(values = c("M = m_1" = "#E69F00", "M = m_2" = "#56B4E9")) +
  labs(
    #title = "Example 2: Gaussian Mixture",
    #subtitle = "E[X₁|X₂] = X₂ is invariant across patterns",
    x = expression(X[2]),
    y = expression(X[1])
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Combine plots
grid.arrange(p1, p2, ncol = 2)



## Second try
#################################################
library(ggplot2)
library(gridExtra)
library(MASS)

set.seed(123)

# Number of samples per pattern
n <- 3000

# ===== EXAMPLE 1: Uniform distributions =====
# Pattern m1: X2* ~ U([0,1])
X2_m1 <- runif(n, 0, 1)
# Pattern m2: X2* ~ U([1,2])
X2_m2 <- runif(n, 1, 2)

# Common noise epsilon ~ U([0,1])
epsilon_m1 <- runif(n, 0, 1)
epsilon_m2 <- runif(n, 0, 1)

# X1* = X2* * epsilon, so X1*|X2* ~ U([0, X2*])
X1_m1 <- X2_m1 * epsilon_m1
X1_m2 <- X2_m2 * epsilon_m2

# Create data frame for Example 1
df1 <- data.frame(
  X1 = c(X1_m1, X1_m2),
  X2 = c(X2_m1, X2_m2),
  Pattern = factor(rep(c("M = m_1", "M = m_2"), each = n))
)

# Plot Example 1 - Pattern m2 with X1 missing
p1 <- ggplot() +
  # Plot pattern m1 (complete data) as points
  geom_point(data = df1[df1$Pattern == "M = m_1", ], 
             aes(x = X2, y = X1, color = Pattern), 
             alpha = 0.5, size = 2) +
  # Plot pattern m2 as circles on x-axis (X1 is missing)
  geom_point(data = df1[df1$Pattern == "M = m_2", ], 
             aes(x = X2, y = 0, color = Pattern), 
             shape = 1, size = 3, alpha = 0.7, stroke = 1.2) +
  
  scale_color_manual(values = c("M = m_1" = "#E69F00", "M = m_2" = "#56B4E9")) +
  labs(
    #title = "Example 1: Uniform Distributions",
    #subtitle = "Pattern m_1: complete data | Pattern m_2: X₁ missing (only X₂ observed)",
    x = expression(X[2]),
    y = expression(X[1])
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(-0.1, 2.1), ylim = c(-0.15, 2.1))

# ===== EXAMPLE 2: Gaussian mixture =====
# Pattern m1: (X1, X2) | M=m1 ~ N(mean1, Sigma1)
mean1 <- c(0, 0)
Sigma1 <- matrix(c(2, 1, 1, 1), nrow = 2)

# Pattern m2: (X1, X2) | M=m2 ~ N(mean2, Sigma2)
mean2 <- c(5, 5)
Sigma2 <- matrix(c(2, 1, 1, 1), nrow = 2)

# Generate samples
samples_m1 <- mvrnorm(n, mu = mean1, Sigma = Sigma1)
samples_m2 <- mvrnorm(n, mu = mean2, Sigma = Sigma2)

# Create data frame for Example 2
df2 <- data.frame(
  X1 = c(samples_m1[, 1], samples_m2[, 1]),
  X2 = c(samples_m1[, 2], samples_m2[, 2]),
  Pattern = factor(rep(c("M = m_1", "M = m_2"), each = n))
)

# Plot Example 2 - Pattern m2 with X1 missing
p2 <- ggplot() +
  # Plot pattern m1 (complete data) as points
  geom_point(data = df2[df2$Pattern == "M = m_1", ], 
             aes(x = X2, y = X1, color = Pattern), 
             alpha = 0.5, size = 2) +
  # Plot pattern m2 as circles on x-axis (X1 is missing)
  geom_point(data = df2[df2$Pattern == "M = m_2", ], 
             aes(x = X2, y = 0, color = Pattern), 
             shape = 1, size = 3, alpha = 0.7, stroke = 1.2) +
  # Add the conditional relationship line: E[X1|X2] = X2
  scale_color_manual(values = c("M = m_1" = "#E69F00", "M = m_2" = "#56B4E9")) +
  labs(
    #title = "Example 2: Gaussian Mixture",
    #subtitle = "Pattern m_1: complete data | Pattern m_2: X₁ missing (only X₂ observed)",
    x = expression(X[2]),
    y = expression(X[1])
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(-5, 5), xlim = c(-5, 7))

# Combine plots
grid.arrange(p1, p2, ncol = 2)










