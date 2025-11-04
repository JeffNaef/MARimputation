
#################################################
##Overlap Condition
#################################################
library(ggplot2)
library(gridExtra)
library(MASS)

set.seed(123)

# Number of samples per pattern
n <- 3000

# # ===== EXAMPLE 1: Uniform distributions =====
# # Pattern m1: X2* ~ U([0,1])
# X2_m1 <- runif(n, 0, 1)
# # Pattern m2: X2* ~ U([1,2])
# X2_m2 <- runif(n, 1, 2)
# 
# # Common noise epsilon ~ U([0,1])
# epsilon_m1 <- runif(n, 0, 1)
# epsilon_m2 <- runif(n, 0, 1)
# 
# # X1* = X2* * epsilon, so X1*|X2* ~ U([0, X2*])
# X1_m1 <- X2_m1 * epsilon_m1
# X1_m2 <- X2_m2 * epsilon_m2
# 
# # Create data frame for Example 1
# df1 <- data.frame(
#   X1 = c(X1_m1, X1_m2),
#   X2 = c(X2_m1, X2_m2),
#   Pattern = factor(rep(c("M = m_1", "M = m_2"), each = n))
# )
# 
# # Plot Example 1 - Pattern m2 with X1 missing
# p1 <- ggplot() +
#   # Plot pattern m1 (complete data) as points
#   geom_point(data = df1[df1$Pattern == "M = m_1", ], 
#              aes(x = X2, y = X1, color = Pattern), 
#              alpha = 0.5, size = 2) +
#   # Plot pattern m2 as circles on x-axis (X1 is missing)
#   geom_point(data = df1[df1$Pattern == "M = m_2", ], 
#              aes(x = X2, y = 0, color = Pattern), 
#              shape = 1, size = 3, alpha = 0.7, stroke = 1.2) +
#   
#   scale_color_manual(values = c("M = m_1" = "#E69F00", "M = m_2" = "#56B4E9")) +
#   labs(
#     #title = "Example 1: Uniform Distributions",
#     #subtitle = "Pattern m_1: complete data | Pattern m_2: X₁ missing (only X₂ observed)",
#     x = expression(X[2]),
#     y = expression(X[1])
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "bottom",
#     plot.title = element_text(face = "bold"),
#     panel.grid.minor = element_blank()
#   ) +
#   coord_cartesian(xlim = c(-0.1, 2.1), ylim = c(-0.15, 2.1))


# ===== EXAMPLE 1: Gaussian mixture =====
# Pattern m1: (X1, X2) | M=m1 ~ N(mean1, Sigma1)
mean1 <- c(0, 0)
Sigma1 <- matrix(c(2, 1, 1, 1), nrow = 2)

# Pattern m2: (X1, X2) | M=m2 ~ N(mean2, Sigma2)
mean2 <- c(10, 10)
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
p1 <- ggplot() +
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
  coord_cartesian(ylim = c(-5, 5), xlim = c(-5, 10))




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



#Plotting:
# Open a PNG device for saving the plot
png(filename = paste0("Overlap_Illustration", ".png"), 
    width = 1000,    # Width in pixels
    height = 600,    # Height in pixels
    res = 120)       # Resolution in dpi

# Combine plots
grid.arrange(p1, p2, ncol = 2)

# Close the PNG device
dev.off()





#################################################
## Example \ref{interesting_new_Example}
#################################################

n<-40000
d=3

X<-simulate_fgm(n=n, alpha=1)
X<-cbind(X,matrix(runif( (d-2)*n ), nrow=n, ncol=d-2 ))

colnames(X)<-paste0("X",1:d)


vectors <- matrix(c(
  rep(0, d),
  0, 1, rep(0,d-2),
  rep(0,d-1), 1,
  1,1, rep(0,d-2)
), nrow = 4, byrow = TRUE)


# Generate random draws
# sample() will generate indices, which we use to select rows from the matrix
M <- vectors[apply(X,1, function(x) sample(1:4, size = 1, prob=c((x[1]+x[2])/3, (1-x[1])/3, (1-x[2])/3, 1/3), replace = TRUE)), ]
X.NA<-X
X.NA[M==1]<-NA

nLm<-nrow(X[ (M[,1]==0 & M[,2]==0 & M[,3]==0) |  (M[,1]==0 & M[,2]==0 & M[,3]==1) ,1:2])

# Create grid for density contours
x1_grid <- seq(0, 1, length.out = 100)
x2_grid <- seq(0, 1, length.out = 100)
density_grid <- outer(x1_grid, x2_grid, function(x1, x2) {
  1 + (2*x1 - 1) * (2*x2 - 1)
})

density_grid2 <- outer(x1_grid, x2_grid, function(x1, x2) {
  1 + (2*x1 - 1) * (2*x2 - 1) * (2/3)*(1+x1)
})

# Plot with density overlay
par(mfrow=c(1,2))

# Plot 1: M=m_4 (random sample)
plot(X[sample(1:n, size=nLm, replace = F), 1:2], 
     pch = 16, cex = 0.5, col = rgb(0, 0, 1, 0.3),
     xlab = "X1", ylab = "X2")
contour(x1_grid, x2_grid, density_grid, 
        add = TRUE, col = "red", lwd = 2, labcex = 0.8)

# Plot 2: M in Lm
plot(X.NA[(M[,1]==0 & M[,2]==0 & M[,3]==0) | (M[,1]==0 & M[,2]==0 & M[,3]==1), 1:2], 
     pch = 16, cex = 0.5, col = rgb(0, 0, 1, 0.3),
     xlab = "X1", ylab = "X2")
contour(x1_grid, x2_grid, density_grid2, 
        add = TRUE, col = "red", lwd = 2, labcex = 0.8)

par(mfrow=c(1,1))

##Alternative:
par(mfrow=c(1,3))

hist(X.NA[(M[,1]==0 & M[,2]==0 & M[,3]==0) | (M[,1]==0 & M[,2]==0 & M[,3]==1), 1], xlab="X_1", main="")
hist(X.NA[(M[,1]==0 & M[,2]==0 & M[,3]==0) | (M[,1]==0 & M[,2]==0 & M[,3]==1), 2], xlab="X_2", main="")

# Plot 2: M in Lm
plot(X.NA[(M[,1]==0 & M[,2]==0 & M[,3]==0) | (M[,1]==0 & M[,2]==0 & M[,3]==1), 1:2], 
     pch = 16, cex = 0.5, col = rgb(0, 0, 1, 0.3),
     xlab = "X_1", ylab = "X_2")
contour(x1_grid, x2_grid, density_grid2, 
        add = TRUE, col = "red", lwd = 2, labcex = 0.8)

par(mfrow=c(1,1))




###Saving the plot

#Plotting:
# Open a PNG device for saving the plot
png(filename = paste0("Example_X1X2", ".png"), 
    width = 1200,    # Width in pixels
    height = 400,    # Height in pixels
    res = 120)       # Resolution in dpi

##Alternative:
par(mfrow=c(1,3))

hist(X.NA[(M[,1]==0 & M[,2]==0 & M[,3]==0) | (M[,1]==0 & M[,2]==0 & M[,3]==1), 1], xlab="X_1", main="", cex.axis=1.5,cex.lab=1.5)
hist(X.NA[(M[,1]==0 & M[,2]==0 & M[,3]==0) | (M[,1]==0 & M[,2]==0 & M[,3]==1), 2], xlab="X_2", main="",cex.axis=1.5,cex.lab=1.5)

# Plot 2: M in Lm
plot(X.NA[(M[,1]==0 & M[,2]==0 & M[,3]==0) | (M[,1]==0 & M[,2]==0 & M[,3]==1), 1:2], 
     pch = 16, cex = 0.5, col = rgb(0, 0, 1, 0.3),
     xlab = "X_1", ylab = "X_2", cex.axis=1.5,cex.lab=1.5)
contour(x1_grid, x2_grid, density_grid2, 
        add = TRUE, col = "red", lwd = 2, labcex = 0.8)

par(mfrow=c(1,1))


# Close the PNG device
dev.off()



#################################################################################
###Example 1: Showing that all patterns need to be used for imputation##
################################################################################
library(mice)
library(miceDRF)
library(drf)
library(MASS)
library(miceDRF)
source("helpers.R")


set.seed(10)
n<-10000

Xstar<-matrix( runif(n=3*n), nrow=n, ncol=3    )
Mindex<-sapply(1:n, function(i) 
  sample(1:3, size=1, replace=F, prob=c(Xstar[i,1]/3, 2/3-Xstar[i,1]/3, 1/3    )))
M<-t(sapply(1:n, function(i) if (Mindex[i]==1){c(0,0,0)}else if (Mindex[i]==2) {c(0,1,0)} else if  (Mindex[i]==3) {c(1,0,0)}) )


X<-Xstar
X[Mindex==2,2]<-NA
X[Mindex==3,1]<-NA

head(cbind(X,M,Mindex))



hist(Xstar[,2])


png(filename = "Example1.png", 
    width = 1200,    # Width in pixels
    height = 500,    # Height in pixels
    res = 120)       # Resolution in dpi

par(mfrow=c(1,3))

# What we want to impute = unconditional distribution
hist(Xstar[Mindex==3,1], xlab="", main="", probability=T)

# Fully observed pattern:
hist(Xstar[Mindex==1,1], xlab="", main="", probability=T)

# Pattern m_2
hist(Xstar[Mindex==2,1], xlab="", main="", probability=T)

## Learning from both patterns
#hist(Xstar[Mindex==1 | Mindex==2,1], xlab="",main="", probability=T)

#left: Distribution we want to impute X_2 \mid M=m_3
#middle: Distribution of X_2 in the fully observed pattern (X_2 \mid M=m_1)
#right: Distribution in all observed patterns (Mixture of X_2 \mid M=m_1 and X_2 \mid M=m_2)
# Close the PNG device
dev.off()


################################################
###Example 2: Showing MAR allows for changes in observed distribution
#####################################
library(MASS)
source("helpers.R")
require(drf)
n<-10000
U<-runif(n)

## In both cases: X_1|X_2 ~ N(X_2, 1)

X2m1<-mvrnorm(n= sum(U <= 1/2),  mu=0, Sigma=1  )
X1m1<-1*X2m1 + rnorm(n= sum(U <= 1/2))
X2m2<-mvrnorm(n= sum(U > 1/2),  mu=5, Sigma=1  )
X1m2<-1*X2m2 + rnorm(n= sum(U > 1/2))

Xstar<-rbind( cbind(X1m1, X2m1), cbind(X1m2, X2m2)   )

par(mfrow=c(1,2))
hist( X2m1, xlab="", main="")
hist( X2m2, xlab="", main="")

#Left: Distribution of the observed X_2 in pattern 1 (X_2 \mid M=m_1)
#Right: Distribution of the observed X_2 in pattern 2 (X_2 \mid M=m_2)





##Try to impute this Example ##
##This is not ideal however, since RF tends to be bad with only one X!
set.seed(10)
n<-2000
U<-runif(n)

## In both cases: X_1|X_2 ~ N(X_2, 1)
##Careful: Compared to the paper the two patterns are reversed:
## m_1=(1,0), m_2=(0,1)
## (X_1, X_2) | M=m_1 ~ N((5,5), Sigma  )
## (X_1, X_2) | M=m_2 ~ N((0,0), Sigma  )
X2m1<-mvrnorm(n= sum(U <= 1/2),  mu=5, Sigma=1  )
X1m1<-1*X2m1 + rnorm(n= sum(U <= 1/2))
X2m2<-mvrnorm(n= sum(U > 1/2),  mu=0, Sigma=1  )
X1m2<-1*X2m2 + rnorm(n= sum(U > 1/2))


Xstar<-rbind( cbind(X1m1, X2m1), cbind(X1m2, X2m2)   )

X.NA<-Xstar
X.NA[1:nrow(X1m1),1]<-NA

imputations<- doimputation(X.NA=X.NA, methods="cart", parallelize=F, m=1, min.node.size=15)
impcart<-imputations$imputations$cart$'1'

# impdrf (reuse code!)
imputations<- doimputation(X.NA=X.NA, methods="DRF", parallelize=F, m=1, min.node.size=15)
impDRF<-imputations$imputations$DRF$'1'


reg<- lm(X1m2~X2m2)
Ysample<-reg$coefficients[1] + reg$coefficients[2]*X2m1+ rnorm(nrow(X2m1),0, sd=var(reg$residuals))
impreg<-X.NA
impreg[1:nrow(X1m1),1]<-Ysample



# DRF<-drf(Y=X1m2, X=X2m2, num.trees=10, min.node.size=15, num.features = 100,
#          compute.oob.predictions = F)
# HhatDRF<-predict(DRF, newdata=X2m1)$weights
# 
# Ygen<-sapply(1:nrow(X1m1), function(i) X1m2[sample(1:nrow(X1m2), size=1, replace = T, HhatDRF[i,])])
# impDRF<-X.NA
# impDRF[1:nrow(X1m1),1]<-Ygen

# ## Robust DRF trial = > Doesnt really work ###
# X2m2rob<- (X2m2 - mean(X2m2))/abs(X2m2-mean(X2m2))
# DRFrob<-drf(Y=X1m2, X=X2m2rob, num.trees=2000, min.node.size=15)
# X2m1rob<- (X2m1 - mean(X2m1))/abs(X2m1-mean(X2m1))
# HhatDRFrob<-predict(DRFrob, newdata=X2m1rob)$weights
# 
# Ygen<-sapply(1:nrow(X1m1), function(i) X1m2[sample(1:nrow(X1m2), size=1, replace = T, HhatDRFrob[i,])])
# impDRFrob<-X.NA
# impDRFrob[1:nrow(X1m1),1]<-Ygen






# engr = engression(X=X2m2,Y=X1m2,num_epochs = 100)
# Ysample = predict(engr,X2m1,type="sample",nsample=1)
# impeng<-X.NA
# impeng[1:nrow(X1m1),1]<-Ysample

#linreg<-lm(X1m2~ X2m2 )
#X1m2rob<- X1m2-cbind(rep(1,nrow(X2m2)), X2m2)%*%linreg$coefficients
#DRFlin<-drf(Y=X1m2rob, X=X2m2, num.trees=2000, min.node.size=15)
#HhatDRFlin<-predict(DRFlin, newdata=X2m1)$weights

#Ygen0<-sapply(1:nrow(X1m1), function(i) X1m2rob[sample(1:nrow(X1m2rob), size=1, replace = T, HhatDRFlin[i,])])
#Ygen <- Ygen0 + cbind(rep(1,nrow(X2m1)), X2m1)%*%linreg$coefficients
#HhatDRFlin<-X.NA
#HhatDRFlin[1:nrow(X1m1),1]<-Ygen

# Open a PNG device for saving the plot
png(filename = "Example2.png", 
    width = 1200,    # Width in pixels
    height = 700,    # Height in pixels
    res = 120)       # Resolution in dpi

par(mfrow=c(2,2))
hist( X1m1, xlab="", ylab="", main="Truth", probability = T, ylim=c(0,1.2), xlim=c(-5,10))
# Almost perfect imputation!
hist( impreg[1:nrow(X1m1),1], xlab="", ylab="", main="mice-norm.nob", probability = T, xlim=c(-5,10), ylim=c(0,1.2))
hist( impcart[1:nrow(X1m1),1], xlab="", ylab="", main="mice-cart", probability = T, xlim=c(-5,10), ylim=c(0,1.2))
hist( impDRF[1:nrow(X1m1),1], xlab="", ylab="", main="mice-DRF", probability = T, xlim=c(-5,10), ylim=c(0,1.2))
#hist( impDRFrob[1:nrow(X1m1),1], xlab="", ylab="", main="DRF Robust", probability = T, xlim=c(-5,10), ylim=c(0,0.3))
#hist( impeng[1:nrow(X1m1),1], xlab="", ylab="", main="Engression", probability = T, xlim=c(-5,10), ylim=c(0,0.3))
# hist( impDRF[1:nrow(X1m1),1],probability = T, add=T, col=rgb(0,1,0,0.3))
# hist( impcart[1:nrow(X1m1),1], probability = T, add=T, col=rgb(1,0,1,0.3))
# hist( impeng[1:nrow(X1m1),1], probability = T, add=T, col=rgb(1,0.5,1,0.3))


dev.off()


