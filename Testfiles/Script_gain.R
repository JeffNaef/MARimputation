#install.packages("reticulate")
library(reticulate)

#Sys.setenv("gain_env" =  path.expand("~/anaconda3/envs/gain_env"))
#use_python(path.expand("~/opt/anaconda3/envs/gain_env/bin/python"))
#use_condaenv(path.expand("~/opt/anaconda3/envs/gain_env"))


##Laptop
Sys.setenv("gain_env" =  "C:/Users/jeffr/anaconda3/envs/gain_env")
#use_python("C:/Users/jeffr/anaconda3/envs/gain_env/bin/python")
use_condaenv("C:/Users/jeffr/anaconda3/envs/gain_env")

## Write 
# conda activate gain_env
## in the terminal!

py_config()

#Required Python Packages
tensorflow <- import("tensorflow")
numpy <- import("numpy")
tqdm <- import("tqdm")
keras <- import("keras")
argparse <- import("argparse")  #pip install argparse
sys<- import("sys")
reticulate::source_python("gain.py") #there will be  warning but don't worry

# Example usage:
n <- 2000
d <- 10
X <- matrix(rnorm(n * d, mean = 0, sd = 2), ncol = d, nrow = n)

# M: Pattern matrix
M <- apply(X, 2, function(x) sample(c(0, 1), size = length(x), replace = TRUE, prob = c(1 - 0.2, 0.2)))
# X.NA: Matrix with missing values
X_NA<-X
X_NA[M == 1] <- NA






# GAIN parameters
gain_parameters =  list(
  batch_size = 64L,
  hint_rate = 0.9,
  alpha = 10, 
  iterations = 10000L
)

# Impute missing values using GAIN
X_imputed <- gain(X_NA, gain_parameters)

#save(X_imputed, file="X_imputed.Rdata")


#In Terminal:
# conda install pytorch::pytorch torchvision torchaudio -c pytorch
# conda install  -c conda-forge scikit-learn 
# conda install  -c conda-forge scipy
# conda install pandas

#Required Python Packages
torch <- import("torch") 
numpy <- import("numpy")
scipy <- import("scipy")
pandas <- import("pandas")
sklearn<- import("sklearn")
torchvision <- import("torchvision")

reticulate::source_python("MIWAE_Pytorch.py") #there will be  warning but don't worry


X_miss = X_NA
X_NA[is.na(X_NA)] <- 0
mask <- is.finite(X_miss)

test = miwae(X_NA, X, mask)

X_imputed_miwae = test[1][[1]]

## without X, I hope this is still correct!

test = miwae(X_NA, X_NA, mask)
X_imputed_miwae2 = test[1][[1]]

