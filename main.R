library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# -----------------------------------------
# The following functions are used to extract parameter estimates from a Stan object.
# -----------------------------------------
get_prod_restricted_prm <- function(prm_vec){
  return(append(1.0/prod(prm_vec), prm_vec))
}

get_mean_restricted_prm <- function(prm_vec){
  return(append(-1*sum(prm_vec), prm_vec))
}

get_multdim_param <- function(prm_vec, N, D){
  for(n in 1:N){
    prm = prm_vec[((n-1)*D+1):((n-1)*D+D)]
    if(n == 1){
      mat = t(data.frame(prm))
    } else {
      mat = rbind(mat, t(data.frame(prm)))
    }
  }  
  return(mat)
}

get_category_prm <- function(prm_vec, N, K){
  for(n in 1:N){
    prm = prm_vec[((n-1)*(K-2)+1):((n-1)*(K-2)+(K-2))]
    prm = append(prm, -1*sum(prm))
    if(n == 1){
      mat = t(data.frame(prm))
    } else {
      mat = rbind(mat, t(data.frame(prm)))
    }
  }  
  return(mat)
}
# -----------------------------------------

# Read the Stan code for estimating the multidimensional generalized MFRM.
stan <- stan_model(file = "mult_gmfrm_uto.stan")

# Read our data and create a list object containing the required information.
dat <- read.table("dat_mult_gmfrm_uto.csv", header=TRUE, sep=",")
colnames(dat) <- c("ExamineeID", "ItemID", "RaterID", "Score")
dat_stan=list(K=4, J=134, I=9, R=18, D=2, N=nrow(dat), 
               ItemID=dat$ItemID, ExamineeID=dat$ExamineeID, RaterID=dat$RaterID, X=dat$Score)

# Run MCMC.
fit <- sampling(stan, data=dat_stan, iter=4000, warmup=2000, chains=1)

# Extract rater parameter estimates.
alpha_r <- get_prod_restricted_prm(summary(fit, par="alpha_r")$summary[,"mean"])
beta_r <- get_mean_restricted_prm(summary(fit, par="beta_r")$summary[,"mean"])

# Extract rubric item parameter estimates.
alpha_i <- get_multdim_param(summary(fit, par="alpha_i")$summary[,"mean"], dat_stan$I, dat_stan$D)
beta_i <- summary(fit, par="beta_i")$summary[,"mean"]
d_im <- get_category_prm(summary(fit, par="beta_ik")$summary[,"mean"], dat_stan$I, dat_stan$K)

# Extract ability estimates.
theta <- get_multdim_param(summary(fit, par="theta")$summary[,"mean"], dat_stan$J, dat_stan$D)

# Calculate WAIC.
waic <- waic(extract_log_lik(fit))$estimates[3,1]
