# ######################################################################
# - This is a R script to fit a meta-analysis using a separate stan file.
# - Created by Tetsuro Oda
# ######################################################################

# ---------------------------------------------------
# Preparation
# ---------------------------------------------------

library(tidyverse)
library(rstan)
source("create_datasets.R")

testname <- "TE"
fibrosis_grade <- "F2"
prior_type <- 1

# Read data
d <- readxl::read_excel("NMA_data_supp.xlsx",sheet = 1)

# Include studies whose cutoffs are not defined to get 90 % of sens or spe
d2<- 
  d %>%
  filter(Flag_Threshold_sesp_90 == 0) 

# Create a dataset for meta-analysis
df <- create_dataset_separate_MA(d2,testname,fibrosis_grade)

# Create a dataset for stan
data_stan <- list(K = nrow(df),
                  J = 2,
                  TP = df$TP,
                  FP = df$FP,
                  FN = df$FN,
                  TN = df$TN,
                  prior_type = prior_type
)

# ---------------------------------------------------
# Fitting
# ---------------------------------------------------

# Compile a model
model <- stan_model(file = "ma_fit_supp.stan")

# Set some options
options(mc.cores = parallelly::availableCores()-1)
rstan_options("auto_write" = TRUE) 
rstan_options(javascript = FALSE)

# Set values for some arguments in sampling()
iter_val <- 30000
warmup_val <- 5000
chains_val <- 4
thin_val <- 10

# MCMC
fit <- sampling(
  model, 
  data = data_stan, 
  seed = 9999,
  iter = iter_val,
  warmup = warmup_val,
  chains = chains_val,
  thin = thin_val,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 12)
)

# Store large outputs
parms <- c("mu",  "sigma", "Omega_b[1,2]")
ms <- rstan::extract(fit)
s1 <- summary(fit, pars = parms)
s2 <- summary(fit, pars = c("pool_P","P"))
s_pred <- summary(fit, pars = c("sesp_pred"))

# ---------------------------------------------------
# Check MCMC samples
# ---------------------------------------------------

# Paired plots
# These plots do not include warm-ups
pairs(fit, pars = parms)
pairs(fit, pars = c("mu[1]","sigma[1]","theta[2,1]"))
pairs(fit, pars = c("mu[2]","sigma[2]","theta[1,2]"))

# Autocorrelations
stan_ac(fit, pars = parms)

# Density plots
stan_dens(fit, pars = parms, separate_chains = T)

# Stan diagnosis
stan_diag(fit)

# Trace plots
stan_trace(fit, pars = parms)

# Rhat histogram and plots
sum(stan_rhat(fit)$data > 1.01, na.rm = TRUE)

# DIC
deviance <- -2*rowSums(ms$log_lik)
(dic <- mean(deviance) + var(deviance)/2)

# WAIC
loo::waic(ms$log_lik)

# Bulk ESS
sum(s1$summary[,"n_eff"]<100*chains_val)
sum(s2$summary[,"n_eff"]<100*chains_val)
sum(s_pred$summary[,"n_eff"]<100*chains_val)

# Result summary
s1$summary 
s2$summary
s_pred$summary
