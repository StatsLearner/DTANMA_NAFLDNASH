# ######################################################################
# - This is a R script to fit network meta-analysis using a separate stan file.
# - Created by Tetsuro Oda
# ######################################################################

# ---------------------------------------------------
# Preparation
# ---------------------------------------------------

library(tidyverse)
library(rstan)
source("create_datasets_supp.R")

fibrosis_grade <- "F2"
prior_type <- 1

# Read data
d <- readxl::read_excel("NMA_data_supp.xlsx",sheet = 1)

# Include studies whose cutoffs are not defined to get 90 % of sen or spe
d2<- 
  d %>%
  filter(Flag_Threshold_sesp_90 == 0) 

# Create a dataset for network meta-analysis
df <- create_dataset_NMA(d2, fibrosis_grade)

# Create a dataset for stan
data_stan <- list(N = nrow(df),
                  J = 2,
                  S = length(unique(df$`Rayyan ID`)),
                  T = length(unique(df$Test)),
                  Test = df$Test_Num,
                  Study = df$Study_Num,
                  Threshold = df$Threshold_Num,
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
model <- stan_model(file = "nma_fit_main_supp.stan")

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
parms <- c("muSe","muSp","sigSe","sigSp","tau0_se",
           "tau1_se","tau0_sp","tau1_sp")

ms <- rstan::extract(fit)
s_params <- summary(fit, pars = parms)
s_p <- summary(fit, pars = c("sepool","sppool"))
s_pred <- summary(fit, pars = c("sesp_pred"))

s_params$summary
s_p$summary
s_pred$summary

s_sediff <- summary(fit, pars = c("se_diff"))
s_spdiff <- summary(fit, pars = c("sp_diff"))

s_sepreddiff <- summary(fit, pars = c("sepred_diff"))
s_sppreddiff <- summary(fit, pars = c("sppred_diff"))

# ---------------------------------------------------
# Check MCMC samples
# ---------------------------------------------------

pairs(fit, pars = c("muSe[1]", "muSp[1]","sigSe","sigSp"))
pairs(fit, pars = c("tau0_se","tau1_se","tau0_sp","tau1_sp"))
pairs(fit, pars = c("gamma0_se[1]","gamma0_sp[1]",
                    "gamma0_se[2]","gamma0_sp[2]",
                    "gamma1_se[1,2]","gamma1_sp[2,3]"))
pairs(fit, pars = "sepool")
pairs(fit, pars = "sppool")

# Autocorrelations
stan_ac(fit, pars = c("gamma0_se[1]","gamma0_sp[1]",
                      "gamma0_se[2]","gamma0_sp[2]",
                      "gamma1_se[1,2]","gamma1_sp[2,3]"))

stan_ac(fit, pars = parms)

# Density plots
stan_dens(fit, pars = c("gamma0_se[1]","gamma0_sp[1]",
                        "gamma0_se[2]","gamma0_sp[2]",
                        "gamma1_se[1,2]","gamma1_sp[2,3]"), 
          separate_chains = T)

# Trace plots
stan_trace(fit, pars = c("gamma0_se[1]","gamma0_sp[1]",
                         "gamma0_se[2]","gamma0_sp[2]",
                         "gamma1_se[1,2]","gamma1_sp[2,3]"))

# Stan diagnosis
stan_diag(fit)

# Rhat 
sum(stan_rhat(fit)$data > 1.01, na.rm = TRUE)

# DIC
deviance <- -2*rowSums(ms$log_lik)
(dic <- mean(deviance) + var(deviance)/2)

# LooIC
loo2 <- 
  loo::loo(ms$log_lik,
           cores = 1,
           r_eff = loo::relative_eff(
             exp(ms$log_lik),
             chain_id = rep(1:chains_val, 
                            each = (iter_val - warmup_val)/thin_val),
             cores = 1
           )
  )

print(loo2)
plot(loo2)
loo::pareto_k_influence_values(loo2)
loo::psis_n_eff_values(loo2)

# WAIC
loo::waic(ms$log_lik)
