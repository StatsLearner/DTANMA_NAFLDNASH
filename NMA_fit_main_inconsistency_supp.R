# ######################################################################
# - This is a R script to fit a NMA inconsistency model.
# - Created by Kengo Nagashima
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
d <- readxl::read_excel("NMA_data_inconsistency.xlsx", sheet = 1)

# Include studies whose cutoffs are not defined to get 90 % of sens or spe
d2<- 
  d %>%
  filter(Flag_Threshold_sesp_90 == 0) 

# Create a dataset for network meta-analysis
df <- create_dataset_NMA(d2, fibrosis_grade)

# Further add flags for inconsistenct parameters
d_inconsistency <-
  df %>%
  mutate(
    fTE = case_when(Test == "TE" ~ 1, TRUE ~ 0),
    fpSWE = case_when(Test == "pSWE" ~ 1, TRUE ~ 0),
    ftSWE = case_when(Test == "2D-SWE" ~ 1, TRUE ~ 0),
    fMRE = case_when(Test == "MRE" ~ 1, TRUE ~ 0)
  ) %>%
  group_by(`Rayyan ID`) %>%
  summarize(TE = max(fTE), pSWE = max(fpSWE),
            tSWE = max(ftSWE), MRE = max(fMRE)) %>%
  mutate(
    id = `Rayyan ID`,
    wtype = case_when(
      TE == 1 & pSWE == 1 & tSWE == 1 & MRE == 0 ~ 1,
      TE == 1 & pSWE == 1 & tSWE == 0 & MRE == 0 ~ 2,
      TE == 0 & pSWE == 1 & tSWE == 1 & MRE == 0 ~ 3,
      TE == 0 & pSWE == 1 & tSWE == 0 & MRE == 1 ~ 4
    )
  )

df <- 
  df %>% mutate(id = `Rayyan ID`)

df <- 
  left_join(df, d_inconsistency, by = "id") %>%
  mutate(
    wtype = case_when(
      Test == "pSWE" ~ 0,
      TRUE ~ wtype
    ))

# Create a dataset for stan
data_stan <- list(N = nrow(df),
                  J = 2,
                  S = length(unique(df$id)),
                  T = length(unique(df$Test)),
                  Test = df$Test_Num,
                  Study = df$Study_Num,
                  Threshold = df$Threshold_Num,
                  TP = df$TP*10^3,
                  FP = df$FP*10^3,
                  FN = df$FN*10^3,
                  TN = df$TN*10^3,
                  prior_type = prior_type,
                  W = df$wtype
)

# ---------------------------------------------------
# Fitting
# ---------------------------------------------------

# Compile a model
model <- stan_model("nma_fit_main_inconsistency_supp.stan")

# Set some options
options(mc.cores = parallelly::availableCores() - 1)
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
parms <- c("muSe", "muSp", "sigSe", "sigSp", 
           "tau0_se", "tau1_se", "tau0_sp", "tau1_sp",
           "omegaSe[1]", "omegaSe[2]", "omegaSp[1]", "omegaSp[2]",
           "sepool[1]", "sepool[2]", "sepool[3]", "sepool[4]",
           "sppool[1]", "sppool[2]", "sppool[3]", "sppool[4]")

ms <- rstan::extract(fit)
s_params <- summary(fit, pars = parms)
s_params$summary[,c("mean","2.5%","50%","97.5%")]

parms <- c(
  "d_Se_TE_pSWE_direct", "d_Sp_TE_pSWE_direct",
  "d_Se_2dSWE_pSWE_direct", "d_Sp_2dSWE_pSWE_direct",
  "d_Se_TE_pSWE_indirect", "d_Sp_TE_pSWE_indirect",
  "d_Se_MRE_pSWE_indirect", "d_Sp_MRE_pSWE_indirect",
  "d_Se_2dSWE_pSWE_indirect", "d_Sp_2dSWE_pSWE_indirect"
)

s_params <- summary(fit, pars = parms)
s_params$summary[,c("mean","2.5%","50%","97.5%")]

save.image(file=paste0("stanfit_NMA_main_supp_inconsistency",
                       fibrosis_grade,"_",prior_type,".RData"))

# ---------------------------------------------------
# Check MCMC samples
# ---------------------------------------------------
parms <- c(parms, c("gamma0_se[1]", "gamma0_sp[1]",
                    "gamma0_se[2]", "gamma0_sp[2]",
                    "gamma1_se[1,2]", "gamma1_sp[2,3]"))

pairs(fit, pars = c("muSe[1]", "muSp[1]", "sigSe", "sigSp"))
pairs(fit, pars = c("tau0_se", "tau1_se", "tau0_sp", "tau1_sp"))
pairs(fit, pars = c("gamma0_se[1]", "gamma0_sp[1]",
                    "gamma0_se[2]", "gamma0_sp[2]",
                    "gamma1_se[1,2]", "gamma1_sp[2,3]"))
pairs(fit, pars = c("omegaSe[1]", "omegaSe[2]", 
                    "omegaSp[1]", "omegaSp[2]"))

# Autocorrelations
stan_ac(fit, pars = parms)

# Density plots
stan_dens(fit, pars = parms, separate_chains = T)

# Trace plots
stan_trace(fit, pars = parms)

# Stan diagnosis
stan_diag(fit)

# Rhat 
sum(stan_rhat(fit)$data > 1.01, na.rm = TRUE)
