# Network meta-analysis for DTA in adults with NAFLD/NASH

## Overview
This repository contains a minimal set of R codes to reproduce our analyses shown in the paper:

## Contents
This repository includes:
- Stan codes for the [separate MA model](ma_fit_supp.stan), the [main NMA model](nma_fit_main_supp.stan), the [simpler model](nma_fit_sub_supp.stan), and the [inconsistency model](nma_fit_main_inconsistency_supp.stan),
- R codes to fit those models with particular options and values for the function arguments,
- the [raw data](NMA_data_supp.xlsx) that we used in the R codes and that includes study characteristics,
- the [data](NMA_data_inconsistency.xlsx) used for the inconsistency model,
- the [data](NMAIncludedStudies.csv) shows a list of studies with a flag of inclusion for MA and NMA, and
- a set of `renv` files to be used to restore the package environment.

