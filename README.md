# Network meta-analysis for DTA in adults with NAFLD

## Overview
This repository contains a minimal set of R codes to reproduce our analyses shown in the paper: amaguchi, R., Oda, T. & Nagashima, K. Comparison of the diagnostic accuracy of shear wave elastography with transient elastography in adult nonalcoholic fatty liver disease: a systematic review and network meta-analysis of diagnostic test accuracy. Abdom Radiol (2024). https://doi.org/10.1007/s00261-024-04546-8

## Contents
This repository includes:
- Stan codes for the [separate MA model](ma_fit_supp.stan), the [main NMA model](nma_fit_main_supp.stan), the [simpler model](nma_fit_sub_supp.stan), and the [inconsistency model](nma_fit_main_inconsistency_supp.stan),
- R codes to fit those models with particular options and values for the function arguments to reproduce the results written in the paper,
- the [raw data](NMA_data_supp.xlsx) that we used in the R codes and that includes study characteristics,
- the [data](NMA_data_inconsistency.xlsx) used for the inconsistency model,
- the [data](NMAIncludedStudies.csv) shows a list of studies with a flag of inclusion for MA and NMA, and
- a set of `renv` files to be used to restore the package environment.

