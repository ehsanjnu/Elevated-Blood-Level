# Elevated-Blood-Level
R scripts for statistical analysis of Elevated Blood Level

# Sociodemographic Determinants and Regional Variation of Pediatric Elevated Blood Lead Levels in Bangladesh

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24239120.svg)](https://doi.org/10.5281/zenodo.24239120)

This repository contains the R scripts used for data preparation, variable derivation, exploratory weighted LASSO screening, and survey-weighted multivariable regression analyses for the study examining elevated blood lead levels (EBLL, ≥ 5 µg/dL) among children aged 12–59 months in Bangladesh using the 2024–2025 Multiple Indicator Cluster Survey (MICS7) dataset.

## Repository Contents

* `MICS7_EBLL_Revised_Analysis.R`: Complete, reproducible R script covering:
  1. Merging household (`hh.sav`) and child (`ch.sav`) datasets and codebook-based variable recoding.
  2. Complex survey design specification (`stratum`, primary sampling unit `HH1`, and biomarker-specific weight `chblweight`) using the `survey` package.
  3. Design-based prevalence estimation (with logit-based 95% confidence intervals) and Rao–Scott chi-square / survey-adjusted Kruskal–Wallis bivariate tests (Table 1).
  4. Exploratory survey-weighted LASSO logistic regression (`glmnet`) with PSU-clustered 10-fold cross-validation (`lambda.min` and `lambda.1se`) and 200 PSU-bootstrap stability resamples (Figure 1).
  5. Primary multivariable survey-weighted logistic regression (`svyglm`, quasibinomial family) on all 17 prespecified covariates, generalized variance inflation factor (GVIF) diagnostics, adjusted prevalence ratios (`quasipoisson`), and Rao–Wu bootstrap standardized prevalences (Table 2).
  6. Sensitivity analyses (alternative thresholds ≥ 3.5 µg/dL and > 5 µg/dL, continuous log-transformed BLL, natural cubic splines for soil lead, and urban/rural stratification).
  7. Division-level geospatial prevalence mapping using `sf` and `ggplot2` (Figure 2).

## Data Availability

The raw MICS7 datasets (`ch.sav` and `hh.sav`) are publicly available upon free registration from the [UNICEF Global MICS Program Repository](https://mics.unicef.org/surveys). Administrative boundary shapefiles (`bgd_admin1`) were obtained from the Humanitarian Data Exchange (HDX).

## Software Environment

All analyses were conducted in **R (version 4.3.1)**. Required packages:
* `haven` (v2.5.4)
* `dplyr` (v1.1.4)
* `survey` (v4.4.2)
* `glmnet` (v4.1.8)
* `gtsummary` (v2.4.0)
* `gt` (v1.1.0)
* `car` (v3.1.3)
* `splines` (v4.3.1)
* `marginaleffects`
* `sf` (v1.0.20)
* `ggplot2` (v3.5.2)
* `patchwork` (v1.3.0)

## Citation

If you use this code, please cite the archived repository:
> Ahmed, E., Methun, M. I. H., Alam, M. M., & Uddin, M. S. G. (2026). *R code for: Sociodemographic determinants and regional variation of pediatric elevated blood lead levels in Bangladesh: A survey-weighted analysis of the MICS7 biomarker survey*. Zenodo. https://doi.org/10.5281/zenodo.24239120
