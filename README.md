# Elevated-Blood-Level
R scripts for statistical analysis of Elevated Blood Level

# Sociodemographic determinants and regional variation of pediatric elevated blood lead levels in Bangladesh: A survey-weighted analysis

This repository contains the complete, reproducible R script used for the statistical analysis and visualization in our manuscript examining elevated blood lead levels (EBLL; ≥ 5 µg/dL) among children aged 12–59 months in Bangladesh using data from the **Multiple Indicator Cluster Survey 7 (MICS7, 2025)**.

## 📁 Repository Structure
* `EBLL_analysis.R` — Main R script covering data merging, survey design specification, weighted prevalence estimation, LASSO variable selection, multivariable logistic/quasi-Poisson regression, marginal standardization, sensitivity analyses, and spatial mapping.
* `data/` *(not included due to data-sharing restrictions)* — Directory for raw MICS7 SPSS datasets (`ch.sav`, `hh.sav`) and Bangladesh administrative boundary shapefiles (`bgd_admin_boundaries.shp`).
* `outputs/` — Automatically generated directory where all publication tables (`.docx`/`.html`), summary statistics (`.csv`), high-resolution figures (`.png`), and `sessionInfo.txt` are saved.

## 📊 Data Availability
The microdata used in this study are publicly available upon registration and request from the **UNICEF MICS portal** ([mics.unicef.org](https://unicef.org)). Administrative boundary shapefiles for Bangladesh can be obtained from **HDX / OCHA** ([data.humdata.org](https://humdata.org)).

## 🛠️ Prerequisites & Required R Packages
The analysis was conducted in **R** (random seeds set at `123`, `2026`, and `2027` for full reproducibility). Install all required packages prior to running the script:

```r
install.packages(c(
  "haven", "dplyr", "survey", "gtsummary", "glmnet",
  "gt", "ggplot2", "sf", "marginaleffects", "car",
  "splines", "patchwork", "scales"
))
```

## 🚀 How to Reproduce the Analysis
1. Clone or download this repository to your local computer.
2. Place the MICS7 SPSS files (`ch.sav` and `hh.sav`) inside `data/raw/` and the shapefile (`bgd_admin_boundaries.shp` and associated files) inside `data/shapefiles/` (or update `in_dir`, `out_dir`, and `shp_path` at the top of the script).
3. Open `EBLL_analysis.R` in **RStudio** and run the script sequentially from **Step 0 to Step 9**.
4. Note: Set `regular_def <- TRUE` in Step 0 to reproduce the behavioural sensitivity analysis (regular contact ≥ 1/week vs. any contact in the past 3 months).

## 📝 Citation
If you use this code in your research, please cite our manuscript and this repository:
>*Sociodemographic determinants and regional variation of pediatric elevated blood lead levels in Bangladesh: A survey-weighted analysis*. GitHub repository: https://github.com/ehsanjnu/Elevated-Blood-Level

