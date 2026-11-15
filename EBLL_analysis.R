# ============================================================
# TITLE: Sociodemographic determinants and regional variation of pediatric elevated blood lead levels in Bangladesh: A survey-weighted analysis


# ============================================================
# STEP 0: PACKAGES & SETTINGS
# ============================================================
library(haven)
library(dplyr)
library(survey)
library(gtsummary)
library(glmnet)
library(gt)
library(ggplot2)
library(sf)
# install.packages(c("marginaleffects","car")) if missing
library(marginaleffects)
library(car)
library(splines)

options(survey.lonely.psu = "adjust")   # singleton PSU in a stratum
set.seed(123)

in_dir   <- "data/raw"          # Place ch.sav and hh.sav here
out_dir  <- "outputs"           # Generated tables and figures will be saved here
shp_path <- "data/shapefiles/bgd_admin_boundaries.shp"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Which definition of behavioural "Yes" for the LT frequency-scale variables?
#   FALSE -> any contact in last 3 months (codes 2,3,4)  [main]
#   TRUE  -> regular contact, at least once a week (codes 3,4) [sensitivity]
regular_def <- FALSE

save_gt <- function(tbl, name) {
  g <- as_gt(tbl)
  tryCatch(gt::gtsave(g, file.path(out_dir, paste0(name, ".docx"))),
           error = function(e) {
             message("docx save failed for ", name, " -> saving html instead")
             gt::gtsave(g, file.path(out_dir, paste0(name, ".html")))
           })
}

# ============================================================
# STEP 1: DATA LOAD & MERGE
# ============================================================
cat("\n--- STEP 1: LOADING AND MERGING DATA ---\n")
ch_data <- read_sav(file.path(in_dir, "ch.sav"))
hh_data <- read_sav(file.path(in_dir, "hh.sav"))

ch_data <- ch_data %>% mutate(HH1 = as.character(HH1), HH2 = as.character(HH2))
hh_data <- hh_data %>% mutate(HH1 = as.character(HH1), HH2 = as.character(HH2))

merged_data <- left_join(ch_data, hh_data, by = c("HH1", "HH2"), suffix = c("_ch", "_hh"))

get_col <- function(df, base_name) {
  variations <- c(
    base_name, toupper(base_name), tolower(base_name),
    paste0(base_name, "_hh"), paste0(toupper(base_name), "_hh"), paste0(tolower(base_name), "_hh"),
    paste0(base_name, "_ch"), paste0(toupper(base_name), "_ch"), paste0(tolower(base_name), "_ch")
  )
  for (v in variations) if (v %in% names(df)) return(df[[v]])
  matched <- names(df)[grepl(paste0("^", base_name, "$"), names(df), ignore.case = TRUE)]
  if (length(matched) > 0) return(df[[matched[1]]])
  warning("get_col(): no match for '", base_name, "' - returning all-NA column.")
  rep(NA_real_, nrow(df))
}

check_col <- function(x, base_name) {
  pct_na <- round(100 * mean(is.na(x)), 1)
  if (all(is.na(x))) warning("Column '", base_name, "' is 100% NA - NOT found correctly.")
  else if (pct_na > 50) cat("NOTE:", base_name, "is", pct_na, "% missing - verify.\n")
  x
}

# ------------------------------------------------------------
# CODEBOOK CHECKS (read the output once, then you can comment out)
# ------------------------------------------------------------
cat("\n=== Value labels of LT variables (CHECK coding) ===\n")
for (v in c("LT13","LT14","LT15","LT16","LT17A","LT17B","LT17C","LT17E",
            "LT18A","LT18B","LT18C","LT18E","LT20","LT27")) {
  x <- get_col(merged_data, v)
  cat("\n", v, ":\n"); print(attr(x, "labels")); print(table(x, useNA = "ifany"))
}
cat("\n=== Candidate STRATA variable names ===\n")
print(grep("strat", names(merged_data), ignore.case = TRUE, value = TRUE))
cat("\n=== Division (HH7) labels ===\n"); print(attr(get_col(merged_data, "HH7"), "labels"))

# ============================================================
# STEP 2: VARIABLE PREPARATION
# ============================================================
cat("\n--- STEP 2: PREPARING VARIABLES ---\n")

# LT variables: auto-detect frequency scale vs Yes/No coding from value labels
bin_lt <- function(x, regular = FALSE) {
  labs <- toupper(names(attr(x, "labels")))
  v <- as.numeric(x)
  if (length(labs) > 0 && any(grepl("NOT AT ALL", labs))) {        # frequency scale
    if (regular) {
      out <- case_when(v %in% 1:2 ~ "No", v %in% 3:4 ~ "Yes", TRUE ~ NA_character_)
    } else {
      out <- case_when(v == 1 ~ "No", v %in% 2:4 ~ "Yes", TRUE ~ NA_character_)
    }
  } else {                                                          # 1=Yes, 2=No  (CHECK)
    out <- case_when(v == 1 ~ "Yes", v == 2 ~ "No", TRUE ~ NA_character_)
  }
  out
}
# "Yes" if any component is Yes; "No" if none Yes and at least one No; else NA
any_yes <- function(...) {
  m <- cbind(...)
  case_when(rowSums(m == "Yes", na.rm = TRUE) > 0 ~ "Yes",
            rowSums(m == "No",  na.rm = TRUE) > 0 ~ "No",
            TRUE ~ NA_character_)
}
yn <- function(x) factor(x, levels = c("No", "Yes"))

# labelled -> factor, with "missing / DK / no response" levels set to NA
clean_f <- function(x) {
  f <- haven::as_factor(x)
  bad <- which(grepl("missing|dk|don't know|no response", as.character(f), ignore.case = TRUE))
  if (length(bad) > 0) f[bad] <- NA
  droplevels(f)
}

analytic_all <- merged_data %>%
  mutate(
    blood_lead_val = as.numeric(check_col(get_col(., "C_Lead"), "C_Lead")),
    # MICS7 codebook: C_Lead 999 = "Lab reject / Barcode mismatch", 0 = "Sample not taken".
    # Both are missing-data codes, not measurements; the "< 90" rule below excludes the 999 codes
    # (all 12 values >= 90 are exactly 999; the largest real measurement is ~85 ug/dL).
    Y_elevated_lead = case_when(
      blood_lead_val >= 5 & blood_lead_val < 90 ~ 1,     # MAIN: >= 5 ug/dL
      blood_lead_val <  5 ~ 0,
      TRUE ~ NA_real_),
    Y_gt5  = case_when(blood_lead_val > 5   & blood_lead_val < 90 ~ 1,
                       blood_lead_val <= 5  ~ 0, TRUE ~ NA_real_),   # sensitivity
    Y_ge35 = case_when(blood_lead_val >= 3.5 & blood_lead_val < 90 ~ 1,
                       blood_lead_val <  3.5 ~ 0, TRUE ~ NA_real_),  # CDC reference value
    Y_elevated_lead_factor = factor(ifelse(Y_elevated_lead == 1, "Elevated (>=5)", "Normal (<5)"),
                                    levels = c("Elevated (>=5)", "Normal (<5)")),
    Test_status = as.numeric(check_col(get_col(., "LT27"), "LT27")),

    # SOCIO-DEMOGRAPHICS
    W_age    = as.numeric(check_col(get_col(., "CAGE"), "CAGE")),
    W_sex    = factor(ifelse(get_col(., "HL4") == 2, "Female", "Male")),
    W_rural  = factor(ifelse(get_col(., "HH6") == 2, "Rural", "Urban"), levels = c("Rural", "Urban")),
    W_wealth = clean_f(get_col(., "WINDEX5")),
    W_meduc  = clean_f(get_col(., "MELEVEL")),     # "Missing/DK" -> NA (n=3)

    # ENVIRONMENT & HOUSING
    Env_floor_raw = as.numeric(get_col(., "HC4")),
    Env_floor = factor(case_when(
      Env_floor_raw %in% c(11, 12) ~ "Unimproved (Earth/Sand/Dung)",
      Env_floor_raw >= 21 ~ "Improved (Wood/Cement/Tiles/Etc)",
      TRUE ~ NA_character_)),
    Env_paint_raw = as.numeric(get_col(., "HC6A")),
    Env_paint = factor(case_when(
      Env_paint_raw == 1 ~ "Flaking/Peeling",
      Env_paint_raw == 2 ~ "Good condition",
      Env_paint_raw == 3 ~ "No paint/Other",
      TRUE ~ NA_character_)),
    Env_water_raw = as.numeric(get_col(., "WS1")),
    Env_water = factor(case_when(
      Env_water_raw %in% c(11, 12, 13, 14, 21, 22, 31, 41, 51, 91) ~ "Improved",
      Env_water_raw %in% c(32, 42, 61, 71, 72, 81, 92, 93, 95, 96, 99) ~ "Unimproved (Surface water/Unprotected)",
      TRUE ~ NA_character_)),
    Env_toilet_raw = as.numeric(get_col(., "WS11")),
    Env_toilet = factor(case_when(
      Env_toilet_raw %in% c(11, 12, 13, 14, 15, 18, 21, 22, 31) ~ "Improved",
      Env_toilet_raw %in% c(23, 24, 32, 41, 51, 95, 96) ~ "Unimproved (Open pit/Hanging/None)",
      TRUE ~ NA_character_)),
    Soil_Lead   = as.numeric(check_col(get_col(., "S_Lead"), "S_Lead")),
    Soil_Lead10 = Soil_Lead / 10,                  # OR reported per 10 ppm

    # BEHAVIOURAL (codebook-based recoding; see bin_lt)
    Beh_pica    = yn(bin_lt(get_col(., "LT20"), regular_def)),
    Beh_amulet  = yn(bin_lt(get_col(., "LT13"), regular_def)),
    Beh_sindoor = yn(bin_lt(get_col(., "LT14"), regular_def)),
    Beh_surma   = yn(any_yes(bin_lt(get_col(., "LT15"), regular_def),
                             bin_lt(get_col(., "LT16"), regular_def))),

    # HAZARDS (residence direction / person in regular contact with child)
    # LT17C/LT18C (repair/disassemble batteries) added.  CHECK labels: 1=Yes,2=No assumed
    Haz_battery = yn(any_yes(bin_lt(get_col(., "LT17A")), bin_lt(get_col(., "LT17B")),
                             bin_lt(get_col(., "LT17C")), bin_lt(get_col(., "LT18A")),
                             bin_lt(get_col(., "LT18B")), bin_lt(get_col(., "LT18C")))),
    Haz_ewaste  = yn(any_yes(bin_lt(get_col(., "LT17E")), bin_lt(get_col(., "LT18E")))),

    # DESIGN VARIABLES
    chblweight = as.numeric(check_col(get_col(., "chblweight"), "chblweight")),
    cluster    = get_col(., "HH1"),
    stratum    = get_col(., "stratum")             # CHECK name (see strata printout above)
  )

# Division: readable names, Barishal as reference
div_chr <- as.character(haven::as_factor(get_col(analytic_all, "HH7")))
div_chr <- tools::toTitleCase(tolower(div_chr))
div_chr <- dplyr::recode(div_chr, "Barisal" = "Barishal", "Chittagong" = "Chattogram")
analytic_all$W_division <- factor(div_chr)
if ("Barishal" %in% levels(analytic_all$W_division))
  analytic_all$W_division <- relevel(analytic_all$W_division, ref = "Barishal")
cat("\nDivision levels:", paste(levels(analytic_all$W_division), collapse = ", "), "\n")

if (all(is.na(analytic_all$stratum)))
  stop("Strata variable not found. Set `stratum = get_col(., \"<name>\")` using the printout above.")

# ------------------------------------------------------------
# Participant flow (reviewer: flow diagram + exclusion counts)
# ------------------------------------------------------------
n_all    <- nrow(analytic_all)
d1 <- analytic_all %>% filter(Test_status == 1)
d2 <- d1 %>% filter(!is.na(blood_lead_val))
d3 <- d2 %>% filter(blood_lead_val < 90)
d4 <- d3 %>% filter(W_age >= 12 & W_age <= 59)
d5 <- d4 %>% filter(!is.na(chblweight))
flow <- data.frame(
  Step = c("Under-5 children in ch file (questionnaire records)",
           "Blood lead test completed (LT27 == 1)",
           "Non-missing blood lead value",
           "Lab result coded 999 (lab reject / barcode mismatch) excluded",
           "Aged 12-59 months",
           "Non-missing chblweight (FINAL analytic sample)"),
  N = c(n_all, nrow(d1), nrow(d2), nrow(d3), nrow(d4), nrow(d5)))
flow$Excluded <- c(NA, -diff(flow$N))
print(flow)
write.csv(flow, file.path(out_dir, "Participant_flow.csv"), row.names = FALSE)
cat("\nLT27 (blood test status) among children with a recorded status:\n")
print(table(analytic_all$Test_status, useNA = "ifany"))
cat("Biomarker non-response among those with a recorded LT27:",
    round(100 * (1 - sum(analytic_all$Test_status == 1, na.rm = TRUE) /
                   sum(!is.na(analytic_all$Test_status))), 1), "%\n")
cat("Zero/very low BLL values (check detection limit):",
    sum(d5$blood_lead_val <= 0), "zeros;", sum(d5$blood_lead_val == 5), "exactly = 5.0\n")

analytic_data <- d5 %>% filter(!is.na(Y_elevated_lead))
cat("\n*** FINAL unweighted analytic n =", nrow(analytic_data), "***\n")
if (nrow(analytic_data) == 0) stop("0 rows: check C_Lead / LT27 / CAGE / chblweight names.")

# ------------------------------------------------------------
# Candidate predictors, missingness, complete-case flag
# ------------------------------------------------------------
cand_vars <- c("W_age", "W_sex", "W_rural", "W_wealth", "W_division", "W_meduc",
               "Env_floor", "Env_paint", "Env_water", "Env_toilet", "Soil_Lead10",
               "Beh_pica", "Beh_amulet", "Beh_sindoor", "Beh_surma",
               "Haz_battery", "Haz_ewaste")

miss_tab <- data.frame(
  Variable = cand_vars,
  N_missing = sapply(cand_vars, function(v) sum(is.na(analytic_data[[v]]))),
  Pct_missing = round(sapply(cand_vars, function(v) 100 * mean(is.na(analytic_data[[v]]))), 2))
print(miss_tab, row.names = FALSE)
write.csv(miss_tab, file.path(out_dir, "Missingness_by_predictor.csv"), row.names = FALSE)

analytic_data$complete_case <- factor(complete.cases(analytic_data[, cand_vars]),
                                      levels = c(TRUE, FALSE), labels = c("Included (complete case)", "Excluded (missing covariate)"))
cat("Complete cases for modelling:", sum(analytic_data$complete_case == "Included (complete case)"),
    "of", nrow(analytic_data), "\n")

# ============================================================
# STEP 3: SURVEY DESIGN, PREVALENCE (with 95% CI)
# ============================================================
cat("\n--- STEP 3: SURVEY DESIGN & PREVALENCE ---\n")
survey_design <- svydesign(ids = ~cluster, strata = ~stratum, weights = ~chblweight,
                           data = analytic_data, nest = TRUE)
cat("Design df:", degf(survey_design), "\n")
cat("Sum of weights:", round(sum(analytic_data$chblweight), 1),
    "| unweighted n:", nrow(analytic_data), "\n")
cat("UNWEIGHTED elevated:", sum(analytic_data$Y_elevated_lead == 1),
    "(", round(100 * mean(analytic_data$Y_elevated_lead == 1), 1), "% )\n")

prev <- svyciprop(~I(Y_elevated_lead == 1), survey_design, method = "logit", df = degf(survey_design))
ci <- attr(prev, "ci")
cat(sprintf("WEIGHTED prevalence (>=5 ug/dL): %.1f%% (95%% CI %.1f-%.1f)\n",
            100 * coef(prev), 100 * ci[1], 100 * ci[2]))
# weights are normalised to sum ~ n, so svymean warns "sample size greater than population
# size" (harmless). deff = "replace" gives the design effect despite that.
print(svymean(~Y_elevated_lead, survey_design, deff = "replace"))
cat("Compare with MICS7 official figure (reviewer quotes 38.34%). If far apart, check age range/exclusions.\n")

# Division-level weighted prevalence with CI (used for Figure 2 and Table S)
div_prev <- svyby(~I(Y_elevated_lead == 1), ~W_division, survey_design, svyciprop,
                  vartype = "ci", method = "logit", df = degf(survey_design))
names(div_prev) <- c("Division", "Prev", "CI_low", "CI_high")
div_prev[, 2:4] <- round(100 * div_prev[, 2:4], 1)
print(div_prev, row.names = FALSE)
write.csv(div_prev, file.path(out_dir, "Division_weighted_prevalence.csv"), row.names = FALSE)

# ============================================================
# STEP 4: TABLE 1 (column %) + TABLE 1b (row % = prevalence with CI)
# ============================================================
cat("\n--- STEP 4: TABLES 1 / 1b ---\n")
t1_vars <- c("W_age", "W_sex", "W_rural", "W_wealth", "W_division", "W_meduc", "Env_floor",
             "Env_paint", "Env_water", "Env_toilet", "Soil_Lead", "Beh_pica", "Beh_amulet",
             "Beh_sindoor", "Beh_surma", "Haz_battery", "Haz_ewaste")
t1_labels <- list(
  W_age ~ "Age (months)", W_sex ~ "Child sex", W_rural ~ "Place of residence",
  W_wealth ~ "Wealth quintile", W_division ~ "Administrative division",
  W_meduc ~ "Mother's education", Env_floor ~ "Main floor material",
  Env_paint ~ "Wall paint condition", Env_water ~ "Main drinking water source",
  Env_toilet ~ "Toilet facility", Soil_Lead ~ "Soil lead concentration (ppm)",
  Beh_pica ~ "Pica (ingesting soil/non-food items), last 3 months",
  Beh_amulet ~ "Wore metallic amulet, last 3 months",
  Beh_sindoor ~ "Contact with sindoor, last 3 months",
  Beh_surma ~ "Use of / contact with surma, last 3 months",
  Haz_battery ~ "Battery manufacture/recycling/repair (nearby or contact person)",
  Haz_ewaste ~ "E-waste recycling (nearby or contact person)")

table1 <- survey_design %>%
  tbl_svysummary(
    by = Y_elevated_lead_factor, include = all_of(t1_vars),
    statistic = list(all_categorical() ~ "{n_unweighted} ({p}%)",
                     all_continuous()  ~ "{median} ({p25}, {p75})"),
    label = t1_labels, missing = "ifany", missing_text = "Missing") %>%
  add_overall() %>% add_p() %>% bold_labels() %>% bold_p(t = 0.05) %>%
  modify_caption("**Table 1. Characteristics by EBLL status (n = unweighted count; % = survey-weighted COLUMN percentage; continuous: weighted median (IQR)). Rao-Scott test.**")
save_gt(table1, "Table1_Column_percent")

# Table 1b: weighted prevalence of EBLL within each category, with design-based CI
cat_vars <- setdiff(t1_vars, c("W_age", "Soil_Lead"))
prev_by <- function(v) {
  keep <- !is.na(survey_design$variables[[v]])
  d <- survey_design[keep, ]
  r <- svyby(~I(Y_elevated_lead == 1), as.formula(paste0("~", v)), d, svyciprop,
             vartype = "ci", method = "logit", df = degf(survey_design))
  data.frame(Variable = v, Category = as.character(r[[1]]),
             Unweighted_n = as.integer(table(d$variables[[v]])[as.character(r[[1]])]),
             Prevalence_pct = round(100 * r[[2]], 1),
             CI_low = round(100 * r[[3]], 1), CI_high = round(100 * r[[4]], 1))
}
table1b <- bind_rows(lapply(cat_vars, prev_by))
write.csv(table1b, file.path(out_dir, "Table1b_Weighted_prevalence_by_category.csv"), row.names = FALSE)
print(table1b, row.names = FALSE)

# Included vs excluded (complete-case) comparison
tbl_cc <- survey_design %>%
  tbl_svysummary(by = complete_case,
                 include = c(W_sex, W_rural, W_wealth, W_division, Y_elevated_lead_factor),
                 statistic = list(all_categorical() ~ "{n_unweighted} ({p}%)"), missing = "ifany") %>%
  add_p() %>% modify_caption("**Table S1. Included vs excluded (missing covariate) children**")
save_gt(tbl_cc, "TableS1_Included_vs_Excluded")

# ============================================================
# STEP 5: WEIGHTED LASSO (exploratory screening)
# ============================================================
cat("\n--- STEP 5: WEIGHTED LASSO ---\n")
lasso_data <- analytic_data %>%
  select(Y_elevated_lead, chblweight, cluster, stratum, all_of(cand_vars)) %>% na.omit()
cat("LASSO n (complete cases):", nrow(lasso_data), "| events:", sum(lasso_data$Y_elevated_lead), "\n")

# reference coding for ALL factors (with "- 1" the first factor, W_sex, got both a Male and a
# Female column, which was confusing in Figure 1c)
x_matrix <- model.matrix(~ ., data = lasso_data[, cand_vars])[, -1]
y_vector <- lasso_data$Y_elevated_lead
w_vector <- lasso_data$chblweight / mean(lasso_data$chblweight)    # normalised weights
cat("Parameters (columns) entering LASSO:", ncol(x_matrix), "\n")

# PSU-based folds: all children of a PSU go to the same fold
cl_id  <- as.integer(factor(paste(lasso_data$stratum, lasso_data$cluster)))
foldid <- sample(rep(1:10, length.out = max(cl_id)))[cl_id]

lasso_model <- cv.glmnet(x_matrix, y_vector, family = "binomial", alpha = 1,
                         weights = w_vector, foldid = foldid, type.measure = "deviance")
cat("lambda.min =", signif(lasso_model$lambda.min, 4),
    "| lambda.1se =", signif(lasso_model$lambda.1se, 4), "\n")

cand_sorted <- cand_vars[order(-nchar(cand_vars))]
col_base <- sapply(colnames(x_matrix), function(cn) cand_sorted[startsWith(cn, cand_sorted)][1])
get_sel <- function(s) {
  cf <- as.matrix(coef(lasso_model, s = s))
  nz <- rownames(cf)[cf[, 1] != 0 & rownames(cf) != "(Intercept)"]
  list(coef = nz, base = unique(col_base[nz]), cf = cf)
}
sel_min <- get_sel("lambda.min"); sel_1se <- get_sel("lambda.1se")
cat("Non-zero coefficients: lambda.min =", length(sel_min$coef),
    "| lambda.1se =", length(sel_1se$coef), "\n")
cat("Variables kept (lambda.min):", paste(sel_min$base, collapse = ", "), "\n")
cat("Variables kept (lambda.1se):", paste(sel_1se$base, collapse = ", "), "\n")

# Figure 1a: CV curve; Figure 1b: coefficient paths; Figure 1c: ranked coefficients
png(file.path(out_dir, "Figure1a_LASSO_CV_curve.png"), width = 2400, height = 1800, res = 300)
plot(lasso_model); title("10-fold (PSU-based) CV deviance", line = 2.5); dev.off()
png(file.path(out_dir, "Figure1b_LASSO_coefficient_paths.png"), width = 2400, height = 1800, res = 300)
plot(lasso_model$glmnet.fit, xvar = "lambda", label = FALSE)
abline(v = log(c(lasso_model$lambda.min, lasso_model$lambda.1se)), lty = 2); dev.off()

pretty_names <- function(x) {
  x <- sub("^W_division", "Division: ", x);  x <- sub("^W_meduc", "Maternal education: ", x)
  x <- sub("^W_wealth", "Wealth quintile: ", x); x <- sub("^W_sex", "Sex: ", x)
  x <- sub("^W_rural", "Residence: ", x);    x <- sub("^W_age$", "Age (months)", x)
  x <- sub("^Env_floor", "Floor: ", x);      x <- sub("^Env_paint", "Wall paint: ", x)
  x <- sub("^Env_water", "Drinking water: ", x); x <- sub("^Env_toilet", "Toilet: ", x)
  x <- sub("^Soil_Lead10$", "Soil lead (per 10 ppm)", x)
  x <- sub("^Beh_pica", "Pica: ", x);        x <- sub("^Beh_amulet", "Metallic amulet: ", x)
  x <- sub("^Beh_sindoor", "Sindoor contact: ", x); x <- sub("^Beh_surma", "Surma use/contact: ", x)
  x <- sub("^Haz_battery", "Battery-related activities: ", x); x <- sub("^Haz_ewaste", "E-waste recycling: ", x)
  x
}
lasso_df <- data.frame(Predictor = sel_min$coef, Coefficient = sel_min$cf[sel_min$coef, 1]) %>%
  mutate(Direction = ifelse(Coefficient > 0, "Higher odds (+)", "Lower odds (-)"),
         Predictor = reorder(pretty_names(Predictor), abs(Coefficient)))
fig1c <- ggplot(lasso_df, aes(Coefficient, Predictor, fill = Direction)) +
  geom_col(alpha = 0.8, color = "black") +
  scale_fill_manual(values = c("Higher odds (+)" = "#d73027", "Lower odds (-)" = "#4575b4")) +
  theme_minimal(base_size = 13) +
  labs(title = "LASSO-selected predictors (lambda.min), weighted",
       subtitle = "Coefficients are on the log-odds scale; factor levels vs reference category",
       x = "LASSO coefficient (log-odds)", y = "", fill = "") +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())
ggsave(file.path(out_dir, "Figure1c_LASSO_ranked_coefficients.png"), fig1c, width = 10, height = 8, dpi = 300)

# ---- Combined Figure 1 (panels A-C): A over B on the left, C on the right ----
# install.packages("patchwork")   # if not installed
library(patchwork)
cvdf <- data.frame(loglam = log(lasso_model$lambda), cvm = lasso_model$cvm,
                   lo = lasso_model$cvlo, hi = lasso_model$cvup)
lam_lines <- log(c(lasso_model$lambda.min, lasso_model$lambda.1se))
pA <- ggplot(cvdf, aes(loglam, cvm)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), colour = "grey70", width = 0) +
  geom_point(colour = "red3", size = 1.4) +
  geom_vline(xintercept = lam_lines, linetype = 3) +
  annotate("text", x = lam_lines, y = max(cvdf$hi), label = c("λ_min", "λ_1se"),
           hjust = -0.1, vjust = 1, size = 3.2) +
  labs(x = "log(λ)", y = "Binomial deviance") + theme_minimal(base_size = 11)
beta <- as.matrix(lasso_model$glmnet.fit$beta)
pathdf <- data.frame(loglam = rep(log(lasso_model$glmnet.fit$lambda), each = nrow(beta)),
                     coef = as.vector(beta), term = rep(rownames(beta), times = ncol(beta)))
pB <- ggplot(pathdf, aes(loglam, coef, group = term)) +
  geom_line(colour = "grey35", linewidth = 0.4) +
  geom_vline(xintercept = lam_lines, linetype = 2) +
  labs(x = "log(λ)", y = "Coefficient (log-odds)") + theme_minimal(base_size = 11)
# same font size as A and B (fig1c itself uses base_size 13)
pC <- fig1c + labs(title = NULL, subtitle = NULL) + theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())
fig1_combined <- ((pA / pB) | pC) + plot_layout(widths = c(1, 1.3)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 14))
ggsave(file.path(out_dir, "Figure1_LASSO_combined.png"), fig1_combined,
       width = 13, height = 9, dpi = 300, bg = "white")

# Selection stability: PSU bootstrap at fixed lambda.min
B <- 200
psu_key <- paste(lasso_data$stratum, lasso_data$cluster)
psu_list <- split(seq_len(nrow(lasso_data)), psu_key)
sel_count <- setNames(rep(0, length(cand_vars)), cand_vars)
for (b in seq_len(B)) {
  idx <- unlist(psu_list[sample(length(psu_list), replace = TRUE)], use.names = FALSE)
  fit <- glmnet(x_matrix[idx, ], y_vector[idx], family = "binomial", alpha = 1,
                weights = w_vector[idx], lambda = lasso_model$lambda.min)
  cf <- as.matrix(coef(fit)); nz <- rownames(cf)[cf[, 1] != 0 & rownames(cf) != "(Intercept)"]
  sel_count[unique(col_base[nz])] <- sel_count[unique(col_base[nz])] + 1
}
stability <- data.frame(Variable = names(sel_count), Selection_pct = round(100 * sel_count / B, 1)) %>%
  arrange(desc(Selection_pct))
print(stability, row.names = FALSE)
write.csv(stability, file.path(out_dir, "LASSO_selection_stability.csv"), row.names = FALSE)

selected_base_vars <- sel_min$base
if (length(selected_base_vars) == 0) selected_base_vars <- cand_vars

# ============================================================
# STEP 6: MULTIVARIABLE SURVEY LOGISTIC MODELS (Table 2)
# ============================================================
cat("\n--- STEP 6: MULTIVARIABLE MODELS ---\n")
lab_list <- list(Soil_Lead10 ~ "Soil lead (per 10 ppm)", W_age ~ "Age (months)")
fit_svy <- function(vars, y = "Y_elevated_lead", des = survey_design, fam = quasibinomial()) {
  svyglm(as.formula(paste(y, "~", paste(vars, collapse = " + "))), design = des, family = fam)
}

# A) LASSO-selected (EXPLORATORY: post-selection CI/p-values are optimistic)
# NOTE: lambda.min kept all 17 candidates (identical to model B), so model A uses the
# sparser lambda.1se selection, which makes the comparison informative.
model_A <- fit_svy(sel_1se$base)
# B) Full candidate model (a-priori; 17 candidates vs thousands of events - no need to penalise)
model_B <- fit_svy(cand_vars)

tA <- tbl_regression(model_A, exponentiate = TRUE, label = lab_list) %>% bold_labels() %>% bold_p(t = 0.05)
tB <- tbl_regression(model_B, exponentiate = TRUE, label = lab_list) %>% bold_labels() %>% bold_p(t = 0.05)
table2 <- tbl_merge(list(tA, tB), tab_spanner = c("**A. LASSO lambda.1se (exploratory)**", "**B. Full candidate model (primary)**")) %>%
  modify_caption("**Table 2. Adjusted odds ratios (survey-weighted quasibinomial). Associations, not causal effects.**")
save_gt(table2, "Table2_aOR_models_A_B")

# Multicollinearity: GVIF on an unpenalised weighted GLM (VIF is not a model-validity test)
glm_vif <- glm(as.formula(paste("Y_elevated_lead ~", paste(cand_vars, collapse = " + "))),
               data = analytic_data, weights = chblweight / mean(chblweight), family = quasibinomial())
vif_out <- as.data.frame(car::vif(glm_vif)); print(vif_out)
write.csv(vif_out, file.path(out_dir, "GVIF_full_model.csv"))

# ============================================================
# STEP 7: PREVALENCE RATIOS & MARGINAL (ADJUSTED) PREVALENCES
# ============================================================
cat("\n--- STEP 7: PR AND MARGINAL ESTIMATES (outcome is common, ~38%) ---\n")
model_PR <- fit_svy(cand_vars, fam = quasipoisson(link = "log"))   # exp(coef) = prevalence ratio
tPR <- tbl_regression(model_PR, exponentiate = TRUE, label = lab_list) %>% bold_labels() %>% bold_p(t = 0.05) %>%
  modify_header(estimate ~ "**PR**") %>% modify_abbreviation("PR = Prevalence ratio") %>%
  modify_caption("**Table S2. Adjusted prevalence ratios (survey quasi-Poisson, log link)**")
save_gt(tPR, "TableS2_Prevalence_ratios")

safe_me <- function(expr, name) tryCatch({
  r <- as.data.frame(expr); write.csv(r, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE); r
}, error = function(e) { message("marginaleffects skipped: ", name, " (replaced by bootstrap block below)"); NULL })

# newdata = the complete cases used to fit model_B, WITH the weight column
nd_cc <- analytic_data[complete.cases(analytic_data[, cand_vars]), ]
for (v in c("W_wealth", "W_meduc", "W_division", "W_rural")) {
  safe_me(avg_predictions(model_B, variables = v, newdata = nd_cc, wts = "chblweight"),
          paste0("Adjusted_prev_", v))
  safe_me(avg_comparisons(model_B, variables = v, newdata = nd_cc, wts = "chblweight",
                          comparison = "differenceavg"), paste0("Adjusted_PD_", v))
  safe_me(avg_comparisons(model_B, variables = v, newdata = nd_cc, wts = "chblweight",
                          comparison = "ratioavg"), paste0("Adjusted_PR_", v))
}
cat("Delta-method adjusted prevalences saved where marginaleffects worked.\n")

# ---- Design-based standardised prevalence, prevalence difference (PD) and ratio (PR) ----
# Rao-Wu bootstrap replicate weights; independent of marginaleffects. For each replicate the
# full model is refitted, every child is set to each level of the variable, predictions are
# averaged with the replicate weights, and PD/PR are taken versus the first (reference) level.
std_vars <- c("W_wealth", "W_meduc", "W_division", "W_rural")
des_cc   <- survey_design[complete.cases(survey_design$variables[, cand_vars]), ]
set.seed(2026)   # fixed seed: bootstrap CIs are now reproducible
rep_cc   <- as.svrepdesign(des_cc, type = "bootstrap", replicates = 500)
form_cc  <- as.formula(paste("Y_elevated_lead ~", paste(cand_vars, collapse = " + ")))
keys <- unlist(lapply(std_vars, function(v) paste(v, levels(des_cc$variables[[v]]), sep = "|")))

theta_std <- function(w, data) {
  # glm() looks up `weights` in `data` and the formula's environment (global), not in this
  # function, so the replicate weights must be stored as a column of `data`
  data$wrep <- w
  fit <- glm(form_cc, data = data, weights = wrep, family = quasibinomial())
  out <- numeric(0)
  for (v in std_vars) {
    lv <- levels(data[[v]])
    for (l in lv) {
      d <- data; d[[v]] <- factor(l, levels = lv)
      out <- c(out, weighted.mean(predict(fit, newdata = d, type = "response"), w))
    }
  }
  out
}
cat("Running bootstrap standardisation (500 replicates; takes several minutes)...\n")
rep_res <- withReplicates(rep_cc, theta_std, return.replicates = TRUE)
est  <- as.numeric(rep_res$theta); names(est) <- keys
reps <- as.matrix(rep_res$replicates); colnames(reps) <- keys

res_list <- list()
for (v in std_vars) {
  k <- grep(paste0("^", v, "\\|"), keys, value = TRUE); ref <- k[1]
  for (kk in k) {
    is_ref <- identical(kk, ref)
    q <- function(x) as.numeric(quantile(x, c(0.025, 0.975), na.rm = TRUE))
    pd_rep <- reps[, kk] - reps[, ref]; pr_rep <- reps[, kk] / reps[, ref]
    res_list[[kk]] <- data.frame(
      Variable = v, Level = sub(".*\\|", "", kk),
      Adj_prev_pct = 100 * est[kk],
      Prev_low = 100 * q(reps[, kk])[1], Prev_high = 100 * q(reps[, kk])[2],
      PD_pct_points = 100 * (est[kk] - est[ref]),
      PD_low = if (is_ref) NA else 100 * q(pd_rep)[1], PD_high = if (is_ref) NA else 100 * q(pd_rep)[2],
      PR = est[kk] / est[ref],
      PR_low = if (is_ref) NA else q(pr_rep)[1], PR_high = if (is_ref) NA else q(pr_rep)[2],
      row.names = NULL)
  }
}
std_tab <- bind_rows(res_list) %>% mutate(across(where(is.numeric), ~ round(.x, 2)))
print(std_tab, row.names = FALSE)
write.csv(std_tab, file.path(out_dir, "Adjusted_standardised_prev_PD_PR.csv"), row.names = FALSE)

# ============================================================
# STEP 8: SENSITIVITY ANALYSES
# ============================================================
cat("\n--- STEP 8: SENSITIVITY ANALYSES ---\n")
# (a) alternative thresholds, full candidate model
m_main <- model_B
m_gt5  <- fit_svy(cand_vars, y = "Y_gt5")
m_ge35 <- fit_svy(cand_vars, y = "Y_ge35")
sens_tbl <- tbl_merge(
  list(tbl_regression(m_main, exponentiate = TRUE, label = lab_list),
       tbl_regression(m_gt5,  exponentiate = TRUE, label = lab_list),
       tbl_regression(m_ge35, exponentiate = TRUE, label = lab_list)),
  tab_spanner = c("**>= 5 (main)**", "**> 5**", "**>= 3.5 (CDC ref.)**")) %>%
  modify_caption("**Table S3. Sensitivity of aORs to the blood-lead threshold**")
save_gt(sens_tbl, "TableS3_Threshold_sensitivity")
cat("Weighted prevalence >5:",
    round(100 * coef(svymean(~Y_gt5, survey_design)), 1), "% | >=3.5:",
    round(100 * coef(svymean(~Y_ge35, survey_design)), 1), "%\n")

# (b) continuous outcome: log(BLL); exp(coef) = ratio of geometric means
des_pos <- subset(survey_design, blood_lead_val > 0)
des_pos <- update(des_pos, log_bll = log(blood_lead_val))
m_cont <- svyglm(as.formula(paste("log_bll ~", paste(cand_vars, collapse = " + "))), design = des_pos)
t_cont <- tbl_regression(m_cont, exponentiate = TRUE, label = lab_list) %>% bold_p(t = 0.05) %>%
  modify_caption("**Table S4. Ratio of geometric mean blood lead (survey-weighted linear model on log scale)**")
save_gt(t_cont, "TableS4_Continuous_logBLL")

# (c) soil-lead functional form: linear vs natural spline (df = 3)
no_soil <- setdiff(cand_vars, "Soil_Lead10")
m_lin <- svyglm(as.formula(paste("Y_elevated_lead ~", paste(c(no_soil, "Soil_Lead10"), collapse = " + "))),
                design = survey_design, family = quasibinomial())
m_spl <- svyglm(as.formula(paste("Y_elevated_lead ~", paste(c(no_soil, "ns(Soil_Lead10, 3)"), collapse = " + "))),
                design = survey_design, family = quasibinomial())
cat("\nSoil lead non-linearity (Wald test of spline terms):\n")
print(regTermTest(m_spl, ~ns(Soil_Lead10, 3)))
iqr <- IQR(analytic_data$Soil_Lead, na.rm = TRUE)
b <- coef(m_lin)["Soil_Lead10"]; se <- sqrt(vcov(m_lin)["Soil_Lead10", "Soil_Lead10"])
cat(sprintf("Soil lead aOR per IQR (%.1f ppm): %.2f (%.2f-%.2f)\n", iqr,
            exp(b * iqr / 10), exp((b - 1.96 * se) * iqr / 10), exp((b + 1.96 * se) * iqr / 10)))

# (c1) soil-lead diagnostics. OR 1.27 per 10 ppm looks large next to PR 1.03 and a
# geometric-mean ratio of 1.04 -> check the tail / outliers and use quartiles
cat("\nSoil lead distribution (ppm):\n"); print(summary(analytic_data$Soil_Lead))
print(quantile(analytic_data$Soil_Lead, c(.01, .05, .5, .95, .99, .999), na.rm = TRUE))
cat("Median soil lead by division:\n"); print(tapply(analytic_data$Soil_Lead, analytic_data$W_division, median, na.rm = TRUE))
sq <- cut(survey_design$variables$Soil_Lead,
          quantile(survey_design$variables$Soil_Lead, 0:4 / 4, na.rm = TRUE),
          include.lowest = TRUE, labels = paste0("Q", 1:4))
des_sq <- update(survey_design, Soil_q = sq)
m_sq <- svyglm(as.formula(paste("Y_elevated_lead ~", paste(c(no_soil, "Soil_q"), collapse = " + "))),
               design = des_sq, family = quasibinomial())
save_gt(tbl_regression(m_sq, exponentiate = TRUE, include = "Soil_q") %>% bold_p(t = 0.05) %>%
          modify_caption("**Table S6. Soil lead quartiles (Q1 = reference), adjusted OR**"),
        "TableS6_Soil_lead_quartiles")

# (c1b) soil lead has a very long right tail (max 674 ppm, 99th pct ~58): try log2 scale and
#       capping at the 99th percentile
p99 <- as.numeric(quantile(analytic_data$Soil_Lead, 0.99, na.rm = TRUE))
cat("99th percentile of soil lead:", p99, "ppm;",
    sum(analytic_data$Soil_Lead > p99, na.rm = TRUE), "children above\n")
des_sl <- update(survey_design, Soil_log2 = log2(Soil_Lead), Soil_w10 = pmin(Soil_Lead, p99) / 10)
m_log2 <- svyglm(as.formula(paste("Y_elevated_lead ~", paste(c(no_soil, "Soil_log2"), collapse = " + "))),
                 design = des_sl, family = quasibinomial())
m_win  <- svyglm(as.formula(paste("Y_elevated_lead ~", paste(c(no_soil, "Soil_w10"), collapse = " + "))),
                 design = des_sl, family = quasibinomial())
t_sl <- tbl_stack(list(
  tbl_regression(m_log2, exponentiate = TRUE, include = "Soil_log2",
                 label = list(Soil_log2 ~ "Soil lead (per doubling)")),
  tbl_regression(m_win, exponentiate = TRUE, include = "Soil_w10",
                 label = list(Soil_w10 ~ "Soil lead capped at 99th pct (per 10 ppm)")))) %>%
  modify_caption("**Table S6b. Soil lead: alternative scalings (adjusted OR)**")
print(as_tibble(t_sl))
save_gt(t_sl, "TableS6b_Soil_lead_scalings")

# (g) full model WITHOUT soil lead: uses the larger sample (soil lead drives most missingness
#     and excluded children differ from included ones - see Table S1)
m_nosoil <- fit_svy(setdiff(cand_vars, "Soil_Lead10"))
cat("n in model without soil lead:", nobs(m_nosoil), "\n")
save_gt(tbl_regression(m_nosoil, exponentiate = TRUE, label = list(W_age ~ "Age (months)")) %>% bold_p(t = 0.05) %>%
          modify_caption("**Table S7. Full model without soil lead (larger analytic sample)**"),
        "TableS7_Model_without_soil_lead")

# (c2) soil-lead dose-response: adjusted (standardised) prevalence across the observed range
soil_grid <- as.numeric(quantile(nd_cc$Soil_Lead10, probs = seq(0.02, 0.98, length.out = 25)))
soil_curve <- tryCatch(
  as.data.frame(avg_predictions(m_spl, variables = list(Soil_Lead10 = soil_grid),
                                newdata = nd_cc, wts = "chblweight")),
  error = function(e) { message("soil curve failed: ", conditionMessage(e)); NULL })
if (FALSE && !is.null(soil_curve)) {   # DISABLED: marginaleffects + ns() gave impossible values (~98%)
  soil_curve$Soil_ppm <- soil_curve$Soil_Lead10 * 10
  write.csv(soil_curve, file.path(out_dir, "Soil_lead_adjusted_curve.csv"), row.names = FALSE)
  fig_soil <- ggplot(soil_curve, aes(Soil_ppm, estimate)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) + geom_line(linewidth = 1) +
    scale_y_continuous(labels = scales::percent) + theme_minimal(base_size = 13) +
    labs(x = "Soil lead (ppm), 2nd-98th percentile", y = "Adjusted prevalence of EBLL",
         title = "Adjusted association between soil lead and EBLL (natural spline, df = 3)")
  ggsave(file.path(out_dir, "FigureS_Soil_lead_curve.png"), fig_soil, width = 8, height = 5.5, dpi = 300)
}

# (c2-fixed) soil-lead dose-response curve.
# The earlier version let marginaleffects recompute ns() on counterfactual data, which produced
# impossible values (~98% prevalence). Here the spline basis is computed ONCE (fixed knots),
# stored as columns, and standardised predictions + CIs come from Rao-Wu bootstrap replicates.
cc      <- survey_design[complete.cases(survey_design$variables[, cand_vars]), ]
dat_cc  <- cc$variables
sl      <- dat_cc$Soil_Lead10
ns_obj  <- splines::ns(sl, df = 3)
sp_cols <- paste0("sp", 1:3)
dat_cc[, sp_cols] <- ns_obj[, 1:3]
grid    <- as.numeric(quantile(sl, probs = seq(0.02, 0.98, length.out = 25)))
ns_grid <- predict(ns_obj, grid)[, 1:3]
form_sp <- as.formula(paste("Y_elevated_lead ~", paste(c(no_soil, sp_cols), collapse = " + ")))
des_sp  <- svydesign(ids = ~cluster, strata = ~stratum, weights = ~chblweight, data = dat_cc, nest = TRUE)
set.seed(2027)
rep_sp  <- as.svrepdesign(des_sp, type = "bootstrap", replicates = 500)

theta_soil <- function(w, data) {
  data$wrep <- w
  fit <- glm(form_sp, data = data, weights = wrep, family = quasibinomial())
  sapply(seq_along(grid), function(i) {
    d <- data
    d[, sp_cols] <- matrix(ns_grid[i, ], nrow(d), 3, byrow = TRUE)
    weighted.mean(predict(fit, newdata = d, type = "response"), w)
  })
}
cat("Soil-lead curve: bootstrap standardisation (500 replicates)...\n")
res_soil <- withReplicates(rep_sp, theta_soil, return.replicates = TRUE)
soil_curve2 <- data.frame(
  Soil_ppm = grid * 10,
  estimate = as.numeric(res_soil$theta),
  conf.low  = apply(res_soil$replicates, 2, quantile, probs = 0.025, na.rm = TRUE),
  conf.high = apply(res_soil$replicates, 2, quantile, probs = 0.975, na.rm = TRUE))
print(round(soil_curve2, 3), row.names = FALSE)
write.csv(soil_curve2, file.path(out_dir, "Soil_lead_adjusted_curve_fixed.csv"), row.names = FALSE)
fig_soil2 <- ggplot(soil_curve2, aes(Soil_ppm, estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) + geom_line(linewidth = 1) +
  geom_vline(xintercept = quantile(analytic_data$Soil_Lead, 0.75, na.rm = TRUE), linetype = 3) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) + theme_minimal(base_size = 13) +
  labs(x = "Soil lead (ppm), 2nd-98th percentile (dotted line = 75th percentile)",
       y = "Adjusted prevalence of EBLL",
       title = "Adjusted association between soil lead and EBLL (natural spline, df = 3)")
ggsave(file.path(out_dir, "FigureS_Soil_lead_curve_fixed.png"), fig_soil2, width = 8, height = 5.5, dpi = 300)

# (d) behaviour definition: re-run script with regular_def <- TRUE (top of file)
# (e) multiple imputation: only if covariate missingness is substantial (see Missingness_by_predictor.csv)

# (f) Surma interaction and urban/rural stratified models (full candidate set)
base_ns <- setdiff(cand_vars, "Beh_surma")
f_int <- function(mod) as.formula(paste("Y_elevated_lead ~",
                   paste(setdiff(base_ns, mod), collapse = " + "), "+ Beh_surma *", mod))
cat("\n--- Wald: Surma x Division ---\n")
print(regTermTest(svyglm(f_int("W_division"), design = survey_design, family = quasibinomial()), ~Beh_surma:W_division))
cat("\n--- Wald: Surma x Wealth ---\n")
print(regTermTest(svyglm(f_int("W_wealth"), design = survey_design, family = quasibinomial()), ~Beh_surma:W_wealth))

# Formal tests behind the stratified results (wealth / education gradients differ by residence?)
f_rint <- function(term) as.formula(paste("Y_elevated_lead ~",
                paste(setdiff(cand_vars, term), collapse = " + "), "+", term, "* W_rural"))
cat("\n--- Wald: Wealth x Residence ---\n")
print(regTermTest(svyglm(f_rint("W_wealth"), design = survey_design, family = quasibinomial()), ~W_wealth:W_rural))
cat("\n--- Wald: Maternal education x Residence ---\n")
print(regTermTest(svyglm(f_rint("W_meduc"), design = survey_design, family = quasibinomial()), ~W_meduc:W_rural))

cat("\nDivision x Residence cells:\n"); print(table(analytic_data$W_division, analytic_data$W_rural))
strat_vars <- setdiff(cand_vars, "W_rural")
m_rural <- fit_svy(strat_vars, des = subset(survey_design, W_rural == "Rural"))
m_urban <- fit_svy(strat_vars, des = subset(survey_design, W_rural == "Urban"))
tbl_strat <- tbl_merge(
  list(tbl_regression(m_rural, exponentiate = TRUE, label = lab_list) %>% bold_p(t = 0.05),
       tbl_regression(m_urban, exponentiate = TRUE, label = lab_list) %>% bold_p(t = 0.05)),
  tab_spanner = c("**Rural**", "**Urban**")) %>%
  modify_caption("**Table S5. Stratified survey-weighted logistic regression by residence**")
save_gt(tbl_strat, "TableS5_Stratified_Urban_Rural")

# ============================================================
# STEP 9: FIGURE 2 - division-level weighted prevalence (REAL data)
# ============================================================
cat("\n--- STEP 9: MAP ---\n")
bd_map <- st_read(shp_path, layer = "bgd_admin1")     # CHECK layer name / shapefile
print(names(bd_map))

# Match shapefile names to survey names (handles common spelling variants)
bd_map$adm1_clean <- dplyr::recode(tools::toTitleCase(tolower(bd_map$adm1_name)),
                            "Barisal" = "Barishal", "Chittagong" = "Chattogram")
ebll_summary <- div_prev %>% rename(Division_Name = Division, EBLL_Prevalence = Prev)
unmatched <- setdiff(ebll_summary$Division_Name, bd_map$adm1_clean)
if (length(unmatched) > 0) warning("Divisions not matched to shapefile: ", paste(unmatched, collapse = ", "))

map_data <- left_join(bd_map, ebll_summary, by = c("adm1_clean" = "Division_Name"))

ebll_map <- ggplot(map_data) +
  geom_sf(aes(fill = EBLL_Prevalence), color = "black", linewidth = 0.3) +
  geom_sf_text(aes(label = paste0(adm1_clean, "\n", sprintf("%.1f", EBLL_Prevalence), "%")),
               size = 3.2, color = "black", fontface = "bold") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Weighted EBLL\nprevalence (%)") +
  theme_void() +
  labs(title = "Division-level prevalence of elevated blood lead (>= 5 ug/dL)",
       subtitle = "Survey-weighted prevalence among children 12-59 months, MICS7 2024-25",
       caption = "Source: Bangladesh MICS7 2024-25, weighted by the biomarker weight (chblweight)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        plot.caption = element_text(hjust = 0.5, size = 9))
print(ebll_map)
ggsave(file.path(out_dir, "Figure_2_EBLL_Map.png"), ebll_map,
       width = 8, height = 8, units = "in", dpi = 300, bg = "white")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))  # reproducibility record
cat("\n=== ALL ANALYSES COMPLETED. Outputs in:", out_dir, "===\n")
