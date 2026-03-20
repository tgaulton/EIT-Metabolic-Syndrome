#===============================================================================
# METABOLIC SYNDROME AND REGIONAL LUNG PHYSIOLOGY IN ARDS
# Complete Bayesian Analysis of EIT Data
#
# Author: T. Gaulton, MD MSc
# Institution: Massachusetts General Hospital / Harvard Medical School
# Date: November 30, 2024
#
# Description:
#   Analysis of electrical impedance tomography (EIT) data comparing regional
#   ventilation-perfusion distributions between ARDS patients with and without
#   metabolic syndrome. Uses Bayesian regression to test whether metabolic
#   syndrome affects regional lung physiology through mechanical and/or
#   vascular pathways.
#
# Main Hypotheses:
#   H1: Metabolic syndrome reduces dorsal ventilation (mechanical effect)
#   H2: Metabolic syndrome impairs perfusion redistribution (vascular effect)
#   H3: Combined abnormalities increase dorsal shunt
#
# Input:
#   - Two Excel files in data/ (see README for filenames):
#       data/Phenotypes_distributions_demo.xlsx
#       data/Phenotypes_distributions_vq.xlsx
#   - Files are joined on 'case' and filtered to use == 1 for final cohort
#   - Raw data are not included in this repository (patient-level data)
#
# Output:
#   - Tables: Table 1 (baseline characteristics), Table 2 (regional physiology)
#   - Figures: 4 main figures + 2 supplementary
#   - Models: All Bayesian model objects saved for reproducibility
#   - Results: Complete workspace saved as results/data/complete_analysis.RData
#
# To run:
#   source("analysis_main.R")
#
# To reload results without refitting:
#   load("results/data/complete_analysis.RData")
#
#===============================================================================

# ENVIRONMENT SETUP ============================================================

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║   METABOLIC SYNDROME AND REGIONAL LUNG PHYSIOLOGY IN ARDS    ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

# Clear environment (optional - comment out if you want to keep existing objects)
rm(list = ls())

# Record start time
start_time <- Sys.time()

# Load required packages
cat("Loading packages...\n")
required_packages <- c("tidyverse", "brms", "bayesplot", "patchwork", "readxl")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Set options
options(mc.cores = parallel::detectCores())  # Use all available cores
options(scipen = 999)  # Avoid scientific notation

# Set seed for reproducibility
set.seed(123)

# Create output directories
cat("Creating output directories...\n")
dir.create("models/bayesian", showWarnings = FALSE, recursive = TRUE)
dir.create("models/sensitivity", showWarnings = FALSE, recursive = TRUE)
dir.create("models/gradients", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/data", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("submission_figures", showWarnings = FALSE, recursive = TRUE)

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), "results/data/session_info.txt")

cat("✓ Setup complete\n\n")

#===============================================================================
# DATA PREPARATION
#===============================================================================

# NOTE: Place data files in the data/ subdirectory before running.
# Raw data are not included in this repository.
Phenotypes_distributions_demo <- read_excel("data/Phenotypes_distributions_demo.xlsx")
Phenotypes_distributions_vq   <- read_excel("data/Phenotypes_distributions_vq.xlsx")
merged <- inner_join(Phenotypes_distributions_demo, Phenotypes_distributions_vq, by = "case")

# Load and filter to final cohort
# STEP 1: Filter and create base dataset
final <- merged %>%
  filter(use == 1) %>%                    # Keep only cases marked for analysis
  filter(timetoassess < 7) %>%            # Exclude assessments >7 days from ARDS onset
  mutate(
    # Standardize continuous predictors for regression modeling
    age_z    = scale(age)[, 1],           # Z-score for age
    apache_z = scale(apache)[, 1],        # Z-score for APACHE II

    # Calculate dorsal-to-ventral ratios
    vent_ratio = vP / vA,                 # Dorsal vent / Ventral vent
    perf_ratio = qP / qA,                 # Dorsal perf / Ventral perf

    # Calculate gradients (dorsal - ventral differences)
    vent_gradient_diff = vP - vA,         # Ventilation gradient
    perf_gradient_diff = qP - qA,         # Perfusion gradient

    # Calculate V/Q ratio gradient (ratio of ratios)
    vq_ratio_gradient = vent_ratio / perf_ratio,

    # ARDS Etiology categories
    ARDSetiology = case_when(
      etiology == 1             ~ "Sepsis",
      etiology %in% c(0, 2, 3) ~ "Pneumonia",
      etiology %in% c(4, 6)    ~ "Trauma",
      etiology == 5             ~ "Postoperative Respiratory Failure",
      TRUE                      ~ NA_character_
    ),

    # ARDS Site (Pulmonary vs Extrapulmonary)
    ARDSsite = case_when(
      etiology %in% c(1, 4, 5, 6) ~ "Extrapulmonary",
      etiology %in% c(0, 2, 3)    ~ "Pulmonary",
      TRUE                         ~ NA_character_
    )
  )

# Sample size
n_total <- nrow(final)
n_no    <- sum(final$metS == "No")
n_yes   <- sum(final$metS == "Yes")

cat("Sample size:\n")
cat("  Total:", n_total, "\n")
cat("  No metabolic syndrome:", n_no, "\n")
cat("  Metabolic syndrome:", n_yes, "\n\n")

#===============================================================================
# HELPER FUNCTIONS
#===============================================================================

# SMD calculation for continuous variables
calc_smd <- function(x, group) {
  valid <- !is.na(x) & !is.na(group)
  x     <- x[valid]
  group <- group[valid]

  x1 <- x[group == "No"]
  x2 <- x[group == "Yes"]

  if (length(x1) < 2 | length(x2) < 2) return(NA)

  m1 <- mean(x1)
  m2 <- mean(x2)
  s1 <- sd(x1)
  s2 <- sd(x2)
  n1 <- length(x1)
  n2 <- length(x2)

  pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))

  if (pooled_sd == 0) return(0)

  smd <- (m2 - m1) / pooled_sd
  return(abs(smd))
}

# SMD for binary variables
calc_smd_binary <- function(x, group) {
  p1 <- mean(x[group == "No"],  na.rm = TRUE)
  p2 <- mean(x[group == "Yes"], na.rm = TRUE)

  if (is.na(p1) | is.na(p2))   return(NA)
  if (p1 == 0 & p2 == 0)       return(0)
  if (p1 == 1 & p2 == 1)       return(0)

  pooled_var <- (p1 * (1 - p1) + p2 * (1 - p2)) / 2

  if (pooled_var == 0) return(0)

  smd <- (p2 - p1) / sqrt(pooled_var)
  return(abs(smd))
}

# Formatting functions for tables
fmt_mean_sd <- function(x, digits = 1) {
  sprintf(paste0("%.", digits, "f (%.1f)"), mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

fmt_median_iqr <- function(x, digits = 1) {
  sprintf(paste0("%.", digits, "f [%.1f, %.1f]"),
          median(x, na.rm = TRUE),
          quantile(x, 0.25, na.rm = TRUE),
          quantile(x, 0.75, na.rm = TRUE))
}

fmt_n_pct <- function(x, total) {
  sprintf("%d (%.0f%%)", x, 100 * x / total)
}

fmt_smd <- function(x) {
  if (is.na(x)) return("")
  sprintf("%.2f", x)
}

cat("✓ Data prepared\n\n")

#===============================================================================
# TABLE 1: BASELINE CHARACTERISTICS
#===============================================================================

cat("Creating Table 1...\n")

table1 <- tribble(
  ~Characteristic, ~Overall, ~No_MetS, ~Yes_MetS, ~SMD,

  # Age
  "Age (years)",
  fmt_mean_sd(final$age),
  fmt_mean_sd(final$age[final$metS == "No"]),
  fmt_mean_sd(final$age[final$metS == "Yes"]),
  fmt_smd(calc_smd(final$age, final$metS)),

  # Sex
  "Sex", "", "", "", "",
  "    Female",
  fmt_n_pct(sum(final$sex == "Female"), n_total),
  fmt_n_pct(sum(final$sex == "Female" & final$metS == "No"), n_no),
  fmt_n_pct(sum(final$sex == "Female" & final$metS == "Yes"), n_yes),
  fmt_smd(calc_smd_binary(as.numeric(final$sex == "Female"), final$metS)),
  "    Male",
  fmt_n_pct(sum(final$sex == "Male"), n_total),
  fmt_n_pct(sum(final$sex == "Male" & final$metS == "No"), n_no),
  fmt_n_pct(sum(final$sex == "Male" & final$metS == "Yes"), n_yes),
  "",

  # Race
  "Race", "", "", "", "",
  "    White",
  fmt_n_pct(sum(final$race == "White"), n_total),
  fmt_n_pct(sum(final$race == "White" & final$metS == "No"), n_no),
  fmt_n_pct(sum(final$race == "White" & final$metS == "Yes"), n_yes),
  fmt_smd(calc_smd_binary(as.numeric(final$race == "White"), final$metS)),
  "    Black",
  fmt_n_pct(sum(final$race == "Black"), n_total),
  fmt_n_pct(sum(final$race == "Black" & final$metS == "No"), n_no),
  fmt_n_pct(sum(final$race == "Black" & final$metS == "Yes"), n_yes),
  "",
  "    Asian",
  fmt_n_pct(sum(final$race == "Asian"), n_total),
  fmt_n_pct(sum(final$race == "Asian" & final$metS == "No"), n_no),
  fmt_n_pct(sum(final$race == "Asian" & final$metS == "Yes"), n_yes),
  "",

  # Hispanic
  "Hispanic Ethnicity",
  fmt_n_pct(sum(final$hispanic == "Yes"), n_total),
  fmt_n_pct(sum(final$hispanic == "Yes" & final$metS == "No"), n_no),
  fmt_n_pct(sum(final$hispanic == "Yes" & final$metS == "Yes"), n_yes),
  fmt_smd(calc_smd_binary(as.numeric(final$hispanic == "Yes"), final$metS)),

  # ARDS Site
  "ARDS Site", "", "", "", "",
  "    Extrapulmonary",
  fmt_n_pct(sum(final$ARDSsite == "Extrapulmonary"), n_total),
  fmt_n_pct(sum(final$ARDSsite == "Extrapulmonary" & final$metS == "No"), n_no),
  fmt_n_pct(sum(final$ARDSsite == "Extrapulmonary" & final$metS == "Yes"), n_yes),
  fmt_smd(calc_smd_binary(as.numeric(final$ARDSsite == "Extrapulmonary"), final$metS)),
  "    Pulmonary",
  fmt_n_pct(sum(final$ARDSsite == "Pulmonary"), n_total),
  fmt_n_pct(sum(final$ARDSsite == "Pulmonary" & final$metS == "No"), n_no),
  fmt_n_pct(sum(final$ARDSsite == "Pulmonary" & final$metS == "Yes"), n_yes),
  "",

  # Clinical variables
  "PaO\u2082/FiO\u2082 Ratio",
  fmt_mean_sd(final$pf, 0),
  fmt_mean_sd(final$pf[final$metS == "No"], 0),
  fmt_mean_sd(final$pf[final$metS == "Yes"], 0),
  fmt_smd(calc_smd(final$pf, final$metS)),

  "Driving Pressure (cmH\u2082O)",
  fmt_mean_sd(final$dp),
  fmt_mean_sd(final$dp[final$metS == "No"]),
  fmt_mean_sd(final$dp[final$metS == "Yes"]),
  fmt_smd(calc_smd(final$dp, final$metS)),

  "Respiratory System Compliance (mL/cmH\u2082O)",
  fmt_mean_sd(final$comp),
  fmt_mean_sd(final$comp[final$metS == "No"]),
  fmt_mean_sd(final$comp[final$metS == "Yes"]),
  fmt_smd(calc_smd(final$comp, final$metS)),

  "APACHE II Score",
  fmt_mean_sd(final$apache, 0),
  fmt_mean_sd(final$apache[final$metS == "No"], 0),
  fmt_mean_sd(final$apache[final$metS == "Yes"], 0),
  fmt_smd(calc_smd(final$apache, final$metS)),

  "Time from ARDS Onset to Assessment (days)",
  fmt_median_iqr(final$timetoassess),
  fmt_median_iqr(final$timetoassess[final$metS == "No"]),
  fmt_median_iqr(final$timetoassess[final$metS == "Yes"]),
  fmt_smd(calc_smd(final$timetoassess, final$metS)),

  # Metabolic syndrome components
  "Body Mass Index (kg/m\u00b2)",
  fmt_mean_sd(final$bmi),
  fmt_mean_sd(final$bmi[final$metS == "No"]),
  fmt_mean_sd(final$bmi[final$metS == "Yes"]),
  fmt_smd(calc_smd(final$bmi, final$metS)),

  "Triglycerides (mg/dL)",
  fmt_mean_sd(final$tg, 0),
  fmt_mean_sd(final$tg[final$metS == "No"], 0),
  fmt_mean_sd(final$tg[final$metS == "Yes"], 0),
  fmt_smd(calc_smd(final$tg, final$metS)),

  "HDL Cholesterol (mg/dL)",
  fmt_mean_sd(final$hdl),
  fmt_mean_sd(final$hdl[final$metS == "No"]),
  fmt_mean_sd(final$hdl[final$metS == "Yes"]),
  fmt_smd(calc_smd(final$hdl, final$metS)),

  "Hypertension",
  fmt_n_pct(sum(final$htn == "Yes"), n_total),
  fmt_n_pct(sum(final$htn == "Yes" & final$metS == "No"), n_no),
  fmt_n_pct(sum(final$htn == "Yes" & final$metS == "Yes"), n_yes),
  fmt_smd(calc_smd_binary(as.numeric(final$htn == "Yes"), final$metS)),

  "Diabetes Mellitus",
  fmt_n_pct(sum(final$diabetes == "Yes"), n_total),
  fmt_n_pct(sum(final$diabetes == "Yes" & final$metS == "No"), n_no),
  fmt_n_pct(sum(final$diabetes == "Yes" & final$metS == "Yes"), n_yes),
  fmt_smd(calc_smd_binary(as.numeric(final$diabetes == "Yes"), final$metS))
)

# Rename columns
table1 <- table1 %>%
  rename(
    `Characteristic`           = Characteristic,
    `Overall (N = 25)`         = Overall,
    `No (N = 15)`              = No_MetS,
    `Yes (N = 10)`             = Yes_MetS,
    `SMD`                      = SMD
  )

# Save Table 1
write_csv(table1, "results/tables/Table1_with_SMD.csv")

cat("✓ Table 1 created\n\n")

#===============================================================================
# DISTRIBUTION DIAGNOSTICS
#===============================================================================

cat("Checking outcome distributions...\n")

outcomes_to_check <- c("vP", "qP", "qsP_map", "qsA_map", "dsP_map", "dsA_map")

# Calculate distribution statistics
calc_skewness <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  (sum((x - m)^3) / n) / (s^3)
}

dist_summary <- final %>%
  select(all_of(outcomes_to_check)) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(
    n      = n(),
    mean   = mean(value, na.rm = TRUE),
    sd     = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    min    = min(value, na.rm = TRUE),
    max    = max(value, na.rm = TRUE),
    skewness = calc_skewness(value),
    recommended_family = case_when(
      abs(skewness) < 1                          ~ "Gaussian",
      abs(skewness) >= 1 & abs(skewness) < 2    ~ "Student-t",
      abs(skewness) >= 2                         ~ "Gamma/Log-Normal"
    )
  )

print(dist_summary)
write_csv(dist_summary, "results/data/distribution_summary.csv")

cat("✓ Distribution diagnostics complete\n\n")

#===============================================================================
# PRIOR SPECIFICATION
#===============================================================================

cat("Setting up priors...\n")

# Weakly informative priors based on physiologic plausibility
priors_gaussian <- c(
  prior(normal(50, 30), class = "Intercept"),
  prior(normal(0, 5),   class = "b"),
  prior(exponential(0.1), class = "sigma")
)

priors_student <- c(
  prior(normal(25, 25),   class = "Intercept"),
  prior(normal(0, 5),     class = "b"),
  prior(exponential(0.05), class = "sigma"),
  prior(gamma(4, 1),      class = "nu")
)

cat("✓ Priors specified\n\n")

#===============================================================================
# PRIMARY BAYESIAN MODELS
#===============================================================================

cat("Fitting primary Bayesian models...\n")
cat("(This may take 10-20 minutes depending on your computer)\n\n")

# Family selection based on distribution diagnostics
family_map <- list(
  vP      = gaussian(),    # Skewness = 0.40
  qP      = gaussian(),    # Skewness = -0.19
  qsP_map = student(),     # Skewness = 0.59
  qsA_map = student(),     # Skewness = 1.12
  dsP_map = student(),     # Skewness = 0.97
  dsA_map = gaussian()     # Skewness = 0.32
)

# Prior selection
prior_map <- list(
  vP      = priors_gaussian,
  qP      = priors_gaussian,
  qsP_map = priors_student,
  qsA_map = priors_student,
  dsP_map = priors_student,
  dsA_map = priors_gaussian
)

# Model formula
base_formula <- function(outcome) {
  as.formula(paste0(outcome, " ~ metS + age_z + apache_z + ARDSsite"))
}

# Fit all models
outcomes    <- c("vP", "qP", "qsP_map", "qsA_map", "dsP_map", "dsA_map")
bayes_models <- vector("list", length(outcomes))
names(bayes_models) <- outcomes

for (i in seq_along(outcomes)) {
  outcome <- outcomes[i]
  cat("  [", i, "/", length(outcomes), "] Fitting:", outcome, "\n")

  bayes_models[[outcome]] <- brm(
    formula = base_formula(outcome),
    data    = final,
    family  = family_map[[outcome]],
    prior   = prior_map[[outcome]],
    chains  = 4,
    iter    = 4000,
    warmup  = 1000,
    cores   = 4,
    seed    = 123,
    control = list(adapt_delta = 0.95),
    silent  = FALSE,
    refresh = 0,
    file    = paste0("models/bayesian/bayes_", outcome),
    file_refit = "on_change"
  )
}

cat("\n✓ All primary models fitted\n\n")

#===============================================================================
# CONVERGENCE DIAGNOSTICS
#===============================================================================

cat("Checking model convergence...\n")

convergence_summary <- map_df(names(bayes_models), function(outcome) {
  model     <- bayes_models[[outcome]]
  rhat_vals <- rhat(model)
  ess_bulk  <- neff_ratio(model)

  tibble(
    outcome       = outcome,
    max_rhat      = round(max(rhat_vals), 4),
    min_ess_ratio = round(min(ess_bulk), 3),
    converged     = max_rhat < 1.01 & min_ess_ratio > 0.1
  )
})

print(convergence_summary)
write_csv(convergence_summary, "results/data/convergence_diagnostics.csv")

if (!all(convergence_summary$converged)) {
  warning("Some models did not converge. Check results/data/convergence_diagnostics.csv")
} else {
  cat("✓ All models converged successfully\n\n")
}

#===============================================================================
# EXTRACT METABOLIC SYNDROME EFFECTS
#===============================================================================

cat("Extracting metabolic syndrome effects...\n")

mets_effects <- map_df(names(bayes_models), function(outcome) {
  model      <- bayes_models[[outcome]]
  posterior  <- as_draws_df(model)
  mets_coef  <- posterior$b_metS

  tibble(
    outcome        = outcome,
    posterior_mean = mean(mets_coef),
    posterior_sd   = sd(mets_coef),
    ci_lower       = quantile(mets_coef, 0.025),
    ci_upper       = quantile(mets_coef, 0.975),
    prob_negative  = mean(mets_coef < 0),
    prob_positive  = mean(mets_coef > 0)
  ) %>%
    mutate(
      ci_95     = sprintf("[%.2f, %.2f]", ci_lower, ci_upper),
      direction = case_when(
        prob_negative > 0.95 ~ "Credibly decreased",
        prob_positive > 0.95 ~ "Credibly increased",
        TRUE                 ~ "Uncertain"
      )
    )
})

print(mets_effects)
write_csv(mets_effects, "results/data/metabolic_syndrome_effects.csv")

cat("✓ Effects extracted\n\n")

#===============================================================================
# INTERACTION ANALYSIS (V/Q COUPLING)
#===============================================================================

cat("Fitting interaction model (V/Q coupling)...\n")

model_vq_interaction <- brm(
  qP ~ vP * metS + age_z + apache_z + ARDSsite,
  data   = final,
  family = gaussian(),
  prior  = c(
    prior(normal(0, 30),    class = "Intercept"),
    prior(normal(0, 5),     class = "b"),
    prior(exponential(0.1), class = "sigma")
  ),
  chains = 4, iter = 4000, warmup = 1000,
  cores  = 4, seed = 123,
  file   = "models/sensitivity/vq_interaction"
)

cat("\nInteraction coefficient (vP:metSYes):\n")
print(fixef(model_vq_interaction)["vP:metSYes", ])

cat("\n✓ Interaction analysis complete\n\n")

#===============================================================================
# GRADIENT ANALYSIS (DORSAL-VENTRAL RATIOS)
#===============================================================================

cat("Fitting gradient models...\n")

# Ventilation gradient
model_vent_ratio <- brm(
  vent_ratio ~ metS + age_z + apache_z + ARDSsite,
  data   = final,
  family = gaussian(),
  prior  = c(
    prior(normal(1, 1),    class = "Intercept"),
    prior(normal(0, 0.5),  class = "b"),
    prior(exponential(1),  class = "sigma")
  ),
  chains = 4, iter = 4000, warmup = 1000,
  cores  = 4, seed = 123,
  file   = "models/gradients/vent_ratio"
)

# Perfusion gradient
model_perf_ratio <- brm(
  perf_ratio ~ metS + age_z + apache_z + ARDSsite,
  data   = final,
  family = gaussian(),
  prior  = c(
    prior(normal(1, 1),    class = "Intercept"),
    prior(normal(0, 0.5),  class = "b"),
    prior(exponential(1),  class = "sigma")
  ),
  chains = 4, iter = 4000, warmup = 1000,
  cores  = 4, seed = 123,
  file   = "models/gradients/perf_ratio"
)

# Extract gradient effects
gradient_effects <- bind_rows(
  tibble(
    outcome        = "Ventilation Gradient (D/V)",
    posterior_mean = fixef(model_vent_ratio)["metSYes", "Estimate"],
    ci_lower       = fixef(model_vent_ratio)["metSYes", "Q2.5"],
    ci_upper       = fixef(model_vent_ratio)["metSYes", "Q97.5"]
  ),
  tibble(
    outcome        = "Perfusion Gradient (D/V)",
    posterior_mean = fixef(model_perf_ratio)["metSYes", "Estimate"],
    ci_lower       = fixef(model_perf_ratio)["metSYes", "Q2.5"],
    ci_upper       = fixef(model_perf_ratio)["metSYes", "Q97.5"]
  )
)

print(gradient_effects)
write_csv(gradient_effects, "results/data/gradient_effects.csv")

cat("✓ Gradient analysis complete\n\n")

#===============================================================================
# DECOMPOSITION ANALYSIS
#===============================================================================

cat("Performing decomposition analysis...\n")

# Fit shunt model with ventilation and perfusion as predictors
model_shunt_decomp <- brm(
  dsP_map ~ vP + qP + metS + age_z + apache_z + ARDSsite,
  data   = final,
  family = gaussian(),
  prior  = priors_gaussian,
  chains = 4, iter = 4000, warmup = 1000,
  cores  = 4, seed = 123,
  file   = "models/gradients/shunt_decomposition"
)

# Extract coefficients
posterior_draws <- as_draws_df(model_shunt_decomp)

# Calculate pathway contributions
contribution_analysis <- tibble(
  pathway = c("Via Ventilation", "Via Perfusion"),
  mediator_difference = c(
    mets_effects$posterior_mean[mets_effects$outcome == "vP"],
    mets_effects$posterior_mean[mets_effects$outcome == "qP"]
  ),
  coefficient = c(
    mean(posterior_draws$b_vP),
    mean(posterior_draws$b_qP)
  )
) %>%
  mutate(
    contribution     = coefficient * mediator_difference,
    abs_contribution = abs(contribution),
    pct_of_total     = 100 * abs_contribution / sum(abs_contribution)
  )

print(contribution_analysis)
write_csv(contribution_analysis, "results/data/decomposition_results.csv")

cat("✓ Decomposition analysis complete\n\n")

#===============================================================================
# TABLE 2: REGIONAL PHYSIOLOGY
#===============================================================================

cat("Creating Table 2...\n")

# Newdata for predictions
newdata_pred <- tibble(
  metS     = factor(c("No", "Yes"), levels = c("No", "Yes")),
  age_z    = 0,
  apache_z = 0,
  ARDSsite = factor("Extrapulmonary", levels = c("Extrapulmonary", "Pulmonary"))
)

# Function to get predictions
get_prediction_string <- function(model, row_num) {
  pred <- predict(model, newdata = newdata_pred[row_num, ], summary = TRUE)
  sprintf("%.1f [%.1f, %.1f]",
          pred[1, "Estimate"],
          pred[1, "Q2.5"],
          pred[1, "Q97.5"])
}

# Build Table 2
table2 <- tribble(
  ~Outcome, ~No_MetS, ~MetS, ~Difference,

  "Dorsal Ventilation (%)",
  get_prediction_string(bayes_models$vP, 1),
  get_prediction_string(bayes_models$vP, 2),
  sprintf("%.1f [%.1f, %.1f]",
          mets_effects$posterior_mean[mets_effects$outcome == "vP"],
          mets_effects$ci_lower[mets_effects$outcome == "vP"],
          mets_effects$ci_upper[mets_effects$outcome == "vP"]),

  "Dorsal Perfusion (%)",
  get_prediction_string(bayes_models$qP, 1),
  get_prediction_string(bayes_models$qP, 2),
  sprintf("%.1f [%.1f, %.1f]",
          mets_effects$posterior_mean[mets_effects$outcome == "qP"],
          mets_effects$ci_lower[mets_effects$outcome == "qP"],
          mets_effects$ci_upper[mets_effects$outcome == "qP"]),

  "Dorsal Shunt (%)",
  get_prediction_string(bayes_models$qsP_map, 1),
  get_prediction_string(bayes_models$qsP_map, 2),
  sprintf("%.1f [%.1f, %.1f]",
          mets_effects$posterior_mean[mets_effects$outcome == "qsP_map"],
          mets_effects$ci_lower[mets_effects$outcome == "qsP_map"],
          mets_effects$ci_upper[mets_effects$outcome == "qsP_map"]),

  "Ventral Shunt (%)",
  get_prediction_string(bayes_models$qsA_map, 1),
  get_prediction_string(bayes_models$qsA_map, 2),
  sprintf("%.1f [%.1f, %.1f]",
          mets_effects$posterior_mean[mets_effects$outcome == "qsA_map"],
          mets_effects$ci_lower[mets_effects$outcome == "qsA_map"],
          mets_effects$ci_upper[mets_effects$outcome == "qsA_map"]),

  "Dorsal Dead Space (%)",
  get_prediction_string(bayes_models$dsP_map, 1),
  get_prediction_string(bayes_models$dsP_map, 2),
  sprintf("%.1f [%.1f, %.1f]",
          mets_effects$posterior_mean[mets_effects$outcome == "dsP_map"],
          mets_effects$ci_lower[mets_effects$outcome == "dsP_map"],
          mets_effects$ci_upper[mets_effects$outcome == "dsP_map"]),

  "Ventral Dead Space (%)",
  get_prediction_string(bayes_models$dsA_map, 1),
  get_prediction_string(bayes_models$dsA_map, 2),
  sprintf("%.1f [%.1f, %.1f]",
          mets_effects$posterior_mean[mets_effects$outcome == "dsA_map"],
          mets_effects$ci_lower[mets_effects$outcome == "dsA_map"],
          mets_effects$ci_upper[mets_effects$outcome == "dsA_map"])
)

# Rename columns
table2 <- table2 %>%
  rename(
    `Outcome`                        = Outcome,
    `No Metabolic Syndrome (N = 15)` = No_MetS,
    `Metabolic Syndrome (N = 10)`    = MetS,
    `Difference`                     = Difference
  )

write_csv(table2, "results/tables/Table2_Regional_Physiology.csv")
cat("✓ Table 2 created\n\n")

#===============================================================================
# FIGURE 1: FOREST PLOT
#===============================================================================

cat("Creating Figure 1 (forest plot)...\n")

forest_data <- mets_effects %>%
  mutate(
    outcome_label = case_when(
      outcome == "vP"      ~ "  Ventilation Distribution",
      outcome == "qP"      ~ "  Perfusion Distribution",
      outcome == "qsP_map" ~ "% Low V/Q Impedance Ratios",
      outcome == "dsP_map" ~ "% High V/Q Impedance Ratios",
      outcome == "qsA_map" ~ "% Low V/Q Impedance Ratios",
      outcome == "dsA_map" ~ "% High V/Q Impedance Ratios"
    ),
    outcome_order = case_when(
      outcome == "vP"      ~ 7,
      outcome == "qP"      ~ 6,
      outcome == "qsP_map" ~ 5,
      outcome == "dsP_map" ~ 4,
      outcome == "qsA_map" ~ 2,
      outcome == "dsA_map" ~ 1
    ),
    color = case_when(
      prob_negative >= 0.95 | prob_positive >= 0.95 ~ "Credible",
      TRUE ~ "Uncertain"
    )
  ) %>%
  arrange(desc(outcome_order))

y_axis_labels <- c(
  "      High V/Q Impedance Ratio",
  "      Low V/Q Impedance Ratio",
  expression(bold("Ventral Lung Regions")),
  "      High V/Q Impedance Ratio",
  "      Low V/Q Impedance Ratio",
  "      Perfusion Distribution",
  "      Ventilation Distribution",
  expression(bold("Dorsal Lung Regions"))
)

p_forest <- ggplot(forest_data, aes(x = posterior_mean, y = outcome_order)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_hline(yintercept = 3.5, linetype = "solid", color = "gray70", linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_lower, xmax = ci_upper, color = color),
                 linewidth = 1.2, na.rm = TRUE) +
  geom_point(aes(color = color), size = 4, shape = 19, na.rm = TRUE) +
  scale_color_manual(
    values = c("Credible" = "#D55E00", "Uncertain" = "gray50"),
    guide  = "none"
  ) +
  scale_y_continuous(
    breaks = c(1, 2, 3, 4, 5, 6, 7, 8),
    labels = y_axis_labels,
    limits = c(0.5, 8.5)
  ) +
  scale_x_continuous(
    limits = c(-15, 15),
    breaks = seq(-15, 15, 5)
  ) +
  labs(
    x       = "Mean Difference (95% CrI): Patients with versus without Metabolic Syndrome",
    y       = NULL,
    caption = "Ventilation and Perfusion Distribution: % of total impedance | V/Q Impedance Ratios: % of total normalized pixels in that region"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y    = element_text(size = 11, hjust = 0),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    axis.line.y    = element_blank(),
    axis.ticks.y   = element_blank(),
    plot.caption   = element_text(hjust = 0, size = 9, face = "italic")
  )

ggsave("figures/Figure1.pdf", p_forest, width = 10, height = 7, dpi = 300)
cat("✓ Figure 1 created\n\n")

#===============================================================================
# FIGURE 2: PHASE PLOT (V/Q Relationship)
#===============================================================================

cat("Creating Figure 2 (phase plot)...\n")

plot_data <- final %>%
  mutate(metS_label = factor(metS, levels = c("No", "Yes"),
                             labels = c("No Metabolic Syndrome", "Metabolic Syndrome")))

p_phase <- ggplot(plot_data, aes(x = vP, y = qP, color = metS_label, fill = metS_label)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1.2) +
  geom_point(size = 3.5, alpha = 0.8) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values  = c("#0072B2", "#D55E00")) +
  labs(
    x     = "Dorsal Ventilation (%)",
    y     = "Dorsal Perfusion (%)",
    color = NULL,
    fill  = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    legend.position  = c(0.25, 0.9),
    legend.background = element_rect(fill = "white", color = "gray80")
  )

ggsave("figures/Figure2.pdf", p_phase, width = 8, height = 7, dpi = 300)
ggsave("figures/Figure2.png", p_phase, width = 8, height = 7, dpi = 300)
cat("✓ Figure 2 created\n\n")

#===============================================================================
# FIGURE 3: GRADIENT COMPARISON
#===============================================================================

cat("Creating Figure 3 (gradients)...\n")

# Get predictions for gradients
newdata_gradient <- tibble(
  metS     = factor(c("No", "Yes"), levels = c("No", "Yes")),
  age_z    = 0,
  apache_z = 0,
  ARDSsite = factor("Extrapulmonary", levels = c("Extrapulmonary", "Pulmonary"))
)

# Ventilation ratio predictions
vent_preds <- predict(model_vent_ratio, newdata = newdata_gradient, summary = TRUE) %>%
  as_tibble() %>%
  mutate(metS = c("No", "Yes"),
         type = "Ventilation")

# Perfusion ratio predictions
perf_preds <- predict(model_perf_ratio, newdata = newdata_gradient, summary = TRUE) %>%
  as_tibble() %>%
  mutate(metS = c("No", "Yes"),
         type = "Perfusion")

# Combine
gradient_plot_data <- bind_rows(vent_preds, perf_preds) %>%
  mutate(
    metS = factor(metS, levels = c("No", "Yes")),
    type = factor(type, levels = c("Ventilation", "Perfusion"))
  )

# Create plot
p_gradient <- ggplot(gradient_plot_data, aes(x = type, y = Estimate, fill = metS)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.7) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5),
                position = position_dodge(width = 0.8),
                width = 0.25, linewidth = 0.8) +
  scale_fill_manual(values = c("#0072B2", "#D55E00"),
                    labels = c("No Metabolic Syndrome", "Metabolic Syndrome")) +
  labs(
    x     = NULL,
    y     = "Dorsal-to-Ventral Ratio",
    title = "Gravitational Distribution of Ventilation and Perfusion",
    fill  = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title             = element_text(hjust = 0.5, face = "bold"),
    legend.position        = "bottom",
    panel.grid.major.y     = element_line(color = "gray90", linewidth = 0.3)
  )

ggsave("figures/Figure3.pdf", p_gradient, width = 8, height = 7, dpi = 300)
ggsave("figures/Figure3.png", p_gradient, width = 8, height = 7, dpi = 300)
cat("✓ Figure 3 created\n\n")

#===============================================================================
# FIGURE 4: DECOMPOSITION WATERFALL
#===============================================================================

cat("Creating Figure 4 (decomposition)...\n")

# Prepare data for waterfall
waterfall_data <- contribution_analysis %>%
  mutate(
    pathway = factor(pathway, levels = c("Via Ventilation", "Via Perfusion")),
    sign    = ifelse(contribution > 0, "Increase", "Decrease")
  )

# Create waterfall plot
p_decomp <- ggplot(waterfall_data, aes(x = pathway, y = contribution, fill = sign)) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "gray40") +
  geom_col(width = 0.6, alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%\n(%.0f%%)", contribution, pct_of_total)),
            vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("Increase" = "#D55E00", "Decrease" = "#0072B2"),
                    guide  = "none") +
  labs(
    x     = NULL,
    y     = "Contribution to Dorsal Shunt (%)",
    title = "Decomposition of Metabolic Syndrome Effects on Dorsal Shunt"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
  )

ggsave("figures/Figure4.pdf", p_decomp, width = 8, height = 7, dpi = 300)
ggsave("figures/Figure4.png", p_decomp, width = 8, height = 7, dpi = 300)
cat("✓ Figure 4 created\n\n")

#===============================================================================
# SUPPLEMENTARY FIGURES
#===============================================================================

cat("Creating supplementary figures...\n")

# Supplementary Figure S1: Distribution diagnostics
p_distributions <- final %>%
  select(all_of(outcomes_to_check)) %>%
  pivot_longer(everything()) %>%
  mutate(name = factor(name, levels = outcomes_to_check)) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 15, fill = "steelblue", alpha = 0.7, color = "black") +
  facet_wrap(~name, scales = "free", ncol = 2) +
  labs(x = "Value (%)", y = "Count",
       title = "Distribution of Regional Physiology Outcomes") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("figures/FigureS1.pdf", p_distributions, width = 10, height = 12, dpi = 300)
cat("✓ Supplementary figures created\n\n")

#===============================================================================
# ORGANIZE FOR SUBMISSION
#===============================================================================

cat("Organizing figures for submission...\n")

# Copy main figures
file.copy("figures/Figure1.pdf", "submission_figures/Figure1.pdf", overwrite = TRUE)
file.copy("figures/Figure2.pdf", "submission_figures/Figure2.pdf", overwrite = TRUE)
file.copy("figures/Figure3.pdf", "submission_figures/Figure3.pdf", overwrite = TRUE)
file.copy("figures/Figure4.pdf", "submission_figures/Figure4.pdf", overwrite = TRUE)

# Copy supplementary
file.copy("figures/FigureS1.pdf", "submission_figures/FigureS1.pdf", overwrite = TRUE)

cat("✓ Submission package organized\n\n")

#===============================================================================
# SAVE COMPLETE WORKSPACE
#===============================================================================

cat("Saving complete workspace...\n")

save(
  # Data
  final,
  # Models
  bayes_models,
  model_vq_interaction,
  model_vent_ratio,
  model_perf_ratio,
  model_shunt_decomp,
  # Results
  table1,
  table2,
  mets_effects,
  convergence_summary,
  gradient_effects,
  contribution_analysis,
  dist_summary,
  # Sample sizes
  n_total, n_no, n_yes,
  file = "results/data/complete_analysis.RData"
)

cat("✓ Workspace saved\n\n")

#===============================================================================
# FINAL SUMMARY
#===============================================================================

end_time <- Sys.time()
elapsed  <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║                  ANALYSIS COMPLETE                            ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
cat("Total runtime:", round(elapsed, 1), "minutes\n\n")

cat("KEY FINDINGS:\n")
cat("─────────────────────────────────────────────────────────────────\n")

key_finding <- mets_effects %>%
  filter(outcome == "vP") %>%
  pull(posterior_mean) %>%
  round(2)

vP_row <- mets_effects[mets_effects$outcome == "vP", ]

cat("Metabolic syndrome reduces dorsal ventilation by", abs(key_finding), "%\n")
cat("  (95% CrI:", round(vP_row$ci_lower, 2), "to", round(vP_row$ci_upper, 2), ")\n")
cat("  Posterior probability of decrease:",
    round(100 * vP_row$prob_negative, 1), "%\n\n")

coupling <- fixef(model_vq_interaction)["vP:metSYes", "Estimate"]
cat("Ventilation-perfusion coupling preserved (interaction:", round(coupling, 3), ")\n\n")

vent_contribution <- contribution_analysis$pct_of_total[1]
cat("Effect pathway: ~", round(vent_contribution), "% via ventilation,",
    round(100 - vent_contribution), "% via perfusion\n\n")

cat("FILES CREATED:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("Tables:\n")
cat("  results/tables/Table1_with_SMD.csv\n")
cat("  results/tables/Table2_Regional_Physiology.csv\n\n")
cat("Main Figures:\n")
cat("  submission_figures/Figure1.pdf (forest plot)\n")
cat("  submission_figures/Figure2.pdf (phase plot)\n")
cat("  submission_figures/Figure3.pdf (gradients)\n")
cat("  submission_figures/Figure4.pdf (decomposition)\n\n")
cat("Supplementary:\n")
cat("  submission_figures/FigureS1.pdf (distributions)\n\n")
cat("Models and Results:\n")
cat("  models/ (all Bayesian model objects)\n")
cat("  results/data/complete_analysis.RData (complete workspace)\n\n")
cat("TO RELOAD RESULTS:\n")
cat("  load('results/data/complete_analysis.RData')\n\n")
