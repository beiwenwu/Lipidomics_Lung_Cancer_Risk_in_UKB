library(dplyr)
library(Matrix)
library(stats)
library(survival)
library(glmnet)
library(purrr)
library(tibble)
library(writexl)

################ Run LASSO-Cox with Relative Metabolites ################
###### Define one common LASSO sample - Restrict every imputation to those same participants ###### 
# The same procedure as used in PCA-Cox in getting data_ratio_complete and common_ids_ratio 

# Create fixed 10-fold cross-validation assignments, so so that the same participants are assigned to the same folds across all imputed datasets
set.seed(123)
nfolds <- 10

# Use the first metabolite-complete imputed dataset as the reference dataset
reference_data <- data_ratio_complete[[1]]

# Create a lookup table linking each participant ID to a fixed CV fold
foldid_lookup_ratio <- data.frame(
  eid = common_ids_ratio,
  foldid = sample(
    rep(seq_len(nfolds), length.out = length(common_ids_ratio))
  )
)

# Check the distribution of lung cancer cases and non-cases across the 10 cross-validation folds
table(
  foldid_lookup_ratio$foldid,
  reference_data[["lungc_status_coxph"]]
)

# Merge fold assignments to the reference dataset by participant ID
data_ratio_complete <- lapply(
  data_ratio_complete,
  function(dat) {
    left_join(
      dat,
      foldid_lookup_ratio,
      by = "eid"
    )
  }
)

# Pack_year is log-transformed to improve model stability
data_ratio_complete <- lapply(
  data_ratio_complete,
  function(dat) {
    dat$log_pack_year <- log1p(dat$pack_year_imputed)
    dat
  }
)

###### Construct the predictor matrix ###### 
extended_covariates_lasso <- c(
  "age",
  "sex",
  "ethnicity_UKB",
  "education_ISCED",
  "tdi",
  "bmi",
  "bmi2",
  "waist_to_hip",
  "smk_status",
  "log_pack_year",
  "family_lc",
  "lung_disease_history",
  "lipid_lowering_meds",
  "summed_MET_minutes_weekly_imputed",
  "fasting_time",
  "centre_group"
)

# Create the LASSO-Cox design matrix
predictor_matrix_design <- function(dat, metabolites, extended_covariates) {
  # Penalized metabolite matrix
  x_met <- Matrix(
    as.matrix(
      dat[, metabolites, drop = FALSE]
    ),
    sparse = TRUE
  )
  
  # Unpenalized covariate matrix
  # Categorical factors become dummy variables
  x_cov <- sparse.model.matrix(
    reformulate(extended_covariates_lasso),
    data = dat
  )
  
  # Remove the intercept automatically created by sparse.model.matrix()
  x_cov <- x_cov[, colnames(x_cov) != "(Intercept)", drop = FALSE]
  
  # Combine the unpenalized covariates and penalized metabolites into a single predictor matrix
  x <- cbind(x_cov,x_met)
  
  # Define the penalty factor:
  # 0 = covariate is not penalized and is retained in the model
  # 1 = metabolite is penalized and eligible for LASSO selection
  penalty_factor <- c(
    rep(0, ncol(x_cov)),
    rep(1, ncol(x_met))
  )
  
  # Label each penalty factor with its corresponding predictor name
  names(penalty_factor) <- colnames(x)
  
  # Return the design matrix and corresponding penalty factors
  list(
    x = x,
    penalty_factor = penalty_factor
  )
}

# Construct the LASSO-Cox design matrix separately for each datasets
predictor_matrix_design_list <- lapply(
  X = data_ratio_complete,
  FUN = predictor_matrix_design,
  metabolites = metabolite_ratio,
  extended_covariates = extended_covariates_lasso
)

###### Run the LASSO-Cox models ###### 
# Fit a LASSO-Cox model separately within each imputed dataset
for (i in seq_along(data_ratio_complete)) {
  
  # Print progress to the console
  message("Fitting imputation ", i)
  
  # Extract the current dataset
  dat <- data_ratio_complete[[i]]
  
  # Extract the corresponding LASSO design matrix
  x <- as.matrix(predictor_matrix_design_list[[i]]$x)
  
  # Create the survival object
  y <- Surv(
    time = dat$lungc_time,
    event = dat$lungc_status_coxph
  )
  
  # Fit a 10-fold cross-validated LASSO-Cox model
  lasso_fits_ratio[[i]] <- cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 1,
    foldid = dat$foldid, # Use the fixed CV fold assignments
    penalty.factor = predictor_matrix_design_list[[i]]$penalty_factor,
    standardize = FALSE,
    type.measure = "deviance" # Select lambda using deviance method
  )
}

# Extract the full LASSO-Cox results from each dataset using lambda.min
lasso_results_ratio_coef_min <- imap_dfr(
  lasso_fits_ratio,
  function(fit, i) {
    # Extract the coefficient vector at lambda.min
    beta <- as.matrix(
      coef(fit, s = "lambda.min")
    )
    # Store the dataset number, lambda.min, variable name, and corresponding coefficient
    tibble(
      imputation = i,
      lambda_min = fit$lambda.min,
      variable = rownames(beta),
      coefficient = as.numeric(beta[, 1])
    ) |>
      # Retain only predictors with non-zero coefficients
      filter(coefficient != 0)
  }
)

# Extract the full LASSO-Cox results from each dataset using lambda.1se
lasso_results_ratio_coef_1se <- imap_dfr(
  lasso_fits_ratio,
  function(fit, i) {
    # Extract the coefficient vector at lambda.1se
    beta <- as.matrix(
      coef(fit, s = "lambda.1se")
    )
    # Store the dataset number, lambda.min, variable name, and corresponding coefficient
    tibble(
      imputation = i,
      lambda_1se = fit$lambda.1se,
      variable = rownames(beta),
      coefficient = as.numeric(beta[, 1])
    ) |>
      # Retain only predictors with non-zero coefficients
      filter(coefficient != 0)
  }
)

###### Extract metabolites selected across 10 datasets ######
extract_selected <- function(fit, imputation, s, metabolites) {
  # Extract the coefficient vector at the chosen lambda
  beta <- as.matrix(
    coef(fit, s = s)
  )
  # Keep only penalized metabolites with non-zero coefficients
  data.frame(
    imputation = imputation,
    metabolite = rownames(beta),
    coefficient = as.numeric(beta[, 1])
  ) |>
    filter(
      metabolite %in% metabolites,
      coefficient != 0
    )
}

# Extract metabolites selected at lambda.1se from each dataset
selected_1se_ratio <- imap_dfr(
  lasso_fits_ratio,
  ~ extract_selected(
    fit = .x,
    imputation = .y,
    s = "lambda.1se",
    metabolites = metabolite_ratio
  )
)

# Extract metabolites selected at lambda.min from each dataset
selected_min_ratio <- imap_dfr(
  lasso_fits_ratio,
  ~ extract_selected(
    fit = .x,
    imputation = .y,
    s = "lambda.min",
    metabolites = metabolite_ratio
  )
)

# Summarize the frequency of each metabolite selected across the 10 datasets at lambda.min
selection_frequency_min_ratio <- selected_min |>
  group_by(metabolite) |>
  summarise(
    n_selected = n_distinct(imputation),
    selection_frequency = n_selected / 10,
    positive_n = sum(coefficient > 0),
    negative_n = sum(coefficient < 0),
    .groups = "drop"
  ) |>
  arrange(
    desc(n_selected)
  )

###### Save the LASSO-Cox results ###### 
write_xlsx(lasso_results_ratio_coef_1se,"lasso_results_ratio_coef_1se.xlsx")
write_xlsx(lasso_results_ratio_coef_min,"lasso_results_ratio_coef_min.xlsx")
saveRDS(lasso_fits_ratio, "lasso_fits_ratio.rds")
write_xlsx(selection_frequency_min_ratio,"lasso_ratio_summary_min.xlsx")

################# Run LASSO-Cox with All Metabolites #################
#### Define one common LASSO sample - Restrict every imputation to those same participants ####
# Remove participants with missing values in any metabolite
data_overall_metabolites <- lapply(
  imputed_data,
  function(dat) {
    dat %>%
      drop_na(all_of(metabolites))
  }
)

# Retrieve the number of participants in each imputed dataset
sapply(data_overall_metabolites, nrow)

# Extract the participant IDs from each metabolite-complete imputed dataset
common_ids_overall_each <- lapply(
  data_overall_metabolites,
  function(dat) {
    dat[["eid"]]
  }
)

# Identify the common IDs across all ten metabolite-complete imputed datasets N=348,207
common_ids_overall <- Reduce(
  intersect,
  common_ids_overall_each
)

# Keep only the common participants across all ten metabolite-complete imputed datasets
data_overall_metabolites_complete <- lapply(
  data_overall_metabolites,
  function(dat) {
    idx <- match(common_ids_overall, dat[["eid"]])
    dat[idx, , drop = FALSE]
  }
)

# Create fixed 10-fold cross-validation assignments, so so that the same participants are assigned to the same folds across all imputed datasets
set.seed(123)
nfolds <- 10

# Use the first metabolite-complete imputed dataset as the reference dataset
reference_data <- data_overall_metabolites_complete[[1]]

# Create a lookup table linking each participant ID to a fixed CV fold
foldid_lookup_overall <- data.frame(
  eid = common_ids_overall,
  foldid = sample(
    rep(seq_len(nfolds), length.out = length(common_ids_overall))
  )
)

# Check the distribution of lung cancer cases and non-cases across the 10 cross-validation folds
table(
  foldid_lookup_overall$foldid,
  reference_data[["lungc_status_coxph"]]
)

# Merge fold assignments to the reference dataset by participant ID
data_overall_metabolites_complete <- lapply(
  data_overall_metabolites_complete,
  function(dat) {
    left_join(
      dat,
      foldid_lookup_overall,
      by = "eid"
    )
  }
)

# Pack_year is log-transformed to improve model stability
data_overall_metabolites_complete <- lapply(
  data_overall_metabolites_complete,
  function(dat) {
    dat$log_pack_year <- log1p(dat$pack_year_imputed)
    dat
  }
)

###### Construct the predictor matrix ###### 
# Construct the LASSO-Cox design matrix separately for each datasets
predictor_matrix_design_list <- lapply(
  X = data_overall_metabolites_complete,
  FUN = predictor_matrix_design,
  metabolites = metabolites,
  extended_covariates = extended_covariates_lasso
)

###### Run the LASSO-Cox models ###### 
# Fit a LASSO-Cox model separately within each imputed dataset
for (i in seq_along(data_overall_metabolites_complete)) {
  
  # Print progress to the console
  message("Fitting imputation ", i)
  
  # Extract the current dataset
  dat <- data_overall_metabolites_complete[[i]]
  
  # Extract the corresponding LASSO design matrix
  x <- as.matrix(predictor_matrix_design_list[[i]]$x)
  
  # Create the survival object
  y <- Surv(
    time = dat$lungc_time,
    event = dat$lungc_status_coxph
  )
  
  # Fit a 10-fold cross-validated LASSO-Cox model
  lasso_fits_overall[[i]] <- cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 1,
    foldid = dat$foldid, # Use the fixed CV fold assignments
    penalty.factor = predictor_matrix_design_list[[i]]$penalty_factor,
    standardize = FALSE,
    type.measure = "deviance" # Select lambda using deviance method
  )
}

# Extract the LASSO-Cox results from each dataset using lambda.min
lasso_results_overall_coef_min <- imap_dfr(
  lasso_fits_overall,
  function(fit, i) {
    # Extract the coefficient vector at lambda.min
    beta <- as.matrix(
      coef(fit, s = "lambda.min")
    )
    # Store the dataset number, lambda.min, variable name, and corresponding coefficient
    tibble(
      imputation = i,
      lambda_min = fit$lambda.min,
      variable = rownames(beta),
      coefficient = as.numeric(beta[, 1])
    ) |>
      # Retain only predictors with non-zero coefficients
      filter(coefficient != 0)
  }
)

# Extract the LASSO-Cox results from each dataset using lambda.1se
lasso_results_overall_coef_1se <- imap_dfr(
  lasso_fits_overall,
  function(fit, i) {
    # Extract the coefficient vector at lambda.1se
    beta <- as.matrix(
      coef(fit, s = "lambda.1se")
    )
    # Store the dataset number, lambda.min, variable name, and corresponding coefficient
    tibble(
      imputation = i,
      lambda_1se = fit$lambda.1se,
      variable = rownames(beta),
      coefficient = as.numeric(beta[, 1])
    ) |>
      # Retain only predictors with non-zero coefficients
      filter(coefficient != 0)
  }
)

###### Extract metabolites selected across 10 datasets ######
# Extract metabolites selected at lambda.1se from each dataset
selected_1se_overall <- imap_dfr(
  lasso_fits_overall,
  ~ extract_selected(
    fit = .x,
    imputation = .y,
    s = "lambda.1se",
    metabolites = metabolites
  )
)

# Extract metabolites selected at lambda.min from each dataset
selected_min_overall <- imap_dfr(
  lasso_fits_overall,
  ~ extract_selected(
    fit = .x,
    imputation = .y,
    s = "lambda.min",
    metabolites = metabolites
  )
)

# Summarize the frequency of each metabolite selected across the 10 datasets at lambda.min
selection_frequency_min_overall <- selected_min_overall |>
  group_by(metabolite) |>
  summarise(
    n_selected = n_distinct(imputation),
    selection_frequency = n_selected / 10,
    positive_n = sum(coefficient > 0),
    negative_n = sum(coefficient < 0),
    .groups = "drop"
  ) |>
  arrange(
    desc(n_selected)
  )

###### Save the LASSO-Cox results ###### 
write_xlsx(lasso_results_overall_coef_1se,"lasso_results_overall_coef_1se.xlsx")
write_xlsx(lasso_results_overall_coef_min,"lasso_results_overall_coef_min.xlsx")
saveRDS(lasso_fits_overall, "lasso_fits_overall.rds")
write_xlsx(selection_frequency_min_overall,"lasso_overall_summary_min.xlsx")

################ Run LASSO-selected Metabolites in a Multivariable CoxPH Model ################
lasso_ratio_metabolites <- c(
  "IDL_PL_%",
  "LA_%",
  "L_HDL_FC_%",
  "L_LDL_CE_%",
  "S_LDL_CE_%",
  "XS_VLDL_FC_%",
  "XXL_VLDL_PL_%",
  "XXL_VLDL_TG_%"
)

lasso_selected_ratio_fits <- lapply(data_ratio_complete, function(dat) {
  
  # Add backticks around metabolite names containing %
  metab_terms <- paste0("`", lasso_ratio_metabolites, "`")
  
  rhs <- paste(
    c(metab_terms, extended_covariates),
    collapse = " + "
  )
  
  form <- as.formula(
    paste0(
      "Surv(lungc_time, lungc_status_coxph) ~ ",
      rhs
    )
  )
  
  coxph(form, data = dat)
})

# Convert the list of Cox models to a mira object
lasso_selected_ratio_mira <- as.mira(lasso_selected_ratio_fits)

# Pool coefficients using Rubin's rules
lasso_selected_ratio_pool <- pool(lasso_selected_ratio_mira)

# Obtain pooled HRs and 95% CIs
lasso_selected_ratio_results <- summary(
  lasso_selected_ratio_pool,
  conf.int = TRUE,
  exponentiate = TRUE
)

# Retain the rows for LASSO-selected metabolites
lasso_selected_ratio_table <- lasso_selected_ratio_results %>%
  slice(1:8) %>%
  select(
    Metabolite = term,
    HR = estimate,
    LCL = `2.5 %`,
    UCL = `97.5 %`,
    p = p.value
  ) %>%
  mutate(
    Metabolite = gsub("`", "", Metabolite, fixed = TRUE)
  ) %>%
  mutate(
    p_fmt = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )

write_xlsx(lasso_selected_ratio_table,"lasso_selected_ratio_table.xlsx")