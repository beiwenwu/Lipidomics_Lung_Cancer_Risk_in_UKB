library(survival)
library(survminer)
library(car)
library(cmprsk)
library(writexl)

# Extract participant ID and metabolites from the CoxPH analysis dataset
metabolite_data <- coxph_data %>%
  select(eid, all_of(metabolites))

# Define covariates included in the fully adjusted CoxPH models
extended_covariates <- c(
  "age",
  "sex",
  "ethnicity_UKB",
  "education_ISCED",
  "tdi",
  "bmi",
  "bmi2",
  "waist_to_hip",
  "smk_status",
  "pack_year_imputed",
  "family_lc",
  "lung_disease_history",
  "lipid_lowering_meds",
  "summed_MET_minutes_weekly_imputed",
  "fasting_time",
  "centre_group"
)

######################### Create a CoxPH Function for one metabolite #########################
run_coxph_analysis_mi <- function(metabolite, imp_mids, covariates_use, metabolite_data, stratified_var = NULL, subgroup = NULL) {
  
  # Create an empty list to store the Cox model fitted in each imputed dataset
  fits <- vector("list", imp_mids$m)
  
  # Construct the right-hand side of the CoxPH formula 
  rhs <- paste(
    c(
      paste0("`", metabolite, "`"),
      paste0("`", covariates_use, "`")
    ),
    collapse = " + "
  )
  
  # Construct the CoxPH model formula
  form <- as.formula(
    paste0("Surv(lungc_time, lungc_status_coxph) ~ ", rhs)
  )
  
  # Fit the metabolite-specific Cox model separately in each imputed datasets
  for (i in seq_len(imp_mids$m)) {
    
    # Extract the ith completed dataset and merge with the metabolite 
    dat_i <- complete(imp_mids, i) %>%
      left_join(
        metabolite_data %>% select(eid, all_of(metabolite)),
        by = "eid"
      )
    
    # Extract the completed dataset to participants in the specified subgroup
    if (!is.null(stratified_var) & !is.null(subgroup)) {
      dat_i <- dat_i %>%
        filter(.data[[stratified_var]] == subgroup)
    }
    
    # Fit the CoxPH model using Efron method for tired events
    fits[[i]] <- withCallingHandlers(
      coxph(
        form,
        data = dat_i,
        ties = "efron",
        x = TRUE
      ),
      warning = function(w) {
        message(
          "Warning | metabolite: ", metabolite,
          " | imputation: ", i,
          " | message: ", conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    )
    
    # Remove the temporary completed dataset to reduce memory usage
    rm(dat_i)
    gc()
  }
  
  # Convert the list of CoxPH fits to a mira object for pooling across imputations
  fit_mira <- as.mira(fits)
  
  # Pool model coefficients and standard errors across imputations using Rubin's rules and report hazard ratios with 95% confidence intervals
  pooled_sum <- summary(
    pool(fit_mira),
    conf.int = TRUE,
    exponentiate = TRUE
  )
  
  # Extract the pooled result corresponding to the metabolite of interest
  met_row <- pooled_sum[1, ]
  
  # Extract the median analytic sample size and number of lung cancer events corresponding to the metabolite of interest across the imputed datasets
  n_total <- median(sapply(fits, function(x) x$n))
  n_events <- median(sapply(fits, function(x) x$nevent))

  # Test the proportional hazards assumption within each imputed dataset using Schoenfeld residuals
  ph_tests <- lapply(fits, function(fit) {
    tryCatch(
      cox.zph(fit),
      error = function(e) NULL
    )
  })
  
  # Extract the proportional hazards test p-value for the metabolite from each imputed dataset
  ph_metabolite_p <- sapply(ph_tests, function(x) {
    if (is.null(x)) NA_real_ else x$table[1, "p"]
  })
  
  # Extract the global proportional hazards test p-value from each imputed dataset
  ph_global_p <- sapply(ph_tests, function(x) {
    if (is.null(x)) NA_real_ else x$table["GLOBAL", "p"]
  })
  
  # Return the summarized results
  data.frame(
    Metabolite = metabolite,
    Subgroup = ifelse(is.null(subgroup), "Overall", subgroup),
    Hazard_Ratio = met_row$estimate,
    CI_Lower = met_row$`2.5 %`,
    CI_Upper = met_row$`97.5 %`,
    P_Value = met_row$p.value,
    N = n_total,
    N_events = n_events,
    PH_P_Global = median(ph_global_p, na.rm = TRUE),
    PH_P_Metabolite = median(ph_metabolite_p, na.rm = TRUE)
  )
}

######################### Create a CoxPH function for all metabolites #########################
run_full_analysis_mi <- function(metabolites, imp_mids, covariates_use, metabolite_data, stratified_var = NULL, subgroup = NULL) {
  
  # Run the CoxPH analysis separately for each metabolite
  all_results <- lapply(metabolites, function(metabolite) {
    
    # Print the name of the metabolite currently being analyzed
    cat("Running:", metabolite, "/n")
    
    # Apply the single-metabolite CoxPH function
    run_coxph_analysis_mi(
      metabolite = metabolite,
      imp_mids = imp_mids,
      covariates_use = covariates_use,
      metabolite_data = metabolite_data,
      stratified_var = stratified_var,
      subgroup = subgroup 
    )
  })
  
  # Combine the metabolite-specific results into one data frame
  combined_results <- do.call(rbind, all_results)
  
  # Adjust p-values for multiple testing using the Benjamini-Hochberg false discovery rate method
  combined_results$P_Adjusted <- p.adjust(combined_results$P_Value, method = "BH")
  combined_results$PH_P_Global_Adjusted <- p.adjust(combined_results$PH_P_Global, method = "BH")
  combined_results$PH_P_Metabolite_Adjusted <- p.adjust(combined_results$PH_P_Metabolite, method = "BH")
  
  # Return the combined results
  list(Adjusted = combined_results)
}

######################### Run Primary CoxPH Analysis #########################
results_overall <- run_full_analysis_mi(
  metabolites = metabolites,
  imp_mids = imp_mids,
  covariates_use = extended_covariates,
  metabolite_data = metabolite_data
)

######################### Run Stratified CoxPH Analysis by Sex #########################
# Do not adjust sex when conduct stratified analysis by sex 
covariates_use <- setdiff(extended_covariates, "sex")

results_female <- run_full_analysis_mi(
  metabolites = metabolites,
  imp_mids = imp_mids,
  covariates_use = covariates_use,
  metabolite_data = metabolite_data,
  stratified_var = "sex",
  subgroup = "Female"
)

results_male <- run_full_analysis_mi(
  metabolites = metabolites,
  imp_mids = imp_mids,
  covariates_use = covariates_use,
  metabolite_data = metabolite_data,
  stratified_var = "sex",
  subgroup = "Male"
)

######################### Run Stratified CoxPH Analysis by smk status #########################
# Do not adjust smk status when conduct stratified analysis by smk status 
covariates_use <- setdiff(extended_covariates, "smk_status")

results_current <- run_full_analysis_mi(
  metabolites = metabolites,
  imp_mids = imp_mids,
  covariates_use = covariates_use,
  metabolite_data = metabolite_data,
  stratified_var = "smk_status",
  subgroup = "Current"
)

results_previous <- run_full_analysis_mi(
  metabolites = metabolites,
  imp_mids = imp_mids,
  covariates_use = covariates_use,
  metabolite_data = metabolite_data,
  stratified_var = "smk_status",
  subgroup = "Previous"
)

# Do not adjust smk status and pack year when conduct stratified analysis in never smokers 
covariates_use <- setdiff(extended_covariates, c("smk_status", "pack_year_imputed"))

results_never <- run_full_analysis_mi(
  metabolites = metabolites,
  imp_mids = imp_mids,
  covariates_use = covariates_use,
  metabolite_data = metabolite_data,
  stratified_var = "smk_status",
  subgroup = "Never"
)

######################### Run Stratified CoxPH Analysis by age group #########################
# Do not adjust age when conduct stratified analysis by age group 
covariates_use <- setdiff(extended_covariates, "age")

results_agegroup1 <- run_full_analysis_mi(
  metabolites = metabolites,
  imp_mids = imp_mids,
  covariates_use = covariates_use,
  metabolite_data = metabolite_data,
  stratified_var = "age_group",
  subgroup = "<=50"
)

results_agegroup2 <- run_full_analysis_mi(
  metabolites = metabolites,
  imp_mids = imp_mids,
  covariates_use = covariates_use,
  metabolite_data = metabolite_data,
  stratified_var = "age_group",
  subgroup = "(50-60]"
)

results_agegroup3 <- run_full_analysis_mi(
  metabolites = metabolites,
  imp_mids = imp_mids,
  covariates_use = covariates_use,
  metabolite_data = metabolite_data,
  stratified_var = "age_group",
  subgroup = ">60"
)

######################### Run Cause-specific CoxPH Analysis by histology #########################
# Include an event_indicator when create the CoxPH Function for a single metabolite
run_coxph_analysis_mi <- function(metabolite, imp_mids, covariates_use, metabolite_data, event_indicator) {
  
  fits <- vector("list", imp_mids$m)
  
  rhs <- paste(
    c(
      paste0("`", metabolite, "`"),
      paste0("`", covariates_use, "`")
    ),
    collapse = " + "
  )
  
  form <- as.formula(
    paste0("Surv(lungc_time, ", event_indicator, ") ~ ", rhs)
  )
  
  for (i in seq_len(imp_mids$m)) {
    
    dat_i <- complete(imp_mids, i) %>%
      left_join(
        metabolite_data %>% select(eid, all_of(metabolite)),
        by = "eid"
      )
    
    fits[[i]] <- withCallingHandlers(
      coxph(
        form,
        data = dat_i,
        ties = "efron",
        x = TRUE
      ),
      warning = function(w) {
        message(
          "Warning | metabolite: ", metabolite,
          " | imputation: ", i,
          # " | subgroup: ", subgroup,
          " | message: ", conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    )
    
    rm(dat_i)
    gc()
  }
  
  fit_mira <- as.mira(fits)
  
  pooled_sum <- summary(
    pool(fit_mira),
    conf.int = TRUE,
    exponentiate = TRUE
  )
  
  met_row <- pooled_sum[1, ]
  
  n_total <- median(sapply(fits, function(x) x$n))
  n_events <- median(sapply(fits, function(x) x$nevent))
  
  ph_tests <- lapply(fits, function(fit) {
    tryCatch(
      cox.zph(fit),
      error = function(e) NULL
    )
  })
  
  ph_metabolite_p <- sapply(ph_tests, function(x) {
    if (is.null(x)) NA_real_ else x$table[1, "p"]
  })
  
  ph_global_p <- sapply(ph_tests, function(x) {
    if (is.null(x)) NA_real_ else x$table["GLOBAL", "p"]
  })
  
  data.frame(
    Metabolite = metabolite,
    Subgroup = event_indicator,
    Hazard_Ratio = met_row$estimate,
    CI_Lower = met_row$`2.5 %`,
    CI_Upper = met_row$`97.5 %`,
    P_Value = met_row$p.value,
    N = n_total,
    N_events = n_events,
    PH_P_Global = median(ph_global_p, na.rm = TRUE),
    PH_P_Metabolite = median(ph_metabolite_p, na.rm = TRUE)
  )
}

# Include an event_indicator when create the CoxPH Function for all metabolite
run_full_analysis_mi <- function(metabolites, imp_mids, covariates_use, metabolite_data, event_indicator) {
  
  all_results <- lapply(metabolites, function(metabolite) {
    
    cat("Running:", metabolite, "/n")
    
    run_coxph_analysis_mi(
      metabolite = metabolite,
      imp_mids = imp_mids,
      covariates_use = covariates_use,
      metabolite_data = metabolite_data,
      event_indicator = event_indicator
    )
  })
  
  combined_results <- do.call(rbind, all_results)
  
  combined_results$P_Adjusted <- p.adjust(combined_results$P_Value, method = "BH")
  combined_results$PH_P_Global_Adjusted <- p.adjust(combined_results$PH_P_Global, method = "BH")
  combined_results$PH_P_Metabolite_Adjusted <- p.adjust(combined_results$PH_P_Metabolite, method = "BH")
  
  list(Adjusted = combined_results)
}

# Run the cause-specific CoxPH by histological subtypes
results_NSCLC <- run_full_analysis_mi(
  metabolites = metabolites,
  imp_mids = imp_mids,
  covariates_use = extended_covariates,
  metabolite_data = metabolite_data,
  event_indicator = "event_NSCLC"
)

results_SCLC <- run_full_analysis_mi(
  metabolites = metabolites,
  imp_mids = imp_mids,
  covariates_use = extended_covariates,
  metabolite_data = metabolite_data,
  event_indicator = "event_SCLC"
)

results_other <- run_full_analysis_mi(
  metabolites = metabolites,
  imp_mids = imp_mids,
  covariates_use = extended_covariates,
  metabolite_data = metabolite_data,
  event_indicator = "event_Other"
)

################ Run Complete-Case CoxPH Analysis ################
run_coxph_analysis <- function(metabolite, coxph_data, covariates_use) {
  
  # Create a subset dataset containing required variables for downstream analysis
  metabolite_df <- subset(coxph_data, select = c("eid",  "lungc_time", "lungc_status_coxph", all_of(covariates_use), metabolite))
  
  # Remove metabolite-specific missing
  metabolite_final <- metabolite_df %>%
    filter(!is.na(metabolite))
  
  # Construct the CoxPH model formula
  formula_adjusted <- as.formula(
    paste0(
      "Surv(lungc_time, lungc_status_coxph) ~ ",
      paste0("`", metabolite, "`"),
      " + ",
      paste0("`", covariates_use, "`", collapse = " + ")
    )
  )
  
  # Fit the CoxPH model using Efron method for tired events
  cox_adjusted <- coxph(formula_adjusted, data = metabolite_final, ties = "efron")
  
  # Check proportional hazards assumption based on the Schoenfeld residuals
  ph_test <- cox.zph(cox_adjusted)
  
  # Extract global and metabolite-specific p-values
  ph_global_p <- ph_test$table["GLOBAL", "p"]
  ph_metabolite_p <- ph_test$table[1, "p"]  # first term = the metabolite itself
  
  # Extract coefficient
  hr_coef <- coef(cox_adjusted)[1]
  
  # Convert the coefficient to a hazard ratio
  hr <- exp(hr_coef)
  
  # Calculate the 95% confidence interval on the hazard-ratio scale
  ci <- exp(confint(cox_adjusted)[1, ])
  
  # Extract the metabolite-specific p-value
  p_value <- summary(cox_adjusted)$coefficients[1, "Pr(>|z|)"]
  
  # Extract the analytic sample size and number of events
  n_total  <- cox_adjusted$n
  n_events <- cox_adjusted$nevent
  
  # Prepare a table of adjusted results
  result_table_adjusted <- data.frame(
    Metabolite = metabolite,
    Model = c("Complete-case"),
    Hazard_Ratio = hr,
    CI_Lower = ci[1],
    CI_Upper = ci[2],
    P_Value = p_value,
    N = n_total,
    N_events = n_events,
    PH_P_Global = ph_global_p,
    PH_P_Metabolite = ph_metabolite_p
  )  
}

run_full_analysis <- function(metabolites, coxph_data, covariates_use) {
  
  # Run the CoxPH analysis separately for each metabolite
  all_results <- lapply(metabolites, function(metabolite) {
    
    # Print the name of the metabolite currently being analyzed
    cat("Running:", metabolite, "/n")
    
    # Apply the single-metabolite CoxPH function
    run_coxph_analysis(
      metabolite = metabolite,
      coxph_data = coxph_data,
      covariates_use = covariates_use)
  })
  
  # Combine the metabolite-specific results into one data frame
  combined_results <- do.call(rbind, all_results)
  
  # Adjust p-values for multiple testing using the Benjamini-Hochberg false discovery rate method
  combined_results$P_Adjusted <- p.adjust(combined_results$P_Value, method = "BH")
  combined_results$PH_P_Global_Adjusted <- p.adjust(combined_results$PH_P_Global, method = "BH")
  combined_results$PH_P_Metabolite_Adjusted <- p.adjust(combined_results$PH_P_Metabolite, method = "BH")
  
  # Return the combined results
  list(Adjusted = combined_results)
}

# Run Complete-case CoxPH Analysis
results_complete_case <- run_full_analysis(
  metabolites = metabolites,
  coxph_data = coxph_data_complete_case,
  covariates_use = extended_covariates_original
)

################ Save CoxPH results ################
write_xlsx(results_overall, "coxph_results_overall.xlsx")
write_xlsx(results_female, "coxph_results_female.xlsx")
write_xlsx(results_male, "coxph_results_male.xlsx")
write_xlsx(results_current, "coxph_results_current.xlsx")
write_xlsx(results_previous, "coxph_results_former.xlsx")
write_xlsx(results_never, "coxph_results_never.xlsx")
write_xlsx(results_agegroup1, "coxph_results_agegroup1.xlsx")
write_xlsx(results_agegroup2, "coxph_results_agegroup2.xlsx")
write_xlsx(results_agegroup3, "coxph_results_agegroup3.xlsx")
write_xlsx(results_NSCLC, "coxph_results_NSCLC.xlsx")
write_xlsx(results_SCLC, "coxph_results_SCLC.xlsx")
write_xlsx(results_other, "coxph_results_other.xlsx")
write_xlsx(results_complete_case, "coxph_results_complete_case.xlsx")

