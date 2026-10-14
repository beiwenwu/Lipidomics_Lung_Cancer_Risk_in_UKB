library(dplyr)
library(mice)
library(writexl)

# Create the dataset used for multiple imputation by retaining the variables to be imputed, their predictors, and the survival outcome variables
impute_data <- coxph_data %>%
  select(
    eid,
    lungc_time,
    lungc_status_coxph,
    age, sex, ethnicity_UKB, education_ISCED, tdi,
    bmi, bmi2, waist_to_hip, 
    
    smk_status, smk_group, pack_year,
    smk_intensity, smk_age_start, smk_age_stop, smk_quit, 
    
    Summed_MET_minutes_weekly,
    MET_walking,MET_moderate,MET_vigorous,
    
    family_lc,lung_disease_history,lipid_lowering_meds,fasting_time,centre_group,
    
    age_group, event_NSCLC, event_SCLC, event_Other
  )

# Summarize missingness
impute_data %>%
  summarise(
    across(everything(),
                   list(n = ~sum(is.na(.)),
                        pct = ~mean(is.na(.))*100)))

################ Convert variables to proper types ################ 
impute_data$smk_group <- factor(impute_data$smk_group, levels = c("never", "current_cigarettes","current_cigar_pipe", "current_occasional","previous"))
impute_data$smk_quit <- factor(impute_data$smk_quit)
impute_data$lungc_status_coxph <- factor(impute_data$lungc_status_coxph)

################ Initialize MICE ################
ini   <- mice(impute_data, maxit = 0)
meth  <- ini$method
pred  <- ini$predictorMatrix
where <- is.na(impute_data)
post  <- ini$post


################ Specify imputation methods ################
meth[] <- "" #default: do not impute others
meth["smk_intensity"] <- "pmm"
meth["smk_age_start"] <- "pmm"
meth["smk_age_stop"]  <- "pmm"
meth["smk_quit"]      <- "logreg"
meth["MET_walking"] <- "pmm"     
meth["MET_moderate"] <- "pmm"     
meth["MET_vigorous"] <- "pmm" 

################ Custom predictors per variable ################
smk_vars <- c("smk_intensity", "smk_age_start", "smk_age_stop", "smk_quit")
met_vars <- c("MET_walking", "MET_moderate", "MET_vigorous")

# Initialize the predictor matrix
pred[,] <- 0

# Allow all available variables to be used as predictors for the smoking-related variables
pred[smk_vars, ] <- 1

# Allow all available variables to be used as predictors for the physical activity variables
pred[met_vars, ] <- 1

# no self-prediction
diag(pred) <- 0

# eid is not imputed, and not be used as predictor
pred[, "eid"] <- 0
pred["eid", ] <- 0

# smk_group is not imputed, and not be used as predictor
pred["smk_group", ] <- 0
pred[,"smk_group"] <- 0

# Do not impute smoking status
pred["smk_status", ] <- 0

# Do not use smoking status as a predictor for the smoking-related variables
pred[smk_vars,"smk_status" ] <- 0

# Retain smoking status as the only smoking-related predictor for physical activity variables, while exclude all other detailed smoking-related variables
pred[met_vars,smk_vars] <- 0
pred[met_vars,"smk_status"] <- 1

# Summed_MET_minutes_weekly is not imputed, and not be used as predictor
pred[, "Summed_MET_minutes_weekly"] <- 0
pred["Summed_MET_minutes_weekly",] <- 0

# pack_year is not imputed, and not be used as predictor
pred[, "pack_year"] <- 0
pred["pack_year",] <- 0

# age_group is not imputed, and not be used as predictor
pred[, "age_group"] <- 0
pred["age_group",] <- 0

# event_NSCLC is not imputed, and not be used as predictor
pred[, "event_NSCLC"] <- 0
pred["event_NSCLC",] <- 0

# event_SCLC is not imputed, and not be used as predictor
pred[, "event_SCLC"] <- 0
pred["event_SCLC",] <- 0

# event_Other is not imputed, and not be used as predictor
pred[, "event_Other"] <- 0
pred["event_Other",] <- 0

# Check predictor-matrix 
pred["smk_intensity",]
pred["smk_age_start",]
pred["smk_age_stop",]
pred["smk_quit",]
pred["MET_walking",]
pred["MET_moderate",]
pred["MET_vigorous",]
pred["pack_year",]

################ Run MICE ################
# When smk_group is never smokers, do not impute smoke-related variables 
where[impute_data$smk_group == "never",
      c("smk_intensity", "smk_age_start", "smk_age_stop", "smk_quit")] <- FALSE

# Run ten imputed datasets with twenty iterations
imp <- mice(
  impute_data,
  m = 10,
  maxit = 20,
  method = meth,
  predictorMatrix = pred,
  where = where,
  printFlag = TRUE,
  seed = 123
)

################ Check MICE results ################ 
imp$loggedEvents
imp$predictorMatrix["smk_intensity", ]
imp$predictorMatrix["smk_status", ]
imp$predictorMatrix["smk_group", ]
imp$predictorMatrix["MET_walking", ]

################ Derive pack-years after imputation ################
# Convert the MICE object to long format
completed_long <- complete(imp, action = "long", include = TRUE)

# Recode mk_age_start, smk_age_stop, smk_duration, pack_year to consider implausible values
completed_long_update <- completed_long %>%
  # Individuals who started and gave up smoking before 16 years of age, their smk_age_start were coded as NA
  mutate(
    smk_age_start_update = case_when(
      smk_age_start <= 16 & smk_age_stop <= 16 ~ NA_real_,
      TRUE ~ smk_age_start
    ),
    
    # Individuals who started and gave up smoking before 16 years of age, their smk_age_stop were coded as NA
    smk_age_stop_update = case_when(
      smk_age_start <= 16 & smk_age_stop <= 16 ~ NA_real_,
      TRUE ~ smk_age_stop
    ),
    
    smk_duration = case_when(
      # when smk_age_stop = smk_age_start, set smk_duration to 0
      smk_age_start_update == smk_age_stop_update ~ 0,
      # when smk_age_stop < smk_age_start, which is implausible, so set smk_duration to NA
      smk_age_stop_update < smk_age_start_update ~ NA_real_,
      # Individuals who gave up smoking for more than 6 months, subtract six months (0.5 year) from smk_duration
      smk_age_stop_update != smk_age_start_update & smk_quit == "Yes" ~ smk_age_stop_update - smk_age_start_update -0.5,
      smk_age_stop_update != smk_age_start_update & smk_quit == "No" ~ smk_age_stop_update - smk_age_start_update,
      smk_age_stop_update != smk_age_start_update & is.na(smk_quit) ~ smk_age_stop_update - smk_age_start_update,
      TRUE ~ NA_real_
    ),
    
    # Recalculate pack_year with newly imputed values
    pack_year_imputed = case_when(
      smk_group == "never" ~ 0,
      smk_duration == 0 ~ 0,
      smk_duration < 0 ~ NA_real_,
      smk_group %in% c("current_cigarettes", "current_cigar_pipe", "current_occasional", "previous") ~ smk_intensity / 20 * smk_duration,
      TRUE ~ NA_real_
    ),
    
    # Recalculate summed_MET_minutes_weekly with newly imputed values
    summed_MET_minutes_weekly_imputed = MET_walking + MET_moderate + MET_vigorous
  )

# Create numeric survival variables before as.mids()
completed_long_update <- completed_long_update %>%
  mutate(
    lungc_time = as.numeric(lungc_time),
    lungc_status_coxph = case_when(
      as.character(lungc_status_coxph) == "1" ~ 1,
      as.character(lungc_status_coxph) == "0" ~ 0,
      TRUE ~ NA_real_
    )
  )

# Reconstruct the multiply imputed dataset as a mice mids object for downstream pooled analyses using Rubin's rules
imp_mids <- as.mids(completed_long_update)

################ Save the MICE ################
setwd("M:/shared/Beiwen Wu/Project 1/UKB_RAP/UKB_Data/MICE")
saveRDS(imp, "imp_m10maxit20.rds")
saveRDS(completed_long_update, "completed_long_m10maxit20.rds")
saveRDS(imp_mids, "imp_mids_m10maxit20.rds")
saveRDS(pred, "predictorMatrix_m10maxit20.rds")
saveRDS(meth, "methods_m10maxit20.rds")
saveRDS(where, "where_m10maxit20.rds")

