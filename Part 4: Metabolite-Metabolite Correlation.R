library(dplyr)
library(tidyr)
library(purrr)
library(mice)

# Retrieve the ten imputed datasets
completed_data1 <- complete(imp_mids, 1)
completed_data2 <- complete(imp_mids, 2)
completed_data3 <- complete(imp_mids, 3)
completed_data4 <- complete(imp_mids, 4)
completed_data5 <- complete(imp_mids, 5)
completed_data6 <- complete(imp_mids, 6)
completed_data7 <- complete(imp_mids, 7)
completed_data8 <- complete(imp_mids, 8)
completed_data9 <- complete(imp_mids, 9)
completed_data10 <- complete(imp_mids, 10)

# Remove remaining missingess in covariates 
completed_data1_update <- completed_data1 %>%
  drop_na(all_of(extended_covariates))
completed_data2_update <- completed_data2 %>%
  drop_na(all_of(extended_covariates))
completed_data3_update <- completed_data3 %>%
  drop_na(all_of(extended_covariates))
completed_data4_update <- completed_data4 %>%
  drop_na(all_of(extended_covariates))
completed_data5_update <- completed_data5 %>%
  drop_na(all_of(extended_covariates))
completed_data6_update <- completed_data6 %>%
  drop_na(all_of(extended_covariates))
completed_data7_update <- completed_data7 %>%
  drop_na(all_of(extended_covariates))
completed_data8_update <- completed_data8 %>%
  drop_na(all_of(extended_covariates))
completed_data9_update <- completed_data9 %>%
  drop_na(all_of(extended_covariates))
completed_data10_update <- completed_data10 %>%
  drop_na(all_of(extended_covariates))

# Merge imputed datasets with metabolite data
completed_data1_update2 <- merge(completed_data1,metabolite_data , by="eid",sort=FALSE,all.x=FALSE)
completed_data2_update2 <- merge(completed_data2,metabolite_data , by="eid",sort=FALSE,all.x=FALSE)
completed_data3_update2 <- merge(completed_data3,metabolite_data , by="eid",sort=FALSE,all.x=FALSE)
completed_data4_update2 <- merge(completed_data4,metabolite_data , by="eid",sort=FALSE,all.x=FALSE)
completed_data5_update2 <- merge(completed_data5,metabolite_data , by="eid",sort=FALSE,all.x=FALSE)
completed_data6_update2 <- merge(completed_data6,metabolite_data , by="eid",sort=FALSE,all.x=FALSE)
completed_data7_update2 <- merge(completed_data7,metabolite_data , by="eid",sort=FALSE,all.x=FALSE)
completed_data8_update2 <- merge(completed_data8,metabolite_data , by="eid",sort=FALSE,all.x=FALSE)
completed_data9_update2 <- merge(completed_data9,metabolite_data , by="eid",sort=FALSE,all.x=FALSE)
completed_data10_update2 <- merge(completed_data10,metabolite_data , by="eid",sort=FALSE,all.x=FALSE)

# Combine the 10 completed imputed datasets into a named list
imputed_data <- list(
  imp1 = completed_data1_update2,
  imp2 = completed_data2_update2,
  imp3 = completed_data3_update2,
  imp4 = completed_data4_update2,
  imp5 = completed_data5_update2,
  imp6 = completed_data6_update2,
  imp7 = completed_data7_update2,
  imp8 = completed_data8_update2,
  imp9 = completed_data9_update2,
  imp10 = completed_data10_update2
)

# Calculate metabolite–metabolite correlation matrices across the imputed datasets and obtain an average correlation matrix
pool_correlations <- function(imputed_data, metabolites, method = "pearson") {
  
  # Calculate a correlation matrix for each metabolite separately within each imputed dataset
  cor_list <- map(
    imputed_data,
    function(dat) {
      dat %>%
        select(all_of(metabolites)) %>%
        cor(
          use = "pairwise.complete.obs", # Use all available observations for each metabolite pair
          method = method # Default to Pearson correlation
        )
    }
  )
  
  # Calculate the mean of each pairwise correlation across the 10 imputed datasets
  cor_avg <- (cor_list[[1]] + cor_list[[2]] + cor_list[[3]] + cor_list[[4]] + cor_list[[5]] +
                cor_list[[6]] + cor_list[[7]] + cor_list[[8]] + cor_list[[9]] + cor_list[[10]]
  ) / 10
  
  # Ensure that the diagonal elements represent perfect self-correlation
  diag(cor_avg) <- 1
  
  # Return the averaged metabolite correlation matrix
  cor_avg
}

################ Run correlation matrics in overall study population ################ 
# Calculate the pooled metabolite–metabolite correlation matrix
pooled_metabolite_cor <- pool_correlations(
  imputed_data = imputed_data,
  metabolites = metabolites,
  method = "pearson"
)

# Extract the correlation submatrix for the absolute metabolites
metabolite_absolute <-c("Total_C", "non_HDL_C", "Remnant_C", "VLDL_C", "Clinical_LDL_C", "LDL_C", "HDL_C", 
                        "Total_TG", "VLDL_TG", "LDL_TG", "HDL_TG", 
                        "Total_PL", "VLDL_PL", "LDL_PL", "HDL_PL", 
                        "Total_CE", "VLDL_CE", "LDL_CE", "HDL_CE", 
                        "Total_FC", "VLDL_FC", "LDL_FC", "HDL_FC", 
                        "Total_L", "VLDL_L", "LDL_L", "HDL_L", 
                        "Total_P", "VLDL_P", "LDL_P", "HDL_P", 
                        "VLDL_PS", "LDL_PS", "HDL_PS", 
                        "Phosphoglycerides", 
                        "Total_Cholines", "Phosphatidylcholines", "Sphingomyelins", 
                        "ApoB", "ApoA1", 
                        "Total_FA", "Unsaturation", "Omega_3", "Omega_6", "PUFA", "MUFA", "SFA", "LA", "DHA", 
                        "Alanine", "Glutamine", "Glycine", "Histidine", "Total_BCAA", "Isoleucine", "Leucine", "Valine", "Phenylalanine", "Tyrosine", 
                        "Glucose", "Lactate", "Pyruvate", "Citrate", 
                        "3_Hydroxybutyrate", "Acetate", "Acetoacetate", "Acetone", 
                        "Creatinine", "Albumin", 
                        "Glycoprotein_Acetyls", 
                        "XXL_VLDL_P", "XXL_VLDL_L", "XXL_VLDL_PL", "XXL_VLDL_C", "XXL_VLDL_CE", "XXL_VLDL_FC", "XXL_VLDL_TG", 
                        "XL_VLDL_P", "XL_VLDL_L", "XL_VLDL_PL", "XL_VLDL_C", "XL_VLDL_CE", "XL_VLDL_FC", "XL_VLDL_TG", 
                        "L_VLDL_P", "L_VLDL_L", "L_VLDL_PL", "L_VLDL_C", "L_VLDL_CE", "L_VLDL_FC", "L_VLDL_TG", 
                        "M_VLDL_P", "M_VLDL_L", "M_VLDL_PL", "M_VLDL_C", "M_VLDL_CE", "M_VLDL_FC", "M_VLDL_TG", 
                        "S_VLDL_P", "S_VLDL_L", "S_VLDL_PL", "S_VLDL_C", "S_VLDL_CE", "S_VLDL_FC", "S_VLDL_TG", 
                        "XS_VLDL_P", "XS_VLDL_L", "XS_VLDL_PL", "XS_VLDL_C", "XS_VLDL_CE", "XS_VLDL_FC", "XS_VLDL_TG", 
                        "IDL_P", "IDL_L", "IDL_PL", "IDL_C", "IDL_CE", "IDL_FC", "IDL_TG", 
                        "L_LDL_P", "L_LDL_L", "L_LDL_PL", "L_LDL_C", "L_LDL_CE", "L_LDL_FC", "L_LDL_TG", 
                        "M_LDL_P", "M_LDL_L", "M_LDL_PL", "M_LDL_C", "M_LDL_CE", "M_LDL_FC", "M_LDL_TG", 
                        "S_LDL_P", "S_LDL_L", "S_LDL_PL", "S_LDL_C", "S_LDL_CE", "S_LDL_FC", "S_LDL_TG", 
                        "XL_HDL_P", "XL_HDL_L", "XL_HDL_PL", "XL_HDL_C", "XL_HDL_CE", "XL_HDL_FC", "XL_HDL_TG", 
                        "L_HDL_P", "L_HDL_L", "L_HDL_PL", "L_HDL_C", "L_HDL_CE", "L_HDL_FC", "L_HDL_TG", 
                        "M_HDL_P", "M_HDL_L", "M_HDL_PL", "M_HDL_C", "M_HDL_CE", "M_HDL_FC", "M_HDL_TG", 
                        "S_HDL_P", "S_HDL_L", "S_HDL_PL", "S_HDL_C", "S_HDL_CE", "S_HDL_FC", "S_HDL_TG",
                        "Glucose_lactate", "Spectrometer_corrected_alanine")

pooled_metabolite_cor_abs <- pooled_metabolite_cor[
  metabolite_absolute,
  metabolite_absolute,
  drop = FALSE
]

# Extract the correlation submatrix for the relative metabolites
metabolite_ratio <-c("TG/PG", 
                     "ApoB/ApoA1", 
                     "Omega_3_%", "Omega_6_%", "PUFA_%", "MUFA_%", "SFA_%", "LA_%", "DHA_%", "PUFA/MUFA", "Omega_6/Omega_3", 
                     "XXL_VLDL_PL_%", "XXL_VLDL_C_%", "XXL_VLDL_CE_%", "XXL_VLDL_FC_%", "XXL_VLDL_TG_%", 
                     "XL_VLDL_PL_%", "XL_VLDL_C_%", "XL_VLDL_CE_%", "XL_VLDL_FC_%", "XL_VLDL_TG_%", 
                     "L_VLDL_PL_%", "L_VLDL_C_%", "L_VLDL_CE_%", "L_VLDL_FC_%", "L_VLDL_TG_%", 
                     "M_VLDL_PL_%", "M_VLDL_C_%", "M_VLDL_CE_%", "M_VLDL_FC_%", "M_VLDL_TG_%", 
                     "S_VLDL_PL_%", "S_VLDL_C_%", "S_VLDL_CE_%", "S_VLDL_FC_%", "S_VLDL_TG_%", 
                     "XS_VLDL_PL_%", "XS_VLDL_C_%", "XS_VLDL_CE_%", "XS_VLDL_FC_%", "XS_VLDL_TG_%", 
                     "IDL_PL_%", "IDL_C_%", "IDL_CE_%", "IDL_FC_%", "IDL_TG_%", 
                     "L_LDL_PL_%", "L_LDL_C_%", "L_LDL_CE_%", "L_LDL_FC_%", "L_LDL_TG_%", 
                     "M_LDL_PL_%", "M_LDL_C_%", "M_LDL_CE_%", "M_LDL_FC_%", "M_LDL_TG_%", 
                     "S_LDL_PL_%", "S_LDL_C_%", "S_LDL_CE_%", "S_LDL_FC_%", "S_LDL_TG_%", 
                     "XL_HDL_PL_%", "XL_HDL_C_%", "XL_HDL_CE_%", "XL_HDL_FC_%", "XL_HDL_TG_%", 
                     "L_HDL_PL_%", "L_HDL_C_%", "L_HDL_CE_%", "L_HDL_FC_%", "L_HDL_TG_%", 
                     "M_HDL_PL_%", "M_HDL_C_%", "M_HDL_CE_%", "M_HDL_FC_%", "M_HDL_TG_%", 
                     "S_HDL_PL_%", "S_HDL_C_%", "S_HDL_CE_%", "S_HDL_FC_%", "S_HDL_TG_%") 

pooled_metabolite_cor_ratio <- pooled_metabolite_cor[
  metabolite_ratio,
  metabolite_ratio,
  drop = FALSE
]

################ Run correlation matrics by sex ################ 
# Restrict each imputed dataset to female participants
imputed_data_female <- map(
  imputed_data,
  ~ .x %>% filter(sex %in% c("Female"))
)

# Calculate the pooled metabolite–metabolite correlation matrix in females
pooled_cor_female <- pool_correlations(
  imputed_data = imputed_data_female,
  metabolites = metabolites,
  method = "pearson"
)

# Extract the correlation submatrix for the absolute metabolites in females
pooled_cor_female_abs <- pooled_cor_female[
  metabolite_absolute,
  metabolite_absolute,
  drop = FALSE
]

# Extract the correlation submatrix for the relative metabolites in females
pooled_cor_female_ratio <- pooled_cor_female[
  metabolite_ratio,
  metabolite_ratio,
  drop = FALSE
]

# Restrict each imputed dataset to male participants
imputed_data_male <- map(
  imputed_data,
  ~ .x %>% filter(sex %in% c("Male"))
)

# Calculate the pooled metabolite–metabolite correlation matrix in males
pooled_cor_male <- pool_correlations(
  imputed_data = imputed_data_male,
  metabolites = metabolites,
  method = "pearson"
)

# Extract the correlation submatrix for the absolute metabolites in males
pooled_cor_male_abs <- pooled_cor_male[
  metabolite_absolute,
  metabolite_absolute,
  drop = FALSE
]

# Extract the correlation submatrix for the relative metabolites in males
pooled_cor_male_ratio <- pooled_cor_male[
  metabolite_ratio,
  metabolite_ratio,
  drop = FALSE
]

################ Run correlation matrics by age groups ################
# Restrict each imputed dataset to participants <=50 years
imputed_data_agegroup1 <- map(
  imputed_data,
  ~ .x %>% filter(age_group %in% c("<=50"))
)

# Calculate the pooled metabolite–metabolite correlation matrix in participants <=50 years
pooled_cor_agegroup1 <- pool_correlations(
  imputed_data = imputed_data_agegroup1,
  metabolites = metabolites,
  method = "pearson"
)

# Extract the correlation submatrix for the absolute metabolites in participants <=50 years
pooled_cor_agegroup1_abs <- pooled_cor_agegroup1[
  metabolite_absolute,
  metabolite_absolute,
  drop = FALSE
]

# Extract the correlation submatrix for the relative metabolites in participants <=50 years
pooled_cor_agegroup1_ratio <- pooled_cor_agegroup1[
  metabolite_ratio,
  metabolite_ratio,
  drop = FALSE
]

# Restrict each imputed dataset to participants >50 and <=60 years
imputed_data_agegroup2 <- map(
  imputed_data,
  ~ .x %>% filter(age_group %in% c("(50-60]"))
)

# Calculate the pooled metabolite–metabolite correlation matrix in participants >50 and <=60 years
pooled_cor_agegroup2 <- pool_correlations(
  imputed_data = imputed_data_agegroup2,
  metabolites = metabolites,
  method = "pearson"
)

# Extract the correlation submatrix for the absolute metabolites in participants >50 and <=60 years
pooled_cor_agegroup2_abs <- pooled_cor_agegroup2[
  metabolite_absolute,
  metabolite_absolute,
  drop = FALSE
]

# Extract the correlation submatrix for the relative metabolites in participants >50 and <=60 years
pooled_cor_agegroup2_ratio <- pooled_cor_agegroup2[
  metabolite_ratio,
  metabolite_ratio,
  drop = FALSE
]

# Restrict each imputed dataset to participants >60 years
imputed_data_agegroup3 <- map(
  imputed_data,
  ~ .x %>% filter(age_group %in% c(">60"))
)

# Calculate the pooled metabolite–metabolite correlation matrix in participants >60 years
pooled_cor_agegroup3 <- pool_correlations(
  imputed_data = imputed_data_agegroup3,
  metabolites = metabolites,
  method = "pearson"
)

# Extract the correlation submatrix for the absolute metabolites in participants >60 years
pooled_cor_agegroup3_abs <- pooled_cor_agegroup3[
  metabolite_absolute,
  metabolite_absolute,
  drop = FALSE
]

# Extract the correlation submatrix for the relative metabolites in participants >60 years
pooled_cor_agegroup3_ratio <- pooled_cor_agegroup3[
  metabolite_ratio,
  metabolite_ratio,
  drop = FALSE
]

################ Run correlation matrics by smoking groups ################
# Restrict each imputed dataset to participants who are current smokers
imputed_data_current <- map(
  imputed_data,
  ~ .x %>% filter(smk_status %in% c("Current"))
)

# Calculate the pooled metabolite–metabolite correlation matrix in participants who are current smokers
pooled_cor_current <- pool_correlations(
  imputed_data = imputed_data_current,
  metabolites = metabolites,
  method = "pearson"
)

# Extract the correlation submatrix for the absolute metabolites in participants who are current smokers
pooled_cor_current_abs <- pooled_cor_current[
  metabolite_absolute,
  metabolite_absolute,
  drop = FALSE
]

# Extract the correlation submatrix for the relative metabolites in participants who are current smokers
pooled_cor_current_ratio <- pooled_cor_current[
  metabolite_ratio,
  metabolite_ratio,
  drop = FALSE
]

# Restrict each imputed dataset to participants who are former smokers
imputed_data_former <- map(
  imputed_data,
  ~ .x %>% filter(smk_status %in% c("Previous"))
)

# Calculate the pooled metabolite–metabolite correlation matrix in participants who are former smokers
pooled_cor_former <- pool_correlations(
  imputed_data = imputed_data_former,
  metabolites = metabolites,
  method = "pearson"
)

# Extract the correlation submatrix for the absolute metabolites in participants who are former smokers
pooled_cor_former_abs <- pooled_cor_former[
  metabolite_absolute,
  metabolite_absolute,
  drop = FALSE
]

# Extract the correlation submatrix for the relative metabolites in participants who are former smokers
pooled_cor_former_ratio <- pooled_cor_former[
  metabolite_ratio,
  metabolite_ratio,
  drop = FALSE
]

# Restrict each imputed dataset to participants who are never smokers
imputed_data_never <- map(
  imputed_data,
  ~ .x %>% filter(smk_status %in% c("Never"))
)

# Calculate the pooled metabolite–metabolite correlation matrix in participants who are never smokers
pooled_cor_never <- pool_correlations(
  imputed_data = imputed_data_never,
  metabolites = metabolites,
  method = "pearson"
)

# Extract the correlation submatrix for the absolute metabolites in participants who are never smokers
pooled_cor_never_abs <- pooled_cor_never[
  metabolite_absolute,
  metabolite_absolute,
  drop = FALSE
]

# Extract the correlation submatrix for the relative metabolites in participants who are never smokers
pooled_cor_never_ratio <- pooled_cor_never[
  metabolite_ratio,
  metabolite_ratio,
  drop = FALSE
]

################ Run correlation matrics by histology ################
## Because the histology-specific analyses used the same full cohort with different event indicators, 
## then the metabolite–metabolite correlation matrix should usually be calculated using the same full analytic cohort, 
## not separately for NSCLC/SCLC/Others.

################ Run correlation matrics in complete cases ################
# Load the complete case data
coxph_data_complete_case 

metabolite_cor_complete_case <-  coxph_data_complete_case %>%
      select(all_of(metabolites)) %>%
      cor(
        use = "pairwise.complete.obs", # Use all available observations for each metabolite pair
        method = "pearson" # Default to Pearson correlation
      )


# Extract the correlation submatrix for the absolute metabolites in complete-cases
metabolite_cor_complete_case_abs <- metabolite_cor_complete_case[
  metabolite_absolute,
  metabolite_absolute,
  drop = FALSE
]

# Extract the correlation submatrix for the relative metabolites in complete-cases
metabolite_cor_complete_case_ratio <- metabolite_cor_complete_case[
  metabolite_ratio,
  metabolite_ratio,
  drop = FALSE
]

################ Save correlation matrics ################
setwd('/Users/beiwen/Library/CloudStorage/OneDrive-UniversityofToronto/PhD Documents/PhD Projects/Objective 1/UKB RAP/Results')
saveRDS(imputed_data,"imputed_data.rds")
saveRDS(pooled_metabolite_cor, "pooled_metabolite_cor.rds")
saveRDS(pooled_metabolite_cor_abs, "pooled_metabolite_cor_abs.rds")
saveRDS(pooled_metabolite_cor_ratio, "pooled_metabolite_cor_ratio.rds")
saveRDS(pooled_cor_female, "pooled_cor_female.rds")
saveRDS(pooled_cor_female_abs, "pooled_cor_female_abs.rds")
saveRDS(pooled_cor_female_ratio, "pooled_cor_female_ratio.rds")
saveRDS(pooled_cor_male, "pooled_cor_male.rds")
saveRDS(pooled_cor_male_abs, "pooled_cor_male_abs.rds")
saveRDS(pooled_cor_male_ratio, "pooled_cor_male_ratio.rds")
saveRDS(pooled_cor_agegroup1, "pooled_cor_agegroup1.rds")
saveRDS(pooled_cor_agegroup1_abs, "pooled_cor_agegroup1_abs.rds")
saveRDS(pooled_cor_agegroup1_ratio, "pooled_cor_agegroup1_ratio.rds")
saveRDS(pooled_cor_agegroup2, "pooled_cor_agegroup2.rds")
saveRDS(pooled_cor_agegroup2_abs, "pooled_cor_agegroup2_abs.rds")
saveRDS(pooled_cor_agegroup2_ratio, "pooled_cor_agegroup2_ratio.rds")
saveRDS(pooled_cor_agegroup3, "pooled_cor_agegroup3.rds")
saveRDS(pooled_cor_agegroup3_abs, "pooled_cor_agegroup3_abs.rds")
saveRDS(pooled_cor_agegroup3_ratio, "pooled_cor_agegroup3_ratio.rds")
saveRDS(pooled_cor_current, "pooled_cor_current.rds")
saveRDS(pooled_cor_current_abs, "pooled_cor_current_abs.rds")
saveRDS(pooled_cor_current_ratio, "pooled_cor_current_ratio.rds")
saveRDS(pooled_cor_former, "pooled_cor_former.rds")
saveRDS(pooled_cor_former_abs, "pooled_cor_former_abs.rds")
saveRDS(pooled_cor_former_ratio, "pooled_cor_former_ratio.rds")
saveRDS(pooled_cor_never, "pooled_cor_never.rds")
saveRDS(pooled_cor_never_abs, "pooled_cor_never_abs.rds")
saveRDS(pooled_cor_never_ratio, "pooled_cor_never_ratio.rds")
saveRDS(metabolite_cor_complete_case, "metabolite_cor_complete_case.rds")
saveRDS(metabolite_cor_complete_case_abs, "metabolite_cor_complete_case_abs.rds")
saveRDS(metabolite_cor_complete_case_ratio, "metabolite_cor_complete_case_ratio.rds")