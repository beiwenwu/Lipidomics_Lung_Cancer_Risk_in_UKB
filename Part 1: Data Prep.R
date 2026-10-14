# clear environment
rm(list = ls())
# clear console
cat("\f")
# close any graphics
if(!is.null(dev.list())) dev.off()

library(dplyr)
library(tidyr)
library(rspiro)
library(stringr)
library(arrow)
library(plyr)
library(reshape2)
library(data.table)
library(writexl)

# Load data
metabolite_data <- read.csv("M:/shared/Beiwen Wu/Project 1/UKB_RAP/UKB_Data/Original/metabolite_data_full_20251210.csv", check.names=FALSE)
cancer_data <- read.csv("M:/shared/Beiwen Wu/Project 1/UKB_RAP/UKB_Data/Original/cancer_data_full_20251210.csv")
covariates_data <- read.csv("M:/shared/Beiwen Wu/Project 1/UKB_RAP/UKB_Data/Original/covariate_data_full_20251210.csv")
withdrawn_ids  <- read.csv("M:/shared/Beiwen Wu/Project 1/UKB_RAP/UKB_Data/Original/w87201_20260310.csv")$X1005970

######################### Covariates data prep ######################### 
# Remove those withdraw from the study
covariates_data <- covariates_data %>%
  filter(!(eid %in% withdrawn_ids))

################ Age ################ 
colnames(covariates_data)[colnames(covariates_data) == "p21022"] <- "age"

# Create age groups
covariates_data <- covariates_data %>%
  mutate(
    age_group = case_when(
      age <= 50 ~ "<=50",
      age > 50 & age <= 60 ~ "(50-60]",
      age > 60 ~ ">60",
      TRUE ~ NA_character_
    ),
    age_group = factor(age_group, levels = c("<=50", "(50-60]", ">60")))

################ Sex ################ 
colnames(covariates_data)[colnames(covariates_data) == "p31"] <- "sex"
covariates_data$sex <- as.factor(covariates_data$sex)

################ Ethnic background ################ 
colnames(covariates_data)[colnames(covariates_data) == "p21000_i0"] <- "ethnicity"
covariates_data$ethnicity <- as.factor(covariates_data$ethnicity)

# Create ethnicity groups
covariates_data <- covariates_data %>%
  mutate(
    ethnicity_UKB = case_when(
      ethnicity %in% c("White", "British", "Irish", "Any other white background") ~ "White",
      ethnicity %in% c("Mixed", "White and Black Caribbean", "White and Black African", 
                       "White and Asian", "Any other mixed background", "Other ethnic group") ~ "Mixed/Others",
      ethnicity %in% c("Asian or Asian British", "Indian", "Pakistani", "Bangladeshi", 
                       "Any other Asian background", "Chinese") ~ "Asian or Asian British",
      ethnicity %in% c("Black or Black British", "Caribbean", "African", "Any other Black background") ~ "Black or Black British",
      # ethnicity == "Other ethnic group" ~ "Other ethnic group",
      TRUE ~ NA_character_ #Empty, Do not know, Prefer not to answer are treated as NA
    ),
    ethnicity_UKB = factor(ethnicity_UKB, levels=c("White", "Black or Black British", "Asian or Asian British","Mixed/Others"))
  )

################ Education ################
colnames(covariates_data)[colnames(covariates_data) == "p6138_i0"] <- "education"
covariates_data$education <- as.factor(covariates_data$education)

# Define education group by mapping UKB qualifications to International Standard Classification of Education (ISCED) scale 
covariates_data <- covariates_data %>%
  mutate(
    education_ISCED = case_when(
      str_detect(education, "College or University degree") ~ "college",
      str_detect(education, "NVQ or HND or HNC or equivalent") ~ "some college",
      str_detect(education, "A levels/AS levels or equivalent") | str_detect(education, "Other professional") ~ "post high school",
      str_detect(education, "O levels") | str_detect(education, "CSEs") ~ "high school",
      education == "None of the above" ~ "less than high school",
      TRUE ~ NA_character_
    ),
    education_ISCED = factor(education_ISCED, levels = c(
      "college", "some college", "post high school", "high school", "less than high school"
    ))
  )

################ Townsend deprivation index ################ 
colnames(covariates_data)[colnames(covariates_data) == "p22189"] <- "tdi"

################ Body size measures ################ 
colnames(covariates_data)[colnames(covariates_data) == "p21001_i0"] <- "bmi"
colnames(covariates_data)[colnames(covariates_data) == "p48_i0"] <- "waist_circumference"
colnames(covariates_data)[colnames(covariates_data) == "p49_i0"] <- "hip_circumference"
colnames(covariates_data)[colnames(covariates_data) == "p50_i0"] <- "height"

# Create BMI^2 term
covariates_data <- covariates_data %>%
  mutate(
    bmi2 = bmi^2)

# Calculate waist to hip ratio
covariates_data$waist_to_hip <- covariates_data$waist_circumference / covariates_data$hip_circumference

################ Smoking  ################ 
covariates_data$p1239_i0 <- factor(covariates_data$p1239_i0, 
                                   levels = c("Yes, on most or all days", "Only occasionally", "No", "Prefer not to answer"))

covariates_data$p1249_i0 <- factor(covariates_data$p1249_i0, 
                                   levels = c("Smoked on most or all days", "Smoked occasionally", "Just tried once or twice", "I have never smoked","Prefer not to answer"))

covariates_data$p2644_i0 <- factor(covariates_data$p2644_i0,levels=c("Yes","No","Do not know","Prefer not to answer"))

covariates_data$p3436_i0 <- as.numeric(covariates_data$p3436_i0) # Empty, Do not know, and Prefer not to answer are treated as NA
covariates_data$p2867_i0 <- as.numeric(covariates_data$p2867_i0)
covariates_data$p2897_i0 <- as.numeric(covariates_data$p2897_i0)
covariates_data$p6194_i0 <- as.numeric(covariates_data$p6194_i0)
covariates_data$p3456_i0 <- as.numeric(covariates_data$p3456_i0)
covariates_data$p2887_i0 <- as.numeric(covariates_data$p2887_i0)
covariates_data$p6183_i0 <- as.numeric(covariates_data$p6183_i0)
covariates_data$p20116_i0 <- ifelse(covariates_data$p20116_i0 == "", NA, covariates_data$p20116_i0)

# Define smoking status
covariates_data <- covariates_data %>%
  mutate(
    smk_status = case_when(
      p1239_i0 %in% c("Yes, on most or all days", "Only occasionally") ~ "Current",
      p1239_i0 == "No" & p1249_i0 == "Smoked on most or all days" ~ "Previous",
      p1239_i0 == "No" & p1249_i0  %in% c("Just tried once or twice","Smoked occasionally") & p2644_i0 == "Yes" ~ "Previous",
      p1239_i0 == "No" & p1249_i0 %in% c("Just tried once or twice","Smoked occasionally") & p2644_i0 == "No" ~ "Never",
      p1239_i0 == "No" & p1249_i0 == "I have never smoked" ~ "Never",
      # is.na(p1239_i0) ~ NA_character_,
      TRUE ~ NA_character_
    ),
    smk_status = factor(smk_status, levels = c("Never", "Current", "Previous")
    ))

# Recalculate pack-year of smoking based on smoking groups
covariates_data <- covariates_data %>%
  mutate(pack_year = case_when(
    # Current daily cigarette smokers
    smk_status == "Current" & p1239_i0 == "Yes, on most or all days" & p3446_i0 %in% c("Manufactured cigarettes", "Hand-rolled cigarettes") ~ 
      ifelse(age == p3436_i0, 0,
             ifelse(p3486_i0 == "Yes, tried and stopped for at least 6 months", p3456_i0/20 * (age - p3436_i0 - 0.5), p3456_i0/20 * (age - p3436_i0))),
    
    # Current daily cigar/pipe smokers
    smk_status == "Current" & p1239_i0 == "Yes, on most or all days" & p3446_i0 == "Cigars or pipes" ~ 
      ifelse(p3436_i0 <= 16 & p6194_i0 <= 16, NA,
             ifelse(p3436_i0 == p6194_i0, 0,
                    ifelse(p3486_i0 == "Yes, tried and stopped for at least 6 months", p6183_i0/20 * (p6194_i0 - p3436_i0 - 0.5), p6183_i0/20 * (p6194_i0 - p3436_i0)))),
    
    # Current occasional smokers
    smk_status == "Current" & p1239_i0 == "Only occasionally" ~ 
      ifelse(p2897_i0 <= 16 & p2867_i0 <= 16, NA,
             ifelse(p2897_i0 == p2867_i0, 0,
                    ifelse(p2907_i0 == "Yes", p2887_i0/20 * (p2897_i0 - p2867_i0 - 0.5), p2887_i0/20 * (p2897_i0 - p2867_i0)))),
    
    # Previous smokers
    smk_status == "Previous" ~ 
      ifelse(p2897_i0 <= 16 & p2867_i0 <= 16, NA,
             ifelse(p2897_i0 == p2867_i0, 0,
                    ifelse(p2907_i0 == "Yes", p2887_i0/20 * (p2897_i0 - p2867_i0 - 0.5), p2887_i0/20 * (p2897_i0 - p2867_i0)))),
    
    # Never smokers
    smk_status == "Never" ~ 0
  ))

# Define smk_group, smk_intensity, smk_age_start, smk_age_stop, smk_quit
covariates_data <- covariates_data %>%
  mutate(
    smk_group = case_when(
      smk_status == "Current" &
        p1239_i0 == "Yes, on most or all days" &
        p3446_i0 %in% c("Manufactured cigarettes", "Hand-rolled cigarettes") ~ "current_cigarettes",

      smk_status == "Current" &
        p1239_i0 == "Yes, on most or all days" &
        p3446_i0 == "Cigars or pipes" ~ "current_cigar_pipe",

      smk_status == "Current" &
        p1239_i0 == "Only occasionally" ~ "current_occasional",

      smk_status == "Previous" ~ "previous",

      smk_status == "Never" ~ "never",

      TRUE ~ NA_character_
    ),

    smk_intensity = case_when(
      smk_group == "current_cigarettes" ~ p3456_i0,
      smk_group == "current_cigar_pipe" ~ p6183_i0,
      smk_group %in% c("current_occasional", "previous") ~ p2887_i0,
      TRUE ~ NA_real_
    ),

    smk_age_start = case_when(
      smk_group %in% c("current_cigarettes", "current_cigar_pipe") ~ p3436_i0,
      smk_group %in% c("current_occasional", "previous") ~ p2867_i0,
      TRUE ~ NA_real_
    ),

    smk_age_stop = case_when(
      smk_group == "current_cigarettes" ~ age,
      smk_group == "current_cigar_pipe" ~ p6194_i0,
      smk_group %in% c("current_occasional", "previous") ~ p2897_i0,
      TRUE ~ NA_real_
    ),

    smk_quit = case_when(
      smk_group %in% c("current_cigarettes", "current_cigar_pipe") ~ p3486_i0,
      smk_group %in% c("current_occasional", "previous") ~ p2907_i0,
      TRUE ~ NA_character_
    ),

    smk_quit = case_when(
      smk_quit %in% c("Yes", "Yes, tried and stopped for at least 6 months") ~ "Yes",
      smk_quit %in% c("No", "Yes, tried but was not able to stop or stopped for less than 6 months") ~ "No",
      TRUE ~ NA_character_
    )
  )

covariates_data$smk_group <- factor(covariates_data$smk_group, levels = c("never", "current_cigarettes","current_cigar_pipe", "current_occasional","previous"))
covariates_data$smk_quit <- factor(covariates_data$smk_quit)

################ Physical activity ################ 
colnames(covariates_data)[colnames(covariates_data) == "p22037_i0"] <- "MET_walking"
colnames(covariates_data)[colnames(covariates_data) == "p22038_i0"] <- "MET_moderate"
colnames(covariates_data)[colnames(covariates_data) == "p22039_i0"] <- "MET_vigorous"
colnames(covariates_data)[colnames(covariates_data) == "p22040_i0"] <- "Summed_MET_minutes_weekly"

################ Fasting Time ################ 
colnames(covariates_data)[colnames(covariates_data) == "p74_i0"] <- "fasting_time"

################ Assessment Centre ################ 
covariates_data$p54_i0 <- as.factor(covariates_data$p54_i0)

# Group UKB assessment centres based on location proximity 
covariates_data <- covariates_data %>%
  mutate(
    centre_group = case_when(
      p54_i0 %in% c("Barts","Croydon","Hounslow","Reading","Oxford","Bristol") ~ "England_South",
      p54_i0 %in% c("Birmingham","Nottingham","Stoke")               ~ "England_Midlands",
      p54_i0 %in% c("Manchester","Liverpool", "Bury","Stockport (pilot)") ~ "England_NorthWest",
      p54_i0 %in% c("Leeds","Sheffield")                             ~ "England_Yorkshire",
      p54_i0 %in% c("Middlesbrough","Newcastle")                    ~ "England_NorthEast",
      p54_i0 %in% c("Cardiff","Swansea","Wrexham")                   ~ "Wales",
      p54_i0 %in% c("Edinburgh","Glasgow")                           ~ "Scotland",
      TRUE                                                      ~ NA_character_
    ),
    centre_group = factor(centre_group)
  )

################ Family history of lung cancer  ################ 
vars <- c("p20107_i0","p20110_i0","p20111_i0")

covariates_data <- covariates_data %>%
  mutate(across(all_of(vars), ~ case_when(
    .x == "" ~ NA,
    TRUE ~ .x
  )))

# Father's LC history
covariates_data <- covariates_data %>%
  mutate(
    father_lc = case_when(
      str_detect(p20107_i0, "Lung cancer") ~ 1,
      !str_detect(p20107_i0, fixed("Do not know (group 1)")) & !str_detect(p20107_i0, fixed("Do not know (group 2)")) & !str_detect(p20107_i0, fixed("Prefer not to answer (group 1)")) & !str_detect(p20107_i0, fixed("Prefer not to answer (group 2)")) ~ 0,
      TRUE ~ NA_real_
    ),
    father_lc = factor(father_lc, levels = c(0,1)
    ))

# Mother's LC history
covariates_data <- covariates_data %>%
  mutate(
    mother_lc = case_when(
      str_detect(p20110_i0, "Lung cancer") ~ 1,
      !str_detect(p20110_i0, fixed("Do not know (group 1)")) & !str_detect(p20110_i0, fixed("Do not know (group 2)")) & !str_detect(p20110_i0, fixed("Prefer not to answer (group 1)")) & !str_detect(p20110_i0, fixed("Prefer not to answer (group 2)")) ~ 0,
      TRUE ~ NA_real_
    ),
    mother_lc = factor(mother_lc, levels = c(0,1)
    ))

#Sibling's LC history
covariates_data <- covariates_data %>%
  mutate(
    sibling_lc = case_when(
      str_detect(p20111_i0, "Lung cancer") ~ 1,
      !str_detect(p20111_i0, fixed("Do not know (group 1)")) & !str_detect(p20111_i0, fixed("Do not know (group 2)")) & !str_detect(p20111_i0, fixed("Prefer not to answer (group 1)")) & !str_detect(p20111_i0, fixed("Prefer not to answer (group 2)")) ~ 0,
      TRUE ~ NA_real_
    ),
    sibling_lc = factor(sibling_lc, levels = c(0,1)
    ))

# Define family lung cancer history based on father, mother, and sibling's histories
covariates_data <- covariates_data %>%
  mutate(
    family_lc = case_when(
      father_lc == 1 | mother_lc == 1 | sibling_lc == 1 ~ 1,
      father_lc == 0 & mother_lc == 0 & sibling_lc == 0 ~ 0,
      is.na(father_lc) & is.na(mother_lc) & is.na(sibling_lc) ~ NA_real_,
      TRUE ~ 0   # mixture of 0 and NA → assume no known history
    ),
    family_lc = factor(family_lc, levels = c(0,1))
  )

################ Medication: Lipid-lowering drugs  ################ 
vars <- c("p6177_i0", "p6153_i0", "p20003_i0_a0", "p20003_i0_a1", "p20003_i0_a2", "p20003_i0_a3", "p20003_i0_a4", "p20003_i0_a5",
          "p20003_i0_a6", "p20003_i0_a7", "p20003_i0_a8", "p20003_i0_a9", "p20003_i0_a10", "p20003_i0_a11", "p20003_i0_a12", "p20003_i0_a13", "p20003_i0_a14",
          "p20003_i0_a15", "p20003_i0_a16", "p20003_i0_a17", "p20003_i0_a18", "p20003_i0_a19", "p20003_i0_a20", "p20003_i0_a21", "p20003_i0_a22", "p20003_i0_a23",
          "p20003_i0_a24", "p20003_i0_a25", "p20003_i0_a26", "p20003_i0_a27", "p20003_i0_a28", "p20003_i0_a29", "p20003_i0_a30", "p20003_i0_a31", "p20003_i0_a32",
          "p20003_i0_a33", "p20003_i0_a34", "p20003_i0_a35", "p20003_i0_a36", "p20003_i0_a37", "p20003_i0_a38", "p20003_i0_a39", "p20003_i0_a40", "p20003_i0_a41",
          "p20003_i0_a42", "p20003_i0_a43", "p20003_i0_a44", "p20003_i0_a45", "p20003_i0_a46", "p20003_i0_a47")

covariates_data <- covariates_data %>%
  mutate(across(all_of(vars), ~ case_when(
    .x == "" ~ NA,
    TRUE ~ .x
  )))

# Touchscreen questionnaire
covariates_data <- covariates_data %>%
  mutate(
    lipid_lowering_meds_touchscreen = case_when(
      str_detect(p6177_i0, fixed("Cholesterol lowering medication")) ~ "Yes",
      str_detect(p6153_i0, fixed("Cholesterol lowering medication")) ~ "Yes",
      str_detect(p6177_i0, fixed("Do not know")) ~ NA_character_,
      str_detect(p6153_i0, fixed("Do not know")) ~ NA_character_,
      str_detect(p6177_i0, fixed("Prefer not to answer")) ~ NA_character_,
      str_detect(p6153_i0, fixed("Prefer not to answer")) ~ NA_character_,
      is.na(p6177_i0) & is.na(p6153_i0) ~ NA_character_,
      TRUE ~ "No"
    )
  )

# Verbal interview
lipid_lowering_meds <- c(
  # Statins (generic + brand)
  "atorvastatin", "lipitor",
  "rosuvastatin", "crestor",
  "simvastatin", "zocor", "simvador",
  "pravastatin", "lipostat",
  "fluvastatin", "lescol",
  "eptastatin",
  
  # Fibrates
  "fenofibrate", "gemfibrozil", "bezafibrate", "clofibrate", "ciprofibrate",
  "bezalip","modalim", 
  
  # Non-statin LDL-lowering
  "ezetimibe","ezetrol",
  "cholestyramine", "colestipol","questran","colestid",
  
  # Omega-3
  "omacor",
  
  # UKB NONSPECIFIC CODE
  "lipid lowering drug"
  
  #supplementt
  # "omega-3/fish oil supplement", "niacin", "nicotinic acid", "inositol nicotinate"
)

vars <- c("p20003_i0_a0", "p20003_i0_a1", "p20003_i0_a2", "p20003_i0_a3", "p20003_i0_a4", "p20003_i0_a5",
          "p20003_i0_a6", "p20003_i0_a7", "p20003_i0_a8", "p20003_i0_a9", "p20003_i0_a10", "p20003_i0_a11", "p20003_i0_a12", "p20003_i0_a13", "p20003_i0_a14",
          "p20003_i0_a15", "p20003_i0_a16", "p20003_i0_a17", "p20003_i0_a18", "p20003_i0_a19", "p20003_i0_a20", "p20003_i0_a21", "p20003_i0_a22", "p20003_i0_a23",
          "p20003_i0_a24", "p20003_i0_a25", "p20003_i0_a26", "p20003_i0_a27", "p20003_i0_a28", "p20003_i0_a29", "p20003_i0_a30", "p20003_i0_a31", "p20003_i0_a32",
          "p20003_i0_a33", "p20003_i0_a34", "p20003_i0_a35", "p20003_i0_a36", "p20003_i0_a37", "p20003_i0_a38", "p20003_i0_a39", "p20003_i0_a40", "p20003_i0_a41",
          "p20003_i0_a42", "p20003_i0_a43", "p20003_i0_a44", "p20003_i0_a45", "p20003_i0_a46", "p20003_i0_a47")


covariates_data <- covariates_data %>%
  mutate(
    lipid_lowering_meds_verbal = case_when(
      if_any(all_of(vars), ~ str_detect(.x, str_c(lipid_lowering_meds, collapse = "|"))) ~ "Yes",
      if_all(all_of(vars), ~ is.na(.x)) ~ NA_character_,
      TRUE ~ "No"
    )
  )

# Define the use of lipid lowering drugs from touchscreen and verbal interview
covariates_data <- covariates_data %>%
  mutate(
    lipid_lowering_meds = case_when(
      lipid_lowering_meds_touchscreen == "Yes" | lipid_lowering_meds_verbal == "Yes" ~ 1,
      lipid_lowering_meds_touchscreen == "No" & lipid_lowering_meds_verbal == "No" ~ 0,
      is.na(lipid_lowering_meds_touchscreen) & is.na(lipid_lowering_meds_verbal) ~ NA_real_,
      TRUE ~ 0   # mixture of 0 and NA → assume no known use
    ),
    lipid_lowering_meds = factor(lipid_lowering_meds, levels = c(0,1)
    ))

################ Lung function ################
vars <- c("p20002_i0_a0", "p20002_i0_a1", "p20002_i0_a2", "p20002_i0_a3", "p20002_i0_a4", "p20002_i0_a5", "p20002_i0_a6", "p20002_i0_a7", 
          "p20002_i0_a8", "p20002_i0_a9", "p20002_i0_a10", "p20002_i0_a11", "p20002_i0_a12", "p20002_i0_a13", "p20002_i0_a14", "p20002_i0_a15", 
          "p20002_i0_a16", "p20002_i0_a17", "p20002_i0_a18", "p20002_i0_a19", "p20002_i0_a20", "p20002_i0_a21", "p20002_i0_a22", "p20002_i0_a23", 
          "p20002_i0_a24", "p20002_i0_a25", "p20002_i0_a26", "p20002_i0_a27", "p20002_i0_a28", "p20002_i0_a29", "p20002_i0_a30", "p20002_i0_a31", 
          "p20002_i0_a32", "p20002_i0_a33")

# Verbal interview
lung_conditions_verbal <- c(
  "asthma",
  "chronic obstructive airways disease/copd",
  "emphysema/chronic bronchitis",
  "bronchitis",
  "emphysema",
  "alpha-1 antitrypsin deficiency")

covariates_data <- covariates_data %>%
  mutate(
    lung_disease_verbal = case_when(
      if_any(all_of(vars), ~ str_detect(.x, str_c(lung_conditions_verbal, collapse = "|"))) ~ "Yes",
      if_all(all_of(vars), ~ is.na(.x)) ~ NA_character_,
      TRUE ~ "No"
    )
  )

# Touchscreen questionnaire
vars <- c("p6150_i0","p6152_i0")

lung_conditions_touchscreen <- c("Emphysema/chronic bronchitis","Asthma")

covariates_data <- covariates_data %>%
  mutate(
    lung_disease_touchscreen = case_when(
      str_detect(p6152_i0, str_c(lung_conditions_touchscreen, collapse = "|")) ~ "Yes",
      str_detect(p6152_i0, "Prefer not to answer") ~ NA_character_,
      is.na(p6152_i0) ~ NA_character_,
      TRUE ~ "No"
    )
  )

# Spirometry measurement
vars <- c("p20031_i0_a0", "p20031_i0_a1", "p20031_i0_a2")

# Recode no flag to universal acceptance (UA)
covariates_data <- covariates_data %>%
  mutate(across(all_of(vars), ~ case_when(
    .x == "" & p3137_i0 != 0 ~ "UA",
    .x == "" & is.na(p3137_i0) ~ "NO BLOWS",
    TRUE ~ .x
  )))

# Flag blows with acceptable start
covariates_data <- covariates_data %>%
  mutate(
    # START flag
    start_flag_blow0 = if_else(str_detect(p20031_i0_a0, "START"), 1, 0),
    start_flag_blow1 = if_else(str_detect(p20031_i0_a1, "START"), 1, 0),
    start_flag_blow2 = if_else(str_detect(p20031_i0_a2, "START"), 1, 0),
    )

# Keep only blows with acceptable start
covariates_data$fvc_start_blow0 <- ifelse(covariates_data$start_flag_blow0==0,covariates_data$p3062_i0_a0,NA)
covariates_data$fvc_start_blow1 <- ifelse(covariates_data$start_flag_blow1==0,covariates_data$p3062_i0_a1,NA)
covariates_data$fvc_start_blow2 <- ifelse(covariates_data$start_flag_blow2==0,covariates_data$p3062_i0_a2,NA)

covariates_data$fev1_start_blow0 <- ifelse(covariates_data$start_flag_blow0==0,covariates_data$p3063_i0_a0,NA)
covariates_data$fev1_start_blow1 <- ifelse(covariates_data$start_flag_blow1==0,covariates_data$p3063_i0_a1,NA)
covariates_data$fev1_start_blow2 <- ifelse(covariates_data$start_flag_blow2==0,covariates_data$p3063_i0_a2,NA)

# Count non-missing FVC start blows
covariates_data$count_start <- rowSums(!is.na(covariates_data[, c("fvc_start_blow0", "fvc_start_blow1", "fvc_start_blow2")]))

# Take max and median for three blows; max and min for two blows.
covariates_data <- covariates_data %>%
  mutate(FVC_max = case_when(
    count_start %in% c(2, 3) ~ pmax(fvc_start_blow0,fvc_start_blow1,fvc_start_blow2,na.rm = TRUE),
    TRUE ~ NA_real_))

covariates_data <- covariates_data %>%
  mutate(FEV1_max = case_when(
    count_start %in% c(2, 3) ~ pmax(fev1_start_blow0,fev1_start_blow1,fev1_start_blow2,na.rm = TRUE),
    TRUE ~ NA_real_))

# Take median for three blows and min for two blows.
covariates_data <- covariates_data %>%
  rowwise() %>%
  mutate(FVC_min = case_when(
    count_start == 3 ~ median(c_across(c(fvc_start_blow0, fvc_start_blow1, fvc_start_blow2)), na.rm = TRUE),
    count_start == 2 ~ min(c_across(c(fvc_start_blow0, fvc_start_blow1, fvc_start_blow2)), na.rm = TRUE),
    TRUE ~ NA_real_))

covariates_data <- covariates_data %>%
  rowwise() %>%
  mutate(FEV1_min = case_when(
    count_start ==3 ~ median(c_across(c(fev1_start_blow0,fev1_start_blow1,fev1_start_blow2)),na.rm = TRUE),
    count_start ==2 ~ min(c_across(c(fev1_start_blow0,fev1_start_blow1,fev1_start_blow2)),na.rm = TRUE),
    TRUE ~ NA_real_))

# Take difference between max and min of FVC and FEV1
covariates_data <- covariates_data %>%
  mutate(
    FVC_diff = round((FVC_max-FVC_min)*1000),
    FEV1_diff = round((FEV1_max-FEV1_min)*1000))

# Lable reproducible spirometry measurement
covariates_data$spirogram_reproducible <- ifelse(covariates_data$FVC_diff<=250 & covariates_data$FEV1_diff<=250, 1, 0)
table(covariates_data$spirogram_reproducible, useNA="always")

# Compute FEV1:FVC ratio
covariates_data$FEV1_to_FVC <- covariates_data$FEV1_max/covariates_data$FVC_max

# Prep for Global Lung Initiative (GLI) 2012 equations
covariates_data$sex_numeric <- ifelse(covariates_data$sex == 'Male' , 1, 2)
covariates_data$height_meter <- covariates_data$height / 100
covariates_data <- covariates_data %>%
  mutate(
    ethnicity_GLI = case_when(
      ethnicity %in% c("White", "British", "Irish", "Any other white background") ~ 1,
      ethnicity %in% c("Black or Black British", "Caribbean", "African", "Any other Black background") ~ 2,
      ethnicity == "Chinese" ~ 3,
      ethnicity %in% c("Asian or Asian British", "Indian", "Pakistani", "Bangladeshi", "Any other Asian background") ~ 4,
      ethnicity %in% c("Mixed", "White and Black Caribbean", "White and Black African", 
                       "White and Asian", "Any other mixed background", "Other ethnic group") ~ 5,
      TRUE ~ NA_real_),
    ethnicity_GLI = factor(ethnicity_GLI)
  )

# Calculate the predicted FEV1 value 
covariates_data <- covariates_data %>%
  mutate(
    FEV1_pred = pred_GLI(age, height=height_meter, gender=sex_numeric, ethnicity=ethnicity_GLI, param = "FEV1")
  )

# Calculate FEV1 to predicted FEV1 ratio
covariates_data$FEV1_to_FEV1_pred <- round(covariates_data$FEV1_max / covariates_data$FEV1_pred,2)

# Define impaired lung function based on spirometry measurement
covariates_data <- covariates_data %>%
  mutate(lung_function_impaired = case_when(
    is.na(FEV1_to_FEV1_pred) | is.na(FEV1_to_FVC) ~ NA_character_,
    spirogram_reproducible == 1 & (FEV1_to_FEV1_pred < 0.8 | FEV1_to_FVC < 0.7) ~ "Yes",
    spirogram_reproducible == 1 & (FEV1_to_FEV1_pred >= 0.8 & FEV1_to_FVC >= 0.7) ~"No",
    TRUE ~ NA_character_
  ),
  lung_function_impaired = factor(lung_function_impaired)
  )

# Define impaired lung function from verbal interview, touchscreen questionnaire, and spirometry measurement
covariates_data <- covariates_data %>%
  mutate(lung_disease_history = case_when(
    lung_disease_verbal == "Yes" | lung_disease_touchscreen == "Yes" | lung_function_impaired == "Yes"~ 1,
    lung_disease_verbal == "No" & lung_disease_touchscreen == "No" & lung_function_impaired == "No" ~ 0,
    is.na(lung_disease_verbal) & is.na(lung_disease_touchscreen) & is.na(lung_function_impaired) ~ NA_real_,
    TRUE ~ 0   # mixture of 0 and NA → assume no known history of chronic lung disease
  ),
  lung_disease_history = factor(lung_disease_history, levels = c(0,1)
  ))

######################### Lung cancer data prep ######################### 
ukb_data <- merge(covariates_data, cancer_data,by="eid",sort=FALSE,all.x=FALSE)

# Convert a tibble to data.table
ukb_data <- setDT(ukb_data)

# Convert data to long format
ukbdata_long <- melt(ukb_data,
                     id.vars = c("eid","p53_i0","p54_i0"),
                     measure.vars = patterns("p40006_","p40013_", "p40005_","p40011_","p40012_","p40021_","p40000_","p40001_","p40020_"),
                     value.name=c("p40006","p40013","p40005","p40011","p40012","p40021","p40000","p40001","p40020"))

# Specify ICD codes for lung cancer
lung_cancer_icd9 <- unique(ukbdata_long$p40013[grep("1622|1623|1624|1625|1628|1629",ukbdata_long$p40013)])
lung_cancer_icd10 <- unique(ukbdata_long$p40006[grep("^C34.0|^C34.1|^C34.2|^C34.3|^C34.8|^C34.9",ukbdata_long$p40006)])

# Identify participants with prevalent lung cancer prior to recruitment or up to 2 years after recruitment N=867
prevalent_lungc <- unique(subset(ukbdata_long, (is.element(p40006,lung_cancer_icd10) | is.element(p40013,lung_cancer_icd9)) & (as.Date(p40005) <= (as.Date(p53_i0)+730)))$eid)

# Remove N=867 participants who had prevalent lung cancer at baseline
ukbdata_long <- collect(filter(ukbdata_long,!(is.element(ukbdata_long$eid,prevalent_lungc))))

# Create a subset of data with those who have new diagnosis of lung cancer 
lungcancer_data <- subset(ukbdata_long, is.element(p40006,lung_cancer_icd10) & (p40012 == 3))

# Find the earliest diagnosis date
lungcancer_data$lc_date_diagnosis_min <- as.Date(tapply(lungcancer_data$p40005,lungcancer_data$eid,min)[as.character(lungcancer_data$eid)],origin= "1970-01-01")

# Order by the diagnosis date
lungcancer_data <- lungcancer_data[order(eid,p40005),]

# Combine all unique lung cancer diagnosis dates for each participant
lungcancer_data <- lungcancer_data %>%
  group_by(eid) %>%
  mutate(lungc_date = paste(unique(p40005), collapse = ", ")) %>%
  ungroup()

# Combine all unique lung cancer ICD-10 codes recorded for each participant
lungcancer_data <- lungcancer_data %>%
  group_by(eid) %>%
  mutate(lungc_icd10 = paste(unique(p40006), collapse = ", ")) %>%
  ungroup()

# Retain all ICD-10 records, including repeated codes, for diagnostic checking
lungcancer_data <- lungcancer_data %>%
  group_by(eid) %>%
  mutate(lungc_icd10_2 = paste(p40006, collapse = ", ")) %>%
  ungroup()

#### Derive histological subtype
# Combine all unique histology codes recorded for each participant
lungcancer_data <- lungcancer_data %>%
  group_by(eid) %>%
  mutate(lungc_histology = paste(unique(p40011), collapse = ", ")) %>%
  ungroup()

# Count the number of distinct diagnosis dates and histology codes per participant
lungcancer_data <- lungcancer_data %>%
  group_by(eid) %>%
  mutate(
    n_date = n_distinct(p40005),
    n_hist = n_distinct(p40011)
  ) %>%
  ungroup()

# Histology codes used when resolving participants with multiple histology records
hist_codes <- c(8000,8010,8020,8046)

# Derive one selected histology code per participant.
# For participants with one histology code, retain that code.
# For participants with diagnoses on multiple dates, prioritize the histology associated with the earliest diagnosis.
# For multiple histologies recorded on the same date, apply rules to select the final histology code.
lungcancer_data <- lungcancer_data %>%
  group_by(eid) %>%
  mutate(
    selected_histology = case_when(
      n_hist == 1 ~ first(p40011),  # Only one histology
      n_date > 1 & n_hist > 1 ~ first(p40011),  # Histology with earliest date
      n_date == 1 & n_hist > 1 & all(p40011 %in% c(8000, 8010)) ~ 8010,  
      n_date == 1 & n_hist > 1 & all(p40011 %in% c(8070, 8140)) ~ 8560,  
      n_date == 1 & n_hist > 1 & all(p40011 %in% c(8070, 8560)) ~ 8560,  
      n_date == 1 & n_hist > 1 & all(p40011 %in% c(8010, 8020)) ~ 8010,  
      n_date == 1 & n_hist > 1 & all(p40011 %in% c(8000, 8046)) ~ 8046,  
      n_date == 1 & n_hist > 1 & all(p40011 %in% c(8070, 8045)) ~ 8045,  
      n_date == 1 & n_hist > 1 & all(p40011 %in% c(8230, 8070)) ~ 8070,  
      n_date == 1 & n_hist > 1 & all(p40011 %in% c(8480, 8013)) ~ 8013,  
      n_date == 1 & n_hist > 1 & any(p40011 %in% hist_codes) ~ p40011[!p40011 %in% hist_codes][1],  
      n_date == 1 & n_hist > 1 & any(p40011 %in% 8140) ~ 8140,  
      TRUE ~ NA_integer_
    ))

# Map histology codes to major lung cancer histological subtypes
lungcancer_data <- lungcancer_data %>%
  group_by(eid) %>%
  mutate(
    lungc_histology_final = case_when(
      selected_histology %in% c(8041, 8042, 8045) ~ "Small cell carcinoma",
      selected_histology %in% c(8070, 8071, 8072, 8074, 8123) ~ "Squamous cell carcinoma",
      selected_histology %in% c(8050, 8140, 8230, 8250, 8251, 8252, 8253, 8254, 8255, 8260, 8310, 8323, 8333, 8440, 8480, 8481, 8490, 8550, 8551, 8574) ~ "Adenocarcinoma",
      selected_histology %in% c(8012, 8013, 8034) ~ "Large cell carcinoma",
      selected_histology %in% c(8000, 8010, 8020) ~ "Unspecified  Carcinoma",
      selected_histology %in% c(8046) ~ "Non–small-cell carcinoma, not otherwise specified",
      selected_histology %in% c(8022, 8032, 8033, 8200, 8240, 8244, 8246, 8249, 8430, 8560, 8800, 8802, 8810, 8890, 8972, 8980, 9040, 9120, 9133) ~ "Other Specified  Carcinoma",
      TRUE ~ NA_character_
    ))

# Group histological subtypes into NSCLC, SCLC, and other lung cancers
lungcancer_data <- lungcancer_data %>%
  group_by(eid) %>%
  mutate(
    lungc_histology_group = case_when(
      lungc_histology_final %in% c("Squamous cell carcinoma", "Adenocarcinoma", "Large cell carcinoma", "Non–small-cell carcinoma, not otherwise specified") ~ "NSCLC",
      lungc_histology_final %in% c("Small cell carcinoma") ~ "SCLC",
      lungc_histology_final %in% c("Other Specified  Carcinoma", "Unspecified  Carcinoma") ~ "Others",
      TRUE ~ NA_character_
    ))

# Retain the first record for each participant
lungcancer_data <- ddply(lungcancer_data,.(eid),function(x) head(x,1)) 

# Create a vector containing participant IDs with incident lung cancer
lungcancer_id<- lungcancer_data$eid

# Identify all ICD-10 cancer codes (codes beginning with "C")
any_cancer_icd10 <- unique(ukbdata_long$p40006[grep("^C", ukbdata_long$p40006)])

#### Identify incident primary cancer (non-lung cancer) diagnosed prior to primary lung cancer cases
# Not lung-cancer [f_40006 in any_cancer_icd10]
# Primary cancer [f_40012 == 3]
# Other lung cancer diagnosed prior to lung-cancer [f_40005 < lc_date_diagnosis_min]
othercancer_data <- subset(merge(ukbdata_long,subset(lungcancer_data,select=c(eid,lc_date_diagnosis_min)),by="eid",sort=FALSE,all=FALSE),
                           is.element(p40006,any_cancer_icd10) & (p40012 == 3 ) & (p40005 < lc_date_diagnosis_min))

# Identify the most recent primary cancer diagnosis date before lung cancer
othercancer_data$nlc_date_diagnosis_max <- as.Date(tapply(othercancer_data$p40005,othercancer_data$eid,max)[as.character(othercancer_data$eid)],origin= "1970-01-01")

# Order by the diagnosis date
othercancer_data <- othercancer_data[order(eid,nlc_date_diagnosis_max),]

# Retain the first record for each participant
othercancer_data <- ddply(othercancer_data,.(eid),function(x) head(x,1))

#### Define lung cancer diagnoses and other primary cancers
#b_1: incident lung cancer with no prior cancer diagnosis
#c_2: incident lung cancer with prior cancer (non-lung cancer) diagnosed > 5 years ago
#d_3: incident lung cancer with prior cancer (non-lung cancer) diagnosed <= 5 years ago
lungcancer_data <- merge(lungcancer_data,subset(othercancer_data,select=c(eid,nlc_date_diagnosis_max)),by="eid",sort=FALSE,all.x=TRUE)
lungcancer_data$bin1 <- ifelse(is.na(lungcancer_data$nlc_date_diagnosis_max),1,
                               ifelse(!is.na(lungcancer_data$nlc_date_diagnosis_max) & (lungcancer_data$lc_date_diagnosis_min - lungcancer_data$nlc_date_diagnosis_max) > 1825,2,3))#1825=365*5

b <- lungcancer_data$eid[lungcancer_data$bin1==1]
c <- lungcancer_data$eid[lungcancer_data$bin1==2]
d <- lungcancer_data$eid[lungcancer_data$bin1==3]

# Merge derived lung cancer information back into the longitudinal UKB dataset
ukbdata_long <- merge(ukbdata_long,subset(lungcancer_data,select=c(eid,lc_date_diagnosis_min,nlc_date_diagnosis_max,bin1,lungc_icd10,lungc_histology,selected_histology,lungc_histology_final,lungc_histology_group)),
                      by="eid",sort=FALSE,all.x=TRUE)

# Identify those with other incident cancer but no lung cancer
incident_oc <- subset(ukbdata_long, (is.element(p40006,any_cancer_icd10)) & (p40012 == 3) & (is.na(lc_date_diagnosis_min)))

#### Finalize cancer-status categories
#0_a cancer-free
#1_b incident lung cancer with no prior cancer diagnosis
#2_c incident lung cancer with prior cancer (non-lung cancer) diagnosed > 5 years ago
#3_d incident lung cancer with prior cancer (non-lung cancer) diagnosed <= 5 years ago
#4_e incident other cancer but no lung cancer
ukbdata_long$bin1 <- ifelse(is.na(ukbdata_long$bin1),0,ukbdata_long$bin1)
ukbdata_long$bin1 <- ifelse(is.element(ukbdata_long$eid,incident_oc$eid),4,ukbdata_long$bin1)
ukbdata_long$bin1 <-as.factor(ukbdata_long$bin1)
a<-unique(ukbdata_long$eid[ukbdata_long$bin1==0])
e<-unique(ukbdata_long$eid[ukbdata_long$bin1==4])

# Retain original UKB records for those without lung cancer at baseline
ukb_data2 <- collect(filter(ukb_data,is.element(ukb_data$eid,ukbdata_long$eid)))

# Match each participant in ukb_data2 to the corresponding row in ukbdata_long
idx <- match(ukb_data2$eid, ukbdata_long$eid)  

# Add the earliest lung cancer diagnosis date from ukbdata_long to ukb_data2
ukb_data2$lungc_date <- ukbdata_long$lc_date_diagnosis_min[idx]
ukb_data2$lungc_date <- as.Date(ukb_data2$lungc_date, origin= "1970-01-01")

# Add lung cancer histology classifications from ukbdata_long to ukb_data2
ukb_data2$lungc_histology_final <- ukbdata_long$lungc_histology_final[idx]
ukb_data2$lungc_histology_group <- ukbdata_long$lungc_histology_group[idx]

# Add lung cancer ICD10 code from ukbdata_long to ukb_data2
ukb_data2$lungc_icd10 <- ukbdata_long$lungc_icd10[idx]

# Group UK Biobank assessment centres by country (England, Wales, or Scotland) 
ukb_data2$centre <- ifelse(is.element(ukb_data2$p54_i0,c("Manchester","Oxford","Stoke", "Reading","Bury", "Newcastle", "Leeds", "Bristol", "Barts", "Nottingham", "Sheffield", "Liverpool", "Middlesbrough", "Hounslow", "Croydon", "Birmingham","Stockport (pilot)")),"England",
                           ifelse(is.element(ukb_data2$p54_i0,c("Cardiff","Swansea","Wrexham")),"Wales",
                                  ifelse(is.element(ukb_data2$p54_i0,c("Glasgow","Edinburgh")),"Scotland","NA")))

# Ensure all date objects are saved as Date
ukb_data2$lungc_date <- as.Date(ukb_data2$lungc_date)
ukb_data2$p40000_i0 <- as.Date(ukb_data2$p40000_i0)
ukb_data2$p53_i0 <- as.Date(ukb_data2$p53_i0)

# Calculate follow-up time on study
ukb_data2$lungc_time <- ifelse(!is.na(ukb_data2$lungc_date),(ukb_data2$lungc_date-ukb_data2$p53_i0)/365, 
                               ifelse(((is.na(ukb_data2$lungc_date)) & (!is.na(ukb_data2$p40000_i0))), (ukb_data2$p40000_i0-ukb_data2$p53_i0)/365,
                                      ifelse(((is.na(ukb_data2$lungc_date)) & (is.na(ukb_data2$p40000_i0)) & ukb_data2$centre=="England"),as.numeric(difftime("2023-5-31", ukb_data2$p53_i0, units="days"))/365,
                                             ifelse(((is.na(ukb_data2$lungc_date)) & (is.na(ukb_data2$p40000_i0)) & ukb_data2$centre=="Wales"),as.numeric(difftime("2016-12-31", ukb_data2$p53_i0, units="days"))/365,
                                                    ifelse(((is.na(ukb_data2$lungc_date)) & (is.na(ukb_data2$p40000_i0)) & ukb_data2$centre=="Scotland"),as.numeric(difftime("2023-9-30", ukb_data2$p53_i0, units="days"))/365,NA))))) 

# Create a default value of NA for lungc_status
ukb_data2$lungc_status <- NA
# Assign lungc_status based on membership in each of the groups (a, b, c, d, e)
#0 cancer-free
#1 lung cancer only 
#2 other incident cancer prior to lung cancer greater than 5 years
#3 other incident cancer prior to lung cancer within 5 years
#4 other incident cancer only
ukb_data2$lungc_status[ukb_data2$eid %in% a] <- 0
ukb_data2$lungc_status[ukb_data2$eid %in% b] <- 1
ukb_data2$lungc_status[ukb_data2$eid %in% c] <- 2
ukb_data2$lungc_status[ukb_data2$eid %in% d] <- 3
ukb_data2$lungc_status[ukb_data2$eid %in% e] <- 4

# Exclude those with other incident cancer prior to lung cancer within 5 years of diagnosis
ukb_data2 <-
  ukb_data2 %>%
  filter(!(lungc_status==3)) %>%
  collect()

######################### Metabolomics data prep ######################### 
# Change column names
metabolites <- c("Total_C", "non_HDL_C", "Remnant_C", "VLDL_C", "Clinical_LDL_C", "LDL_C", "HDL_C", 
                 "Total_TG", "VLDL_TG", "LDL_TG", "HDL_TG", 
                 "Total_PL", "VLDL_PL", "LDL_PL", "HDL_PL", 
                 "Total_CE", "VLDL_CE", "LDL_CE", "HDL_CE", 
                 "Total_FC", "VLDL_FC", "LDL_FC", "HDL_FC", 
                 "Total_L", "VLDL_L", "LDL_L", "HDL_L", 
                 "Total_P", "VLDL_P", "LDL_P", "HDL_P", 
                 "VLDL_PS", "LDL_PS", "HDL_PS", 
                 "Phosphoglycerides", 
                 "TG/PG", 
                 "Total_Cholines", "Phosphatidylcholines", "Sphingomyelins", 
                 "ApoB", "ApoA1", "ApoB/ApoA1", 
                 "Total_FA", "Unsaturation", "Omega_3", "Omega_6", "PUFA", "MUFA", "SFA", "LA", "DHA", 
                 "Omega_3_%", "Omega_6_%", "PUFA_%", "MUFA_%", "SFA_%", "LA_%", "DHA_%", "PUFA/MUFA", "Omega_6/Omega_3", 
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
                 "S_HDL_PL_%", "S_HDL_C_%", "S_HDL_CE_%", "S_HDL_FC_%", "S_HDL_TG_%",
                 "Glucose_lactate", "Spectrometer_corrected_alanine")

colnames(metabolite_data) [-1] <- metabolites

# Create a function for single-metabolite preprocessing
run_metabolites_preprocessing <- function(metabolites, metabolite_data) {
  
  # Keep eid + the metabolite
  metabolite_df <- subset(metabolite_data, select = c("eid", metabolites))
  
  # Remove missing values in the target metabolite
  metabolite_df <- metabolite_df %>%
    filter(!is.na(.data[[metabolites]]))
  
  # Log-transform
  metabolite_df[[metabolites]] <- log1p(metabolite_df[[metabolites]])
  
  # Remove outliers (Z-score method)
  metabolite_df <- metabolite_df %>%
    filter(abs(!!rlang::sym(metabolites) - mean(!!rlang::sym(metabolites), na.rm = TRUE)) /
             sd(!!rlang::sym(metabolites), na.rm = TRUE) < 4)
  
  # Scale (coerce matrix -> numeric)
  metabolite_df[[metabolites]] <- as.numeric(scale(metabolite_df[[metabolites]]))
  
  return(metabolite_df)
}

# Create a function for all metabolites and combine by eid
run_full_metabolites_preprocessing <- function(metabolites, metabolite_data) {
  
  # Process each metabolite and store results
  res_list <- vector("list", length(metabolites))
  names(res_list) <- metabolites
  for (i in seq_along(metabolites)) {
    res_list[[i]] <- run_metabolites_preprocessing(metabolites[i], metabolite_data)
  }
  
  # Merge all on eid (row counts may differ after filtering)
  combined_metabolites_preprocessed <- Reduce(function(x, y) {
    merge(x, y, by = "eid", all = TRUE)
  }, res_list)
  
  return(combined_metabolites_preprocessed)
}

# Run across all metabolites 
metabolites_preprocessed_data <- run_full_metabolites_preprocessing(metabolites, metabolite_data)

######################### CoxPH data prep #########################
# Merge covariates data, lung cancer outcome data, and metabolites data
ukb_data3 <- merge(ukb_data2, metabolites_preprocessed_data, by="eid",sort=FALSE,all.x=FALSE)

# Update lung cancer status for CoxPH analysis
#2 other incident cancer prior to lung cancer greater than 5 years - Make lung cancer status as 1
#4 other incident cancer only Make lung cancer status as 0
coxph_data <- ukb_data3
coxph_data$lungc_status_coxph <- ifelse(coxph_data$lungc_status == 2,1,
                                        ifelse(coxph_data$lungc_status == 4,0,coxph_data$lungc_status))

# Make appropriate lung cancer status for cause-specific analysis 
coxph_data <- coxph_data %>%
  mutate(
    event_NSCLC = ifelse(lungc_status_coxph == 1 & lungc_histology_group == "NSCLC", 1, 0),
    event_SCLC  = ifelse(lungc_status_coxph == 1 & lungc_histology_group == "SCLC", 1, 0),
    event_Other = ifelse(lungc_status_coxph == 1 & lungc_histology_group == "Others", 1, 0)
  )

# Define covariates included in the fully adjusted CoxPH models
extended_covariates_original <- c("age", "sex", "ethnicity_UKB", "education_ISCED","tdi","bmi","bmi2","waist_to_hip","smk_status","pack_year","family_lc","lung_disease_history","lipid_lowering_meds","Summed_MET_minutes_weekly","fasting_time","centre_group")

# Create the complete-case analysis dataset by excluding participants with missing values in any covariates
coxph_data_complete_case <- coxph_data %>%
  drop_na(all_of(extended_covariates_original))

# Save data
setwd("M:/shared/Beiwen Wu/Project 1/UKB_RAP/UKB_Data/Processed")
write.csv(metabolites_preprocessed_data, "metabolite_data_preprocessed.csv", row.names = FALSE)
saveRDS(covariates_data,"covariates_data_501927.rds")
saveRDS(ukb_data, "ukb_data_501927.rds")
saveRDS(ukb_data2, "ukb_data2_500454.rds")
saveRDS(ukb_data3, "ukb_data3_486658.rds")
saveRDS(coxph_data, "coxph_data_486658.rds")
saveRDS(coxph_data_complete_case, "coxph_data_complete_301425.rds")


