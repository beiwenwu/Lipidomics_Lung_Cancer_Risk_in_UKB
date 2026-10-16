library(dplyr)
library(tidyr)
library(stringr)
library(forestplot)
library(ggplot2)
library(mice)
library(writexl)
library(patchwork)

################ Run PCA-Cox with Relative Metabolites ################
# Define one common PCA sample - Restrict every imputation to those same participants ########## 
# Remove participants with missing values in any relative metabolite
data_ratio <- lapply(
  imputed_data,
  function(dat) {
    dat %>%
      drop_na(all_of(metabolite_ratio))
  }
)

# Retrieve the number of participants in each imputed dataset
sapply(data_ratio, nrow)

# Extract the participant IDs from each metabolite-complete imputed dataset
common_ids_ratio_each <- lapply(
  data_ratio,
  function(dat) {
    dat[["eid"]]
  }
)

# Identify the common IDs across all ten metabolite-complete imputed datasets
common_ids_ratio <- Reduce(
  intersect,
  common_ids_ratio_each
)

# Keep only the common participants across all ten metabolite-complete imputed datasets
data_ratio_complete <- lapply(
  data_ratio,
  function(dat) {
    idx <- match(common_ids_ratio, dat[["eid"]])
    dat[idx, , drop = FALSE]
  }
)

# Extract the relative metabolites measures for those participants with common IDs
metabolite_ratio_data <- data_ratio_complete$imp1 %>%
  select(eid,all_of(metabolite_ratio))

# Perform principal component analysis on the relative metabolite measures
pca_result_ratio <- prcomp(metabolite_ratio_data[,-1], scale.= FALSE)

# Display a summary of the PCA results, including standard deviation, proportion of variance explained, cumulative proportion of variance explained
summary(pca_result_ratio)

# Extracts the individual participant PC scores
pca_scores_ratio <- data.frame(eid = metabolite_ratio_data$eid, pca_result_ratio$x)

# Calculate the percentage of total variance explained by each principal component
explained_variance_ratio <- pca_result_ratio$sdev^2 / sum(pca_result_ratio$sdev^2) * 100

# Extract loading: How much each metabolite contributes to each principal component
loadings_ratio <- pca_result_ratio$rotation

# Take the absolute values of the loading: magnitude of the loading regardless of the direction 
abs_loadings_ratio <- abs(loadings_ratio)

# Sort metabolites by the absolute magnitude of their loading
pc1_ratio_contributions <- sort(abs_loadings_ratio[, 1], decreasing = TRUE)
pc2_ratio_contributions <- sort(abs_loadings_ratio[, 2], decreasing = TRUE)
pc3_ratio_contributions <- sort(abs_loadings_ratio[, 3], decreasing = TRUE)
pc4_ratio_contributions <- sort(abs_loadings_ratio[, 4], decreasing = TRUE)
pc5_ratio_contributions <- sort(abs_loadings_ratio[, 5], decreasing = TRUE)
pc6_ratio_contributions <- sort(abs_loadings_ratio[, 6], decreasing = TRUE)
pc7_ratio_contributions <- sort(abs_loadings_ratio[, 7], decreasing = TRUE)
pc8_ratio_contributions <- sort(abs_loadings_ratio[, 8], decreasing = TRUE)
pc9_ratio_contributions <- sort(abs_loadings_ratio[, 9], decreasing = TRUE)

# Get the top 10 metabolites contributing to PCs:
top_pc1_ratio <- head(pc1_ratio_contributions,10)
top_pc2_ratio <- head(pc2_ratio_contributions,10)
top_pc3_ratio <- head(pc3_ratio_contributions,10)
top_pc4_ratio <- head(pc4_ratio_contributions,10)
top_pc5_ratio <- head(pc5_ratio_contributions,10)
top_pc6_ratio <- head(pc6_ratio_contributions,10)
top_pc7_ratio <- head(pc7_ratio_contributions,10)
top_pc8_ratio <- head(pc8_ratio_contributions,10)
top_pc9_ratio <- head(pc9_ratio_contributions,10)

# Extract the metabolite loading matrix for the first 9 principal components
L_ratio <- pca_result_ratio$rotation[, paste0("PC", 1:9), drop = FALSE]

# Identify the top-k metabolites with the largest absolute loadings for each principal component and store them in a tidy table
k <- 10
top_loads_ratio <- lapply(seq_len(ncol(L_ratio)), function(x) {
  pc <- colnames(L_ratio)[x] # get the name of the current PC
  v <- L_ratio[, x] # extract all metabolite loadings for the current PC
  ord <- order(abs(v), decreasing = TRUE) # take the absolute value of each loading and rank from largest to smallest
  sel <- ord[seq_len(k)] # select the first k metabolites from the ranking
  labs <- names(v)[sel] # extract the names of the selected metabolites
  tibble(PC = pc,
         top_metabolites = paste(labs, collapse = ", ")) # create one row for the current PC, combining all 10 metabolite names separated by commas
}) %>% bind_rows() # combine the 9 individual tibbles into a single data frame

# Merge PCA scores into each metabolite-complete imputed dataset and retain the variables required for the PCA-Cox analyses
pcs <- paste0("PC", 1:9)
pca_data_ratio <- lapply(
  data_ratio_complete,
  function(dat) {
    dat %>%
      left_join(
        pca_scores_ratio,
        by = "eid"
      ) %>%
      select(
        eid,
        all_of(pcs),
        lungc_time,
        lungc_status_coxph,
        all_of(extended_covariates)
      ) %>%
      arrange(eid)
  }
)

################ Model 1: Single PC for relative metabolites ################
fit_one_pc_mi <- function(pc, data_list, covariates_use) {
  
  # Construct the Cox model formula
  formula_adjusted <- as.formula(
    paste0(
      "Surv(lungc_time, lungc_status_coxph) ~ ",
      paste0("`", pc, "`"),
      " + ",
      paste0("`", covariates_use, "`", collapse = " + ")
    )
  )
  
  # Fit the same model in every imputed dataset
  fits <- lapply(
    data_list,
    function(dat) {
      coxph(
        formula = formula_adjusted,
        data = dat,
        ties = "efron",
        x = TRUE
      )
    }
  )
  
  # Convert the model list to a mira object and pool with Rubin's rules
  pooled_fit <- pool(
    as.mira(fits)
  )
  
  # Extract pooled HR, CI, and p-value for the PC
  pooled_result <- summary(
    pooled_fit,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    filter(term == pc) %>%
    transmute(
      PC  = pc,
      HR  = estimate,
      LCL = `2.5 %`,
      UCL = `97.5 %`,
      p   = p.value
    )
  
  pooled_result
}

# Run PCA-Cox analysis by fitting one PC at a time
pca_ratio_1pc <- bind_rows(
  lapply(
    pcs,
    fit_one_pc_mi,
    data_list = pca_data_ratio,
    covariates_use = extended_covariates
  )
) %>%
  mutate(
    FDR = p.adjust(p, method = "BH")
  ) %>%
  arrange(PC)

# Create a table for PCA-Cox results along with top contributing metabolites
tbl_ratio_1pc <- pca_ratio_1pc %>%
  mutate(
    HR_CI = sprintf("%.3f (%.3f–%.3f)", HR, LCL, UCL),
    p_adj = ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR))
  ) %>%
  left_join(top_loads_ratio, by = "PC") %>%
  arrange(PC)

################ Model 2: Sequential PCs for relative metabolites ################ 
fit_sequential_pcs_mi <- function(pcs, data_list, covariates_use) {
  
  # Loop through all the PCs, run a function for each one
  map_dfr(seq_along(pcs), function(i) {
    
    # PC newly evaluated in this model
    current_pc <- pcs[i]
    
    # Current PC plus all preceding PCs
    pc_terms <- pcs[seq_len(i)]
    
    # Construct the right-hand side of CoxPH formula
    rhs <- paste(
      c(pc_terms, covariates_use), collapse = " + "
    )
    
    # Construct the full CoxPH formula
    formula_adjusted <- as.formula(
      paste0(
        "Surv(lungc_time, lungc_status_coxph) ~ ", rhs
      )
    )
    
    # Fit the same sequential model in every imputed dataset
    fits <- lapply(
      data_list,
      function(dat) {
        coxph(
          formula = formula_adjusted,
          data = dat,
          ties = "efron"
        )
      }
    )
    
    # Convert the model list to a mira object and pool with Rubin's rules
    pooled_fit <- pool(
      as.mira(fits)
    )
    
    # Extract pooled HR, CI, and p-value only for the PC only the newly in the current loop
    summary(
      pooled_fit,
      conf.int = TRUE,
      exponentiate = TRUE
    ) %>%
      filter(term == current_pc) %>%
      transmute(
        model_number = i,
        PC = term,
        PCs_in_model = paste(pc_terms, collapse = " + "),
        HR = estimate,
        LCL = `2.5 %`,
        UCL = `97.5 %`,
        p = p.value
      )
    
  })
}

# Run PCA-Cox analysis by fitting the current PC plus all preceding PCs
pca_ratio_seq <- fit_sequential_pcs_mi(
  pcs = pcs,
  data_list = pca_data_ratio,
  covariates_use = extended_covariates
) %>%
  mutate(
    FDR = p.adjust(p, method = "BH")
  ) %>%
  arrange(PC)

# Create a table for PCA-Cox results
tbl_ratio_seq <- pca_ratio_seq %>%
  mutate(
    HR_CI = sprintf("%.3f (%.3f–%.3f)", HR, LCL, UCL),
    p_adj = ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR))
  )

################ Model 3: All PCs for relative metabolites ################ 
fit_all_pcs_mi <- function(data_list, pcs, covariates_use) {
  
  # Construct the CoxPH formula
  formula_adjusted <- as.formula(
    paste0(
      "Surv(lungc_time, lungc_status_coxph) ~ ",
      paste0("`", pcs, "`", collapse = " + "),
      " + ",
      paste0("`", covariates_use, "`", collapse = " + ")
    )
  )
  
  # Fit the same model in every imputed dataset
  fits <- lapply(
    data_list,
    function(dat) {
      coxph(
        formula = formula_adjusted,
        data = dat,
        ties = "efron"
      )
    }
  )
  
  # Convert the model list to a mira object and pool with Rubin's rules
  pooled_fit <- pool(
    as.mira(fits)
  )
  
  # Extract pooled HR, CI, and p-value for all the PCs
  summary(
    pooled_fit,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    filter(term %in% pcs) %>%
    transmute(
      PC = term,
      HR = estimate,
      LCL = `2.5 %`,
      UCL = `97.5 %`,
      p = p.value
    )
}

# Run PCA-Cox analysis by fitting all the PCs in the same model
pca_ratio_allpc <- fit_all_pcs_mi(
  data_list = pca_data_ratio,
  pcs = pcs,
  covariates_use = extended_covariates
) %>%
  mutate(
    FDR = p.adjust(p, method = "BH")
  )

# Create a table for PCA-Cox results
tbl_ratio_allpcs <- pca_ratio_allpc %>%
  mutate(
    HR_CI = sprintf("%.3f (%.3f–%.3f)", HR, LCL, UCL),
    p_adj = ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR))
  )

################ Save PCA-Cox Results for Relative Metabolites ################
write_xlsx(tbl_ratio_1pc, "pca_ratio_1pc.xlsx")
write_xlsx(tbl_ratio_seq, "pca_ratio_seq.xlsx")
write_xlsx(tbl_ratio_allpcs, "pca_ratio_allpcs.xlsx")