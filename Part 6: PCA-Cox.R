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
setwd('/Users/beiwen/Library/CloudStorage/OneDrive-UniversityofToronto/PhD Documents/PhD Projects/Objective 1/UKB RAP/Results/PCA')
write_xlsx(tbl_ratio_1pc, "pca_ratio_1pc.xlsx")
write_xlsx(tbl_ratio_seq, "pca_ratio_seq.xlsx")
write_xlsx(tbl_ratio_allpcs, "pca_ratio_allpcs.xlsx")

################ Plotting: Scree Plot ################ 
tiff("screeplot_pca_ratio.tiff",
     width = 280, height = 160, units = "mm", res = 600)

plot(
  cumsum(explained_variance_ratio[1:15]),
  type = "b",
  pch = 19,
  cex = 0.4,
  xlab = "Principal Component",
  ylab = "Cumulative Explained Variance",
  main = "Scree plot of relative metabolites in principal component analysis"
)

axis(1, at = 1:15, labels = c(1:15))

# Add 90% line
abline(h = 90, col = "red", lty = 2, lwd = 2)

dev.off()

################ Plotting: Forest Plot ################ 
# Combine results from the three PCA modeling strategies into one dataframe
pca_table_ratio <- bind_rows(
  # Model 1: Each PC fitted separately
  tbl_ratio_1pc %>%
    mutate(model = "Model 1"),
  
  # Model 2: Sequentially PC fitted
  tbl_ratio_seq %>%
    mutate(model = "Model 2"),
  
  # Model 3: All PCs fitted simultaneously
  tbl_ratio_allpcs %>%
    mutate(model = "Model 3")
  
) %>%
  
  # Retain variables required for the final results table
  select(
    PC, model, HR, LCL, UCL, FDR, p_adj
  ) %>%
  
  # Arrange results by principal component and model
  arrange(model, PC)

pca_table_ratio$PC <- factor(pca_table_ratio$PC,
                       levels = paste0("PC", 9:1))  # reverses order

pca_table_ratio$model <- factor(pca_table_ratio$model,
                          levels = paste0("Model ", 3:1))  # reverses order


pca_ratio_forest <- ggplot(pca_table_ratio, aes(x = HR, y = PC, xmin = LCL, xmax = UCL, shape = model)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbarh(
    aes(group = model),
    position = position_dodge(width = 0.6),
    height = 0.2
  ) +
  geom_point(
    aes(group = model),
    position = position_dodge(width = 0.6),
    size = 2.5
  )  +
  scale_shape_manual(
    values = c(16, 17, 15),
    breaks = c("Model 1", "Model 2", "Model 3")
  ) +
  labs(
    x = "Hazard ratio (95% CI)",
    y = "Principle component (PCs)",
    shape = "Model",
    title = "Associations of PCs with lung cancer risk"
  ) +  
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
  )

pca_ratio_forest

################ Plotting PCA loadings ################ 
loadings_ratio_df <- as.data.frame(pca_result_ratio$rotation)

# Change the rowname for display purpose
rownames(loadings_ratio_df)[rownames(loadings_ratio_df) == "Omega_3_%"] <- "Omega_3_FA_%"
rownames(loadings_ratio_df)[rownames(loadings_ratio_df) == "Omega_6_%"] <- "Omega_6_FA_%"
rownames(loadings_ratio_df)[rownames(loadings_ratio_df) == "Omega_6/Omega_3"] <- "Omega_6_FA/Omega_3_FA"
loadings_ratio_df$metabolite <- rownames(loadings_ratio_df)

# Convert to long format
loadings_ratio_df_long <- loadings_ratio_df %>%
  dplyr::select(metabolite, all_of(pcs)) %>%
  tidyr::pivot_longer(
    cols = all_of(pcs),
    names_to = "PC",
    values_to = "loading"
  ) 

# Set the order of metabolite group for display, 
loadings_ratio_df_long <- loadings_ratio_df_long %>%
  mutate(
    # identify lipid component
    component = case_when(
      str_detect(metabolite, "Apo|TG/PG|Omega_6_FA/Omega_3_FA|PUFA/MUFA") ~ "Lipid ratios",
      str_detect(metabolite, "Omega|DHA|PUFA|MUFA|SFA|LA") ~ "Fatty acids",
      str_detect(metabolite, "_C_%") ~ "C%",
      str_detect(metabolite, "_FC_%") ~ "FC%",
      str_detect(metabolite, "_CE_%") ~ "CE%",
      str_detect(metabolite, "_TG_%") ~ "TG%",
      str_detect(metabolite, "_PL_%") ~ "PL%"
    ),
    component = factor(component, levels = c("C%", "FC%", "CE%", "TG%", "PL%","Fatty acids","Lipid ratios")
    ),
    # Identify the direction of loadings
    direction = if_else(loading >= 0, "Positive", "Negative"),
    abs_loading = abs(loading)
  ) %>%
  arrange(component)

# Set the metabolite order for plotting
met_order <- loadings_ratio_df_long %>%
  distinct(component, metabolite) %>%
  arrange(component) %>%
  pull(metabolite)

# Reverse met_order, so the desired first metabolite appears at the top of the y-axis
loadings_ratio_df_long <- loadings_ratio_df_long %>%
  mutate(
    metabolite = factor(metabolite, levels = rev(met_order))
  )

# Plot the PCA loadings
factor_loadings_ratio <-
  ggplot(loadings_ratio_df_long, aes(x = PC, y = metabolite)) +
  geom_point(aes(size = abs_loading, color = direction), alpha = 0.8) +
  scale_size(name = "Loadings", range = c(0.1, 5)) +
  scale_color_manual(
    values = c("Positive" = "red3", "Negative" = "darkgreen"),
    name = "Direction"
  ) +
  labs(
    x = "Principal components (PCs)",
    y = "Metabolites",
    title = "The factor loadings of principal components"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(
      hjust = 0.5
    )
  )

# Combine the forest plot and the loading plot side-by-side
p <- pca_ratio_forest | factor_loadings_ratio

tiff("forest_factor_loadings_ratio.tiff",
     width = 400, height = 200, units = "mm", res = 900)
p

dev.off()

################ Run PCA in Absolute Metabolites ################
pca_data_absolute <- readRDS('/Users/beiwen/Library/CloudStorage/OneDrive-UniversityofToronto/PhD Documents/PhD Projects/Objective 1/UKB RAP/Results/V5-change hypertension_bmi2_ethnicity/without_hypertension/m10_maxit20/pca/pca_data_abs.rds')

# Remove participants with missing values in any absolute metabolite
data_absolute <- lapply(
  imputed_data,
  function(dat) {
    dat %>%
      drop_na(all_of(metabolite_absolute))
  }
)

# Retrieve the number of participants in each imputed dataset
sapply(data_absolute, nrow)

# Extract the participant IDs from each metabolite-complete imputed dataset
common_ids_each <- lapply(
  data_absolute,
  function(dat) {
    dat[["eid"]]
  }
)

# Identify the common IDs across all ten metabolite-complete imputed datasets
common_ids <- Reduce(
  intersect,
  common_ids_each
)

# Keep only the common participants across all ten metabolite-complete imputed datasets
data_absolute_complete <- lapply(
  data_absolute,
  function(dat) {
    idx <- match(common_ids, dat[["eid"]])
    dat[idx, , drop = FALSE]
  }
)

# Extract the absolute metabolites measures for those participants with common IDs
metabolite_absolute_data <- data_absolute_complete$imp1 %>%
  select(eid,all_of(metabolite_absolute))

# Perform principal component analysis on the absolute metabolite measures
pca_result_absolute <- prcomp(metabolite_absolute_data[,-1], scale.= FALSE)

# Display a summary of the PCA results, including standard deviation, proportion of variance explained, cumulative proportion of variance explained
summary(pca_result_absolute)

# Extracts the individual participant PC scores
pca_scores_absolute <- data.frame(eid = metabolite_absolute_data$eid, pca_result_absolute$x)

# Calculate the percentage of total variance explained by each principal component
explained_variance_absolute <- pca_result_absolute$sdev^2 / sum(pca_result_absolute$sdev^2) * 100

# Extract loading: How much each metabolite contributes to each principal component
loadings_absolute <- pca_result_absolute$rotation

# Take the absolute values of the loading: magnitude of the loading regardless of the direction 
abs_loadings_absolute <- abs(loadings_absolute)

# Sort metabolites by the absolute magnitude of their loading
pc1_absolute_contributions <- sort(abs_loadings_absolute[, 1], decreasing = TRUE)
pc2_absolute_contributions <- sort(abs_loadings_absolute[, 2], decreasing = TRUE)
pc3_absolute_contributions <- sort(abs_loadings_absolute[, 3], decreasing = TRUE)
pc4_absolute_contributions <- sort(abs_loadings_absolute[, 4], decreasing = TRUE)
pc5_absolute_contributions <- sort(abs_loadings_absolute[, 5], decreasing = TRUE)
pc6_absolute_contributions <- sort(abs_loadings_absolute[, 6], decreasing = TRUE)
pc7_absolute_contributions <- sort(abs_loadings_absolute[, 7], decreasing = TRUE)
pc8_absolute_contributions <- sort(abs_loadings_absolute[, 8], decreasing = TRUE)
pc9_absolute_contributions <- sort(abs_loadings_absolute[, 9], decreasing = TRUE)

# Get the top 10 metabolites contributing to PCs:
top_pc1_absolute <- head(pc1_absolute_contributions,10)
top_pc2_absolute <- head(pc2_absolute_contributions,10)
top_pc3_absolute <- head(pc3_absolute_contributions,10)
top_pc4_absolute <- head(pc4_absolute_contributions,10)
top_pc5_absolute <- head(pc5_absolute_contributions,10)
top_pc6_absolute <- head(pc6_absolute_contributions,10)
top_pc7_absolute <- head(pc7_absolute_contributions,10)
top_pc8_absolute <- head(pc8_absolute_contributions,10)
top_pc9_absolute <- head(pc9_absolute_contributions,10)

# Extract the metabolite loading matrix for the first 9 principal components
L_absolute <- pca_result_absolute$rotation[, paste0("PC", 1:9), drop = FALSE]

# Identify the top-k metabolites with the largest absolute loadings for each principal component and store them in a tidy table
k <- 10
top_loads_absolute <- lapply(seq_len(ncol(L_absolute)), function(x) {
  pc <- colnames(L_absolute)[x] # get the name of the current PC
  v <- L_absolute[, x] # extract all metabolite loadings for the current PC
  ord <- order(abs(v), decreasing = TRUE) # take the absolute value of each loading and rank from largest to smallest
  sel <- ord[seq_len(k)] # select the first k metabolites from the ranking
  labs <- names(v)[sel] # extract the names of the selected metabolites
  tibble(PC = pc,
         top_metabolites = paste(labs, collapse = ", ")) # create one row for the current PC, combining all 10 metabolite names separated by commas
}) %>% bind_rows() # combine the 9 individual tibbles into a single data frame

# Merge PCA scores into each metabolite-complete imputed dataset and retain the variables required for the PCA-Cox analyses
pcs <- paste0("PC", 1:9)
pca_data_absolute <- lapply(
  data_absolute_complete,
  function(dat) {
    dat %>%
      left_join(
        pca_scores_absolute,
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

################ Model 1: Single PC for absolute metabolites ################
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
pca_absolute_1pc <- bind_rows(
  lapply(
    pcs,
    fit_one_pc_mi,
    data_list = pca_data_absolute,
    covariates_use = extended_covariates
  )
) %>%
  mutate(
    FDR = p.adjust(p, method = "BH")
  ) %>%
  arrange(PC)

# Create a table for PCA-Cox results along with top contributing metabolites
tbl_absolute_1pc <- pca_absolute_1pc %>%
  mutate(
    HR_CI = sprintf("%.3f (%.3f–%.3f)", HR, LCL, UCL),
    p_adj = ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR))
  ) %>%
  left_join(top_loads_absolute, by = "PC") %>%
  arrange(PC)

################ Model 2: Sequential PCs for absolute metabolites ################ 
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
pca_absolute_seq <- fit_sequential_pcs_mi(
  pcs = pcs,
  data_list = pca_data_absolute,
  covariates_use = extended_covariates
) %>%
  mutate(
    FDR = p.adjust(p, method = "BH")
  ) %>%
  arrange(PC)

# Create a table for PCA-Cox results
tbl_absolute_seq <- pca_absolute_seq %>%
  mutate(
    HR_CI = sprintf("%.3f (%.3f–%.3f)", HR, LCL, UCL),
    p_adj = ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR))
  )

################ Model 3: All PCs for absolute metabolites ################ 
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
pca_absolute_allpc <- fit_all_pcs_mi(
  data_list = pca_data_absolute,
  pcs = pcs,
  covariates_use = extended_covariates
) %>%
  mutate(
    FDR = p.adjust(p, method = "BH")
  )

# Create a table for PCA-Cox results
tbl_absolute_allpcs <- pca_absolute_allpc %>%
  mutate(
    HR_CI = sprintf("%.3f (%.3f–%.3f)", HR, LCL, UCL),
    p_adj = ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR))
  )

################ Save PCA-Cox Results for absolute Metabolites ################
setwd('/Users/beiwen/Library/CloudStorage/OneDrive-UniversityofToronto/PhD Documents/PhD Projects/Objective 1/UKB RAP/Results/PCA')
write_xlsx(tbl_absolute_1pc, "pca_absolute_1pc.xlsx")
write_xlsx(tbl_absolute_seq, "pca_absolute_seq.xlsx")
write_xlsx(tbl_absolute_allpcs, "pca_absolute_allpcs.xlsx")

################ Plotting: Scree Plot ################ 
tiff("screeplot_pca_abs.tiff",
     width = 280, height = 160, units = "mm", res = 600)

plot(
  cumsum(explained_variance_absolute[1:15]),
  type = "b",
  pch = 19,
  cex = 0.4,
  xlab = "Principal Component",
  ylab = "Cumulative Explained Variance",
  main = "Scree plot of absolute metabolites in principal component analysis"
)
axis(1, at = 1:15, labels = c(1:15))

# Add 90% line
abline(h = 90, col = "red", lty = 2, lwd = 2)

dev.off()

################ Plotting: Forest Plot ################ 
# Combine results from the three PCA modeling strategies into one dataframe
pca_table_absolute <- bind_rows(
  # Model 1: Each PC fitted separately
  tbl_absolute_1pc %>%
    mutate(model = "Model 1"),
  
  # Model 2: Sequentially PC fitted
  tbl_absolute_seq %>%
    mutate(model = "Model 2"),
  
  # Model 3: All PCs fitted simultaneously
  tbl_absolute_allpcs %>%
    mutate(model = "Model 3")
  
) %>%
  
  # Retain variables required for the final results table
  select(
    PC, model, HR, LCL, UCL, FDR, p_adj
  ) %>%
  
  # Arrange results by principal component and model
  arrange(model, PC)

pca_table_absolute$PC <- factor(pca_table_absolute$PC,
                             levels = paste0("PC", 9:1))  # reverses order

pca_table_absolute$model <- factor(pca_table_absolute$model,
                                levels = paste0("Model ", 3:1))  # reverses order


pca_absolute_forest <- ggplot(pca_table_absolute, aes(x = HR, y = PC, xmin = LCL, xmax = UCL, shape = model)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbarh(
    aes(group = model),
    position = position_dodge(width = 0.6),
    height = 0.2
  ) +
  geom_point(
    aes(group = model),
    position = position_dodge(width = 0.6),
    size = 2.5
  )  +
  scale_shape_manual(
    values = c(16, 17, 15),
    breaks = c("Model 1", "Model 2", "Model 3")
  ) +
  labs(
    x = "Hazard ratio (95% CI)",
    y = "Principle component (PCs)",
    shape = "Model",
    title = "Associations of PCs with lung cancer risk"
  ) +  
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
  )

pca_absolute_forest

################ Plotting PCA loadings ################ 
loadings_absolute_df <- as.data.frame(pca_result_absolute$rotation)

loadings_absolute_df$metabolite <- rownames(loadings_absolute_df)

# Convert to long format
loadings_absolute_df_long <- loadings_absolute_df %>%
  dplyr::select(metabolite, all_of(pcs)) %>%
  tidyr::pivot_longer(
    cols = all_of(pcs),
    names_to = "PC",
    values_to = "loading"
  ) 

# Set the order of metabolite group for display, 
loadings_absolute_df_long <- loadings_absolute_df_long %>%
  mutate(
    # identify lipid component
    component = case_when(
      str_detect(metabolite, "Apo|TG/PG|Omega_6_FA/Omega_3_FA|PUFA/MUFA") ~ "Lipid ratios",
      str_detect(metabolite, "Omega|DHA|PUFA|MUFA|SFA|LA") ~ "Fatty acids",
      str_detect(metabolite, "_C") ~ "C",
      str_detect(metabolite, "_FC") ~ "FC",
      str_detect(metabolite, "_CE") ~ "CE",
      str_detect(metabolite, "_TG") ~ "TG",
      str_detect(metabolite, "_PL") ~ "PL"
    ),
    component = factor(component, levels = c("C", "FC", "CE", "TG", "PL","Fatty acids","Lipid ratios")
    ),
    # Identify the direction of loadings
    direction = if_else(loading >= 0, "Positive", "Negative"),
    abs_loading = abs(loading)
  ) %>%
  arrange(component)

# Set the metabolite order for plotting
met_order <- loadings_absolute_df_long %>%
  distinct(component, metabolite) %>%
  arrange(component) %>%
  pull(metabolite)

# Reverse met_order, so the desired first metabolite appears at the top of the y-axis
loadings_absolute_df_long <- loadings_absolute_df_long %>%
  mutate(
    metabolite = factor(metabolite, levels = rev(met_order))
  )

# Plot the PCA loadings
factor_loadings_absolute <-
  ggplot(loadings_absolute_df_long, aes(x = PC, y = metabolite)) +
  geom_point(aes(size = abs_loading, color = direction), alpha = 0.8) +
  scale_size(name = "Loadings", range = c(0.1, 5)) +
  scale_color_manual(
    values = c("Positive" = "red3", "Negative" = "darkgreen"),
    name = "Direction"
  ) +
  labs(
    x = "Principal components (PCs)",
    y = "Metabolites",
    title = "The factor loadings of principal components"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(
      hjust = 0.5
    )
  )

# Combine the forest plot and the loading plot side-by-side
p <- pca_absolute_forest | factor_loadings_absolute

tiff("forest_factor_loadings_abs.tiff",
     width = 400, height = 200, units = "mm", res = 900)

p

dev.off()