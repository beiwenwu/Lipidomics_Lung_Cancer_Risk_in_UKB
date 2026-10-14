library(dplyr)
library(tidyr)
library(readxl)
library(circlize)
library(Hmisc)

# Step 1: Define the sector order — Out_label categories such as C%, FC%, CE%, TG%, PL%, FA%, and lipid ratios.
# Step 2: Define subgroups within sectors — Mid_label, such as VLDL, IDL, LDL, HDL, or none.
# Step 3: Assign positions to individual metabolites — using inner_slot, n_inner, cum_start, and ultimately xpos.
# Step 4: Initialize the circular layout — using circos.clear(), circos.par(), and the sector ranges.
# Step 5: Create the HR track — background shading, HR bars, and 95% confidence intervals.
# Step 6: Draw reference lines, Out_label, Mid_label, Inner_label
# Step 7: Load the metabolite correlation matrix, convert to long format, retain one unique metabolite pairs, and retain strong correlations such as |r| >= 0.90.
# Step 8: Map correlated metabolites to their Circos coordinates.
# Step 9: Draw the correlation links with using color for direction and width for strength. 
# Step 10: Close/save the plot.

# Below are sample code for overall analysis.
# For plotting stratified analyses results, load appropriate result and correlation matrix, use appropriate plot title and file name.

################ Make Circos Plots for Relative Metabolites ################
# Load CoxPH results for Relative Metabolites 
# Manually Added Out_label, Mid_label, and Inner_label to excel sheets
coxph_results_overall_ratio 
coxph_results_female_ratio
coxph_results_male_ratio
coxph_results_current_ratio
coxph_results_previous_ratio
coxph_results_never_ratio
coxph_results_agegroup1_ratio
coxph_results_agegroup2_ratio
coxph_results_agegroup3_ratio
coxph_results_NSCLC_ratio
coxph_results_SCLC_ratio
coxph_results_other_ratio
coxph_results_complete_case_ratio

# Define the display order of metabolite groups for out layer
group_order <- c("Cholesterol to Total Lipids %","Free Cholesterol to Total Lipids %", "Cholesteryl Esters to Total Lipids %", "Triglycerides to Total Lipids %",
                 "Phospholipids to Total Lipids %", "% in Total Fatty Acids (FA)", "Lipid Ratios")

# Define the display order of major lipoprotein classes for mid layer
mid_order <- c("VLDL","IDL","LDL","HDL")

# Prepare the relative metabolite results for plotting
# Use appropriate coxph_results_subgroup_ratio name
core_ratio <- coxph_results_overall_ratio %>%
  
  # Keep only the metabolite groups included in the predefined display order
  filter(Out_label %in% group_order) %>%
  
  # Rename the hazard ratio for simpler plotting
   rename(HR = Hazard_Ratio) %>%
  
  # Create variables indicating statistical significance, direction of association, and plotting colour
  mutate(
    signif = P_Adjusted < 0.05,
    dir    = ifelse(HR >= 1, "Up", "Down"),
    col      = ifelse(!signif, "#9CA3AF", ifelse(dir == "Up", "#D55E00", "#0072B2")),
  ) %>%
  
  # Set the plotting order of metabolite groups and lipoprotein classes
  mutate(
    Out_label = factor(Out_label, levels = group_order),
    Mid_label = factor(Mid_label, levels = mid_order)
  ) %>%
  
  # Assign an internal plotting position 
  group_by(Out_label, Mid_label) %>%
  mutate(inner_slot = row_number(), # Position of each metabolite within the subgroup
         n_inner    = n()) %>% # Total number of lipoprotein subclasses within each lipoprotein class, VLDL=6, IDL=1, LDL=3, HDL=4,FA_pct=7,FA_ratio=4
  ungroup()

# Add "none" to Mid_label for those non-lipoprotein metabolites
core_ratio <- core_ratio %>%
  mutate(
    Mid_label = if_else(is.na(Mid_label), "none", Mid_label),
    Mid_label = factor(Mid_label, levels = c(mid_order, "none"))
  )

# Define the spacing between adjacent lipoprotein classes
gap <- 0.6

# Extract the size and starting position of each lipoprotein class
mid_sizes <- core_ratio %>%
  distinct(Out_label, Mid_label, n_inner) %>%
  arrange(Out_label, Mid_label) %>%
  group_by(Out_label) %>%
  mutate(
    mid_index = match(as.character(Mid_label), c(mid_order, "none")), # Converts the mid-label into a numeric order, VLDL=1, IDL=2, LDL=3, HDL=4, none=5
    cum_start = cumsum(lag(n_inner + gap, default = 0)) # define where does each lipoprotein subgroup start within its Out_label sector, VLDL=0, IDL=6+0.6, LDL=6.6+1+0.6,HDL=8.2+3+0.6=11.8
  )

# Assign a global sequential index to each lipoprotein subgroup, this is used to alternate background shading
mid_sizes <- mid_sizes |>
  mutate(Out_label = factor(Out_label, levels = levels(core_ratio$Out_label))) |>
  arrange(Out_label, cum_start) |>
  mutate(global_i = row_number())

# Merge subgroup positioning information back into the metabolite data and calculate the final plotting position for each metabolite
core_ratio <- core_ratio %>%
  left_join(mid_sizes, by = c("Out_label","Mid_label")) %>%
  mutate(xpos = cum_start + inner_slot)

# Open a 9 × 9 inch, 900-dpi TIFF graphics device for saving the plot
# Change the name according to the coxph_results_ratio_subgroup name
tiff("coxph_results_overall_ratio.tiff", width = 9, height = 9, units = "in", res = 900)

# Clear any existing Circos plot and reset the plotting environment
circos.clear()

# Remove outer plot margins, while leaving a small top margin for the title
par(mar = c(0,0,2,0))

# Add horizontal padding so the first and last metabolites have enough space for HR bars and CI caps
sector_ranges <- core_ratio %>%
  group_by(Out_label) %>%
  summarise(
    xmin = min(xpos, na.rm = TRUE) - 0.5,
    xmax = max(xpos, na.rm = TRUE) + 0.5,
    .groups = "drop"
  )

# Set global Circos plot parameters
circos.par(
  cell.padding = c(0, 0, 0, 0), # Remove padding inside sectors/tracks
  start.degree = 90, # Start the circular plot at the top
  gap.after = c(rep(6, nlevels(core_ratio$Out_label) - 1), 10), # Add gaps between metabolite sectors
  track.margin = c(0.01, 0.01), # Set small margins between tracks
  canvas.xlim = c(-1.2, 1.2),  # Expand horizontal plotting area
  canvas.ylim = c(-1.2, 1.2),   # Expand vertical plotting area
  points.overflow.warning = FALSE
)

# Initialize sectors
circos.initialize(
  factors = sector_ranges$Out_label,
  xlim = as.matrix(
    sector_ranges[, c("xmin", "xmax")]
  )
)
# Set fixed y-axis limits
ymin <- 0.7
ymax <- 1.3


# Create the Circos track containing subgroup shading, HR estimates, and CI
circos.trackPlotRegion(
  ylim = c(ymin, ymax), bg.border = NA, track.height = 0.33, #sets the HR scale and ring height
  panel.fun = function(x, y) {
    # Identify the current Out_label sector
    si   <- as.character(get.cell.meta.data("sector.index"))     
    # Obtain the y-axis limits for the current sector
    ylim <- get.cell.meta.data("ylim")
    y1   <- ylim[1]; y2 <- ylim[2]
    
    # Extract lipoprotein subgroup information for the current sector
    ms <- mid_sizes |>
      filter(Out_label == si) |>
      arrange(cum_start)
    
    # Add alternating background shading for the lipoprotein subgroups
    if (nrow(ms) > 0) {
      for (i in seq_len(nrow(ms))) {
        x1 <- ms$cum_start[i] + 0.5
        x2 <- ms$cum_start[i] + ms$n_inner[i] + 0.5

        fill <- if (si == "Lipid Ratios") {"white"} 
                else if ((ms$global_i[i] %% 2) == 1) {"grey95"} else {"white"}        
        
        circos.rect(x1, y1, x2, y2, col = fill, border = NA)
      }
    }

    # Select metabolite-specific CoxPH results for the current sector
    idx <- core_ratio$Out_label == si
    if (!any(idx)) return()
    
    # Extract plotting positions, HRs, confidence intervals, and colors
    x0 <- as.numeric(core_ratio$xpos[idx])
    hr <- as.numeric(core_ratio$HR[idx])
    lo <- as.numeric(core_ratio$CI_Lower[idx])
    hi <- as.numeric(core_ratio$CI_Upper[idx])
    cc <- core_ratio$col[idx]
    
    # Draw vertical 95% CI lines
    circos.segments(x0, lo, x0, hi, col = "#111111", lwd = 0.6)
    
    # Draw horizontal caps at the lower and upper CI limits
    w <- 0.45
    circos.segments(x0 - w, lo, x0 + w, lo, col = "#111111", lwd = 0.6)
    circos.segments(x0 - w, hi, x0 + w, hi, col = "#111111", lwd = 0.6)
    
    # Draw the HR bar and use HR = 1 as the reference value
    yb <- pmin(1, hr); yt <- pmax(1, hr)
    bar_w <- 0.6
    circos.rect(x0 - bar_w/2, yb, x0 + bar_w/2, yt, col = cc, border = NA)
  }
)

# Draw horizontal HR reference lines (0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3) across every out_label sector
for (g in levels(core_ratio$Out_label)) {
  
  # Get the actual x-axis limits of the current sector
  xlim_g <- get.cell.meta.data(
    "xlim",
    sector.index = g)
  
  for (h in c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)) {
    circos.lines(xlim_g, c(h, h), sector.index = g,
      col = if (h == 1) "#6B7280" else "#9CA3AF",
      lty = if (h == 1) 1 else 3,
      lwd = 0.8
    )
  }
}


# Identify the first Out_label sector
first_sector <- levels(core_ratio$Out_label)[1]

xlim_first <- get.cell.meta.data(
  "xlim",
  sector.index = first_sector
)

# Add HR reference-value labels to the first Out_label sector only
for (h in c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)) {
    circos.text(x = xlim_first[1], y = h, labels = as.character(h),
                sector.index = first_sector,
                facing = "inside", adj = c(1, -0.2),
                cex = 0.6,
                col = if (h==1) "#6B7280" else "#9CA3AF")
  }


# Draw Out labels: Loops over each sector in the current track and adds its sector name
circos.trackPlotRegion(track.index = get.current.track.index(),
                       panel.fun = function(x, y) {
                         label  <- as.character(get.cell.meta.data("sector.index"))  
                         xcenter <- get.cell.meta.data("xcenter") # finds the horizontal center of that sector so the label can be centered
                         ybase   <- get.cell.meta.data("ylim")[2] * 1.41 #places the label outside the outer boundary of the HR track
                         # Draws the sector label curved around the outside of the circle
                         circos.text(
                           x = xcenter,
                           y = ybase,
                           labels = label,
                           facing = "bending.out",
                           niceFacing = TRUE,
                           adj = c(0.5, 0),
                           cex = 0.6,
                           font = 2
                         )
                       },
                       bg.border = NA
)

# Draw Mid labels
circos.trackPlotRegion(
  track.index = get.current.track.index(),
  bg.border = NA,
  panel.fun = function(x, y) {
    si    <- as.character(get.cell.meta.data("sector.index"))
    ylim  <- get.cell.meta.data("ylim")
    ytop  <- ylim[2]
    y_mid <- ytop * 1.29   
    
    ms <- mid_sizes |>
      filter(Out_label == si) |>
      arrange(Mid_label)
    
    if (nrow(ms) == 0) return()
    
    # center of each Mid block
    x_centers <- ms$cum_start + ms$n_inner / 2
    
    for (i in seq_len(nrow(ms))) {

      lab_i <- as.character(ms$Mid_label[i])
      if (lab_i == "none") next
      
      circos.text(
        x = x_centers[i],
        y = y_mid,
        labels = lab_i,
        facing = "bending.out",
        niceFacing = TRUE,
        adj = c(0.5, 0),
        cex = 0.55,
        font = 2,
        col  = "#374151"
      )
    }
  }
)

# Draw Inner labels
circos.trackPlotRegion(
  track.index = get.current.track.index(),   
  bg.border = NA,
  panel.fun = function(x, y) {
    si    <- as.character(get.cell.meta.data("sector.index"))
    ylim  <- get.cell.meta.data("ylim")
    ytop  <- ylim[2]
    
    y_lab <- ytop * 1.19
    
    # Select only metabolites belonging to the current sector and order them by Mid_label, cum_start, inner_slot
    sub <- core_ratio |>
      filter(Out_label == si) |>
      arrange(Mid_label, cum_start, inner_slot)
    
    # stop if the sector has no more metabolites
    if (nrow(sub) == 0) return()
    
    # Loop through every metabolite in the current sector, draw inner label
    for (k in seq_len(nrow(sub))) {
      lab_txt <- paste0(sub$Inner_label[k])
      circos.text(
        x = as.numeric(sub$xpos[k]),
        y = y_lab,
        labels = lab_txt,
        sector.index = si,
        facing = "clockwise",
        niceFacing = TRUE, # Automatically flip labels when needed so they remain readable rather than upside down
        adj = c(0.5, 0), # the label is horizontally centered around its xpos
        cex = 0.5, # font size
        col = "#6B7280"
      )
    }
  }
)

# Create a title for the plot
title("CoxPH Results for Relative Metabolites",font = 3)

# Convert the metabolite correlation matrix to long format
r_long <- as.data.frame(as.table(pooled_metabolite_cor_ratio)) %>%
  rename(from = Var1, to = Var2, r = Freq) %>%
  mutate(from = as.character(from),
                to   = as.character(to)) %>%
  filter(from < to)  # Retain only one copy of each metabolite pair and remove self-correlations

# Calculate the absolute value of correlation magnitude for each metabolite pair
r_long <- r_long %>%
  mutate(r_abs = abs(r))

# Create a lookup table linking each metabolite to its Circos sector
feat_map <- core_ratio %>%
  transmute(
    Metabolite = as.character(Metabolite),
    Group      = as.character(Out_label),
    xpos       = xpos
  )

# Create the set of strong correlation links to be drawn in the Circos plot.
links <- r_long %>%
  mutate(from = as.character(from),
         to   = as.character(to)) %>%
  left_join(feat_map, by = c("from" = "Metabolite")) %>%
   rename(from_group = Group, from_xpos = xpos) %>%
  left_join(feat_map, by = c("to" = "Metabolite")) %>%
   rename(to_group = Group, to_xpos = xpos) %>%
  filter(!is.na(from_group), !is.na(to_group)) %>%
  filter(r_abs >= 0.9)   # Retain only strongly correlated metabolite pairs

# Adjust color transparency 
adj_col <- function(col, alpha = 0.35) grDevices::adjustcolor(col, alpha.f = alpha)

# Draw the correlation links 
for (i in seq_len(nrow(links))) {
  fsec <- links$from_group[i];  fpos <- links$from_xpos[i]
  tsec <- links$to_group[i];    tpos <- links$to_xpos[i]
  if (is.na(fpos) || is.na(tpos)) next
  col_i <- if (links$r[i] >= 0) adj_col("#E41A1C", 0.30) else adj_col("#08306B", 0.30) # color = direction
  lwd_i <- 0.1 + 2.5 * links$r_abs[i] # line width = strength
  circos.link(
    sector.index1 = fsec, point1 = fpos,
    sector.index2 = tsec, point2 = tpos,
    col = col_i, lwd = lwd_i,
    h.ratio = 0.65, rou = 0.6
  )
}

# Legend
legend("bottomright",
       legend = c("Significant HR > 1", "Significant HR < 1","Non-significant HR", "Correlation >= 0.9", "Correlation <= -0.9"),
       fill = c("#D55E00","#0072B2","#9CA3AF", "#E41A1C","#08306B"),
       border = NA, bty = "n", cex = 0.85)


circos.clear()
dev.off()

################ Make Circos Plots for Absolute Metabolites ################
# Load CoxPH results for Absolute Metabolites 
# Manually Added Out_label, Mid_label, and Inner_label to excel sheets
coxph_results_overall_abs 
coxph_results_female_abs
coxph_results_male_abs
coxph_results_current_abs
coxph_results_previous_abs
coxph_results_never_abs
coxph_results_agegroup1_abs
coxph_results_agegroup2_abs
coxph_results_agegroup3_abs
coxph_results_NSCLC_abs
coxph_results_SCLC_abs
coxph_results_other_abs
coxph_results_complete_case_abs

##### Part 1: Draw Lipid-related metabolites ##### 
# Define the display order of metabolite groups for out layer
group_order <- c("Particle concentration","Total lipids","Cholesterol", "Free cholesterol", "Cholesteryl esters", "Triglycerides", "Phospholipids", "Particle size")

# Define the display order of major lipoprotein classes for mid layer
mid_order <- c("VLDL","IDL","LDL","HDL","Others")

# Prepare the absolute metabolite results for plotting
# Use appropriate coxph_results_abs_subgroup name
core_abs_part1 <- coxph_results_overall_abs %>%
  
  # Keep only the metabolite groups included in the predefined display order
  filter(Out_label %in% group_order) %>%
  
  # Rename the hazard ratio for simpler plotting
  rename(HR = Hazard_Ratio) %>%

  # Create variables indicating statistical significance, direction of association, and plotting colour
  mutate(
    signif = P_Adjusted < 0.05,
    dir    = ifelse(HR >= 1, "Up", "Down"),
    col      = ifelse(!signif, "#9CA3AF", ifelse(dir == "Up", "#D55E00", "#0072B2")),
  ) %>%
  
  # Set the plotting order of metabolite groups and lipoprotein classes
  mutate(
    Out_label = factor(Out_label, levels = group_order),
    Mid_label = factor(Mid_label, levels = mid_order)
  ) %>%
  
  # Assign an internal plotting position 
  group_by(Out_label, Mid_label) %>%
  mutate(inner_slot = row_number(), # Position of each metabolite within the subgroup  
         n_inner    = n()) %>% # Total number of lipoprotein subclasses within each lipoprotein class, VLDL=7,IDL=1,LDL=4,HDL=5,Others=3
  ungroup()

# Add "none" to Mid_label for those does not have a mid-label
core_abs_part1 <- core_abs_part1 %>%
  mutate(
    Mid_label = if_else(is.na(Mid_label), "none", Mid_label),
    Mid_label = factor(Mid_label, levels = c(mid_order, "none"))
  )

# Define the spacing between adjacent lipoprotein classes
gap <- 0.6

# Extract the size and starting position of each lipoprotein class
mid_sizes <- core_abs_part1 %>%
  distinct(Out_label, Mid_label, n_inner) %>%
  arrange(Out_label, Mid_label) %>%
  group_by(Out_label) %>%
  mutate(
    mid_index = match(as.character(Mid_label), c(mid_order, "none")), # Converts the mid-label into a numeric order, VLDL=1, IDL=2, LDL=3, HDL=4, Others=5, none=6
    cum_start = cumsum(lag(n_inner + gap, default = 0)) # define where does each mid-label start within its Out_label sector, VLDL=0, IDL=7+0.6, LDL=7.6+1+0.6,HDL=9.2+4+0.6,Others/TG/TP/PL=13.8+5+0.6=19.4,TC=19.4+3+0.6=23
  ) %>%
  ungroup()

# Assign a global sequential index to each lipoprotein subgroup, this is used to alternate background shading
mid_sizes <- mid_sizes |>
  mutate(Out_label = factor(Out_label, levels = levels(core_abs_part1$Out_label))) |>
  arrange(Out_label, cum_start) |>
  mutate(global_i = row_number())

# Merge subgroup positioning information back into the metabolite data and calculate the final plotting position for each metabolite
core_abs_part1 <- core_abs_part1 %>%
  left_join(mid_sizes, by = c("Out_label","Mid_label")) %>%
  mutate(xpos = cum_start + inner_slot)   # final position of each metabolite within its Out_label sector

# Open a 9 × 9 inch, 900-dpi TIFF graphics device for saving the plot
# Change the name according to the coxph_results_XXX_abs name
tiff("coxph_results_overall_abs_part1.tiff", width = 9, height = 9, units = "in", res = 900)

# Clear any existing Circos plot and reset the plotting environment
circos.clear()

# Remove outer plot margins, while leaving a small top margin for the title
par(mar = c(0,0,2,0))

# Add horizontal padding so the first and last metabolites have enough space for HR bars and CI caps
sector_ranges <- core_abs_part1 %>%
  group_by(Out_label) %>%
  summarise(
    xmin = min(xpos, na.rm = TRUE) - 0.5,
    xmax = max(xpos, na.rm = TRUE) + 0.5,
    .groups = "drop"
  )

# Set global Circos plot parameters
circos.par(
  cell.padding = c(0, 0, 0, 0), # Remove padding inside sectors/tracks
  start.degree = 90, # Start the circular plot at the top
  gap.after = c(rep(6, nlevels(core_abs_part1$Out_label) - 1), 10), # Add gaps between metabolite sectors
  track.margin = c(0.01, 0.01), # Set small margins between tracks
  canvas.xlim = c(-1.2, 1.2),  # Expand horizontal plotting area
  canvas.ylim = c(-1.2, 1.2),   # Expand vertical plotting area
  points.overflow.warning = FALSE
)

# Initialize sectors
circos.initialize(
  factors = sector_ranges$Out_label,
  xlim = as.matrix(
    sector_ranges[, c("xmin", "xmax")]
  )
)

# Set fixed y-axis limits
ymin <- 0.7
ymax <- 1.3

# Create the Circos track containing subgroup shading, HR estimates, and CI
circos.trackPlotRegion(
  ylim = c(ymin, ymax), bg.border = NA, track.height = 0.33, #sets the HR scale and ring height
  panel.fun = function(x, y) {
    # Identify the current Out_label sector
    si   <- as.character(get.cell.meta.data("sector.index"))     
    # Obtain the y-axis limits for the current sector
    ylim <- get.cell.meta.data("ylim")
    y1   <- ylim[1]; y2 <- ylim[2]
    
    # Extract lipoprotein subgroup information for the current sector
    ms <- mid_sizes |>
      filter(Out_label == si) |>
      arrange(cum_start)
    
    # Add alternating background shading for the lipoprotein subgroups
    if (nrow(ms) > 0) {
      for (i in seq_len(nrow(ms))) {
        x1 <- ms$cum_start[i] + 0.5
        x2 <- ms$cum_start[i] + ms$n_inner[i] + 0.5
        
        fill <- if ((ms$global_i[i] %% 2) == 1) {"grey95"} else {"white"}        
        
        circos.rect(x1, y1, x2, y2, col = fill, border = NA)
      }
    }
    
    # Select metabolite-specific CoxPH results for the current sector
    idx <- core_abs_part1$Out_label == si
    if (!any(idx)) return()
    
    # Extract plotting positions, HRs, confidence intervals, and colors
    x0 <- as.numeric(core_abs_part1$xpos[idx])
    hr <- as.numeric(core_abs_part1$HR[idx])
    lo <- as.numeric(core_abs_part1$CI_Lower[idx])
    hi <- as.numeric(core_abs_part1$CI_Upper[idx])
    cc <- core_abs_part1$col[idx]
    
    # Draw vertical 95% CI lines
    circos.segments(x0, lo, x0, hi, col = "#111111", lwd = 0.6)
    
    # Draw horizontal caps at the lower and upper CI limits
    w <- 0.45
    circos.segments(x0 - w, lo, x0 + w, lo, col = "#111111", lwd = 0.6)
    circos.segments(x0 - w, hi, x0 + w, hi, col = "#111111", lwd = 0.6)
    
    # Draw the HR bar and use HR = 1 as the reference value
    yb <- pmin(1, hr); yt <- pmax(1, hr)
    bar_w <- 0.6
    circos.rect(x0 - bar_w/2, yb, x0 + bar_w/2, yt, col = cc, border = NA)
  }
)

# Draw horizontal HR reference lines (0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3) across every out_label sector
for (g in levels(core_abs_part1$Out_label)) {
  
  # Get the actual x-axis limits of the current sector
  xlim_g <- get.cell.meta.data(
    "xlim",
    sector.index = g)
  
  for (h in c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)) {
    circos.lines(xlim_g, c(h, h), sector.index = g,
                 col = if (h == 1) "#6B7280" else "#9CA3AF",
                 lty = if (h == 1) 1 else 3,
                 lwd = 0.8
    )
  }
}

# Identify the first Out_label sector
first_sector <- levels(core_abs_part1$Out_label)[1]

xlim_first <- get.cell.meta.data(
  "xlim",
  sector.index = first_sector
)

# Add HR reference-value labels to the first Out_label sector only
for (h in c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)) {
  circos.text(x = xlim_first[1], y = h, labels = as.character(h),
              sector.index = first_sector,
              facing = "inside", adj = c(1, -0.2),
              cex = 0.6,
              col = if (h==1) "#6B7280" else "#9CA3AF")
}

# Draw Out labels: Loops over each sector in the current track and adds its sector name
circos.trackPlotRegion(track.index = get.current.track.index(),
                       panel.fun = function(x, y) {
                         label  <- as.character(get.cell.meta.data("sector.index"))  
                         xcenter <- get.cell.meta.data("xcenter") # finds the horizontal center of that sector so the label can be centered
                         ybase   <- get.cell.meta.data("ylim")[2] * 1.41 #places the label outside the outer boundary of the HR track
                         # Draws the sector label curved around the outside of the circle
                         circos.text(
                           x = xcenter,
                           y = ybase,
                           labels = label,
                           facing = "bending.out",
                           niceFacing = TRUE,
                           adj = c(0.5, 0),
                           cex = 0.6,
                           font = 2
                         )
                       },
                       bg.border = NA
)

# Draw Mid labels
circos.trackPlotRegion(
  track.index = get.current.track.index(),
  bg.border = NA,
  panel.fun = function(x, y) {
    si    <- as.character(get.cell.meta.data("sector.index"))
    ylim  <- get.cell.meta.data("ylim")
    ytop  <- ylim[2]
    y_mid <- ytop * 1.29   
    
    ms <- mid_sizes |>
      filter(Out_label == si) |>
      arrange(Mid_label)
    
    if (nrow(ms) == 0) return()
    
    # center of each Mid block
    x_centers <- ms$cum_start + ms$n_inner / 2
    
    for (i in seq_len(nrow(ms))) {
      
      lab_i <- as.character(ms$Mid_label[i])
      if (lab_i == "none") next
      
      circos.text(
        x = x_centers[i],
        y = y_mid,
        labels = lab_i,
        facing = "bending.out",
        niceFacing = TRUE,
        adj = c(0.5, 0),
        cex = 0.55,
        font = 2,
        col  = "#374151"
      )
    }
  }
)

# Draw Inner labels
circos.trackPlotRegion(
  track.index = get.current.track.index(),   
  bg.border = NA,
  panel.fun = function(x, y) {
    si    <- as.character(get.cell.meta.data("sector.index"))
    ylim  <- get.cell.meta.data("ylim")
    ytop  <- ylim[2]
    
    y_lab <- ytop * 1.15
    
    sub <- core_abs_part1 |>
      filter(Out_label == si) |>
      arrange(Mid_label, cum_start, inner_slot)
    
    if (nrow(sub) == 0) return()
    
    for (k in seq_len(nrow(sub))) {
      lab_txt <- paste0(sub$Inner_label[k])
      circos.text(
        x = as.numeric(sub$xpos[k]),
        y = y_lab,
        labels = lab_txt,
        sector.index = si,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0.5, 0),
        cex = 0.5,
        col = "#6B7280"
      )
    }
  }
)

# Create a title for the plot
title("CoxPH Results for Concentrations and Compositions of Main and Sub Classes of Lipoproteins",font = 3)

# Convert the metabolite correlation matrix to long format
r_long <- as.data.frame(as.table(pooled_metabolite_cor_abs)) %>%
  rename(from = Var1, to = Var2, r = Freq) %>%
  mutate(from = as.character(from),
         to   = as.character(to)) %>%
  filter(from < to)  # Retain only one copy of each metabolite pair and remove self-correlations

# Calculate the absolute value of correlation magnitude for each metabolite pair
r_long <- r_long %>%
  mutate(r_abs = abs(r))

# Create a lookup table linking each metabolite to its Circos sector
feat_map <- core_abs_part1 %>%
  transmute(
    Metabolite = as.character(Metabolite),
    Group      = as.character(Out_label),
    xpos       = xpos
  )

# Create the set of strong correlation links to be drawn in the Circos plot.
links <- r_long %>%
  mutate(from = as.character(from),
         to   = as.character(to)) %>%
  left_join(feat_map, by = c("from" = "Metabolite")) %>%
  rename(from_group = Group, from_xpos = xpos) %>%
  left_join(feat_map, by = c("to" = "Metabolite")) %>%
  rename(to_group = Group, to_xpos = xpos) %>%
  filter(!is.na(from_group), !is.na(to_group)) %>%
  filter(r_abs >= 0.99)   # Retain only strongly correlated metabolite pairs

# Adjust color transparency 
adj_col <- function(col, alpha = 0.35) grDevices::adjustcolor(col, alpha.f = alpha)

# Draw the correlation links 
for (i in seq_len(nrow(links))) {
  fsec <- links$from_group[i];  fpos <- links$from_xpos[i]
  tsec <- links$to_group[i];    tpos <- links$to_xpos[i]
  if (is.na(fpos) || is.na(tpos)) next
  col_i <- if (links$r[i] >= 0) adj_col("#E41A1C", 0.30) else adj_col("#08306B", 0.30) # color = direction
  lwd_i <- 0.1 + 2.5 * links$r_abs[i] # line width = strength
  circos.link(
    sector.index1 = fsec, point1 = fpos,
    sector.index2 = tsec, point2 = tpos,
    col = col_i, lwd = lwd_i,
    h.ratio = 0.65, rou = 0.6
  )
}

# Legend
legend("bottomright",
       legend = c("Significant HR < 1","Non-significant HR", "Correlation >= 0.9"),
       fill = c("#0072B2","#9CA3AF", "#E41A1C"),
       border = NA, bty = "n", cex = 0.85)

circos.clear()
dev.off()

##### Part 2: Draw other metabolites ##### 
group_order <- c(
  "Other lipids", "Apolipoproteins (Apo)", "Fatty acids (FA)", "Amino acids (AA)", "Glycolysis related metabolites",
  "Ketone bodies", "Fluid balance & Inflammation"
)

core_abs_part2 <- coxph_results_overall_abs %>%
  
  # Keep only the metabolite groups included in the predefined display order
  filter(Group %in% group_order) %>%
  
  # Rename the hazard ratio for simpler plotting
  rename(HR = Hazard_Ratio) %>%
  
  # Create variables indicating statistical significance, direction of association, and plotting colour
  mutate(
    signif = P_Adjusted < 0.05,
    dir    = ifelse(HR >= 1, "Up", "Down"),
    col      = ifelse(!signif, "#9CA3AF", ifelse(dir == "Up", "#D55E00", "#0072B2")),
  ) %>%
  
  # Set the plotting order of metabolite groups and lipoprotein classes
  mutate(
    Group = factor(Group, levels = group_order)
  ) %>%
  
  # Assign the final plotting position 
  group_by(Group) %>%
  mutate(xpos = row_number()) %>%
  ungroup()

# Open a 9 × 9 inch, 900-dpi TIFF graphics device for saving the plot
# Change the name according to the coxph_results_XXX_abs name
tiff("coxph_results_overall_abs_part2.tiff", width = 9, height = 9, units = "in", res = 900)

# Clear any existing Circos plot and reset the plotting environment
circos.clear()

# Remove outer plot margins, while leaving a small top margin for the title
par(mar = c(0,0,2,0))

# Add horizontal padding so the first and last metabolites have enough space for HR bars and CI caps
sector_ranges <- core_abs_part2 %>%
  group_by(Group) %>%
  summarise(
    xmin = min(xpos, na.rm = TRUE) - 0.5,
    xmax = max(xpos, na.rm = TRUE) + 0.5,
    .groups = "drop"
  )

# Set global Circos plot parameters
circos.par(
  cell.padding = c(0, 0, 0, 0), # Remove padding inside sectors/tracks
  start.degree = 90, # Start the circular plot at the top
  gap.after = c(rep(6, nlevels(core_abs_part2$Group) - 1), 10), # Add gaps between metabolite sectors
  track.margin = c(0.01, 0.01), # Set small margins between tracks
  canvas.xlim = c(-1.2, 1.2),  # Expand horizontal plotting area
  canvas.ylim = c(-1.2, 1.2),   # Expand vertical plotting area
  points.overflow.warning = FALSE
)

# Initialize sectors
circos.initialize(
  factors = core_abs_part2$Group, 
  xlim = as.matrix(
    sector_ranges[, c("xmin", "xmax")]
    )
  )

# Set fixed y-axis limits
ymin <- 0.7
ymax <- 1.3

# Create the Circos track containing subgroup shading, HR estimates, and CI
circos.trackPlotRegion(
  ylim = c(ymin, ymax), bg.border = NA, track.height = 0.33, #sets the HR scale and ring height
  panel.fun = function(x, y) {
    # Identify the current group sector
    si   <- as.character(get.cell.meta.data("sector.index"))     
    # Obtain the y-axis limits for the current sector
    ylim <- get.cell.meta.data("ylim")
    y1   <- ylim[1]; y2 <- ylim[2]
    xlim <- get.cell.meta.data("xlim")
    x1 <- xlim[1] - 0.5; x2 <- xlim[2] + 0.5
        
    # alternate by sector index
    sector_levels <- unique(as.character(core_abs_part2$Group))
    sector_i <- match(si, sector_levels)
    fill <- if ((sector_i %% 2) == 1) "white" else "grey95"
        
    circos.rect(x1, y1, x2, y2, col = fill, border = NA)
    
    # Select metabolite-specific CoxPH results for the current sector
    idx <- core_abs_part2$Group == si
    if (!any(idx)) return()
    
    # Extract plotting positions, HRs, confidence intervals, and colors
    x0 <- as.numeric(core_abs_part2$xpos[idx])
    hr <- as.numeric(core_abs_part2$HR[idx])
    lo <- as.numeric(core_abs_part2$CI_Lower[idx])
    hi <- as.numeric(core_abs_part2$CI_Upper[idx])
    cc <- core_abs_part2$col[idx]
    
    # Draw vertical 95% CI lines
    circos.segments(x0, lo, x0, hi, col = "#111111", lwd = 0.6)
    
    # Draw horizontal caps at the lower and upper CI limits
    w <- 0.2
    circos.segments(x0 - w, lo, x0 + w, lo, col = "#111111", lwd = 0.6)
    circos.segments(x0 - w, hi, x0 + w, hi, col = "#111111", lwd = 0.6)
    
    # Draw the HR bar and use HR = 1 as the reference value
    y_bottom <- pmin(1, hr)
    y_top    <- pmax(1, hr)
    bar_w    <- 0.3
    
    circos.rect(x0 - bar_w/2, y_bottom, x0 + bar_w/2, y_top,
                col = cc, border = NA)
  }
)

# Draw horizontal HR reference lines (0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3) across every out_label sector
for (g in levels(core_abs_part2$Group)) {
  
  # Get the actual x-axis limits of the current sector
  xlim_g <- get.cell.meta.data(
    "xlim",
    sector.index = g)
  
  for (h in c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)) {
    circos.lines(xlim_g, c(h, h), sector.index = g,
                 col = if (h == 1) "#6B7280" else "#9CA3AF",
                 lty = if (h == 1) 1 else 3,
                 lwd = 0.8
    )
  }
}

# Identify the first Out_label sector
first_sector <- levels(core_abs_part2$Group)[1]

xlim_first <- get.cell.meta.data(
  "xlim",
  sector.index = first_sector
)

# Add HR reference-value labels to the first group sector only
for (h in c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)) {
  circos.text(x = xlim_first[1], y = h, labels = as.character(h),
              sector.index = first_sector,
              facing = "inside", adj = c(1, -0.2),
              cex = 0.6,
              col = if (h==1) "#6B7280" else "#9CA3AF")
}

# Draw Out labels: Loops over each sector in the current track and adds its sector name
circos.trackPlotRegion(track.index = get.current.track.index(),
                       panel.fun = function(x, y) {
                         label  <- as.character(get.cell.meta.data("sector.index")) 
                         xcenter <- get.cell.meta.data("xcenter") # finds the horizontal center of that sector so the label can be centered
                         ybase   <- get.cell.meta.data("ylim")[2] * 1.41 #places the label outside the outer boundary of the HR track
                         # Draws the sector label curved around the outside of the circle
                         circos.text(
                           x = xcenter,
                           y = ybase,
                           labels = label,
                           facing = "bending.out",
                           niceFacing = TRUE,
                           adj = c(0.5, 0),
                           cex = 0.6,
                           font = 2
                         )
                       },
                       bg.border = NA
)

# Draw Inner labels
circos.trackPlotRegion(
  track.index = get.current.track.index(),   
  bg.border = NA,
  panel.fun = function(x, y) {
    si    <- as.character(get.cell.meta.data("sector.index"))
    ylim  <- get.cell.meta.data("ylim")
    ytop  <- ylim[2]
    
    y_lab <- ytop * 1.19
    
    sub <- core_abs_part2 |>
      filter(Group == si) |>
      arrange(xpos)
    
    if (nrow(sub) == 0) return()
    
    for (k in seq_len(nrow(sub))) {
      lab_txt <- paste0(sub$Description[k])
      circos.text(
        x = as.numeric(sub$xpos[k]),
        y = y_lab,
        labels = lab_txt,
        sector.index = si,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0.5, 0),
        # cex = cex_i,
        cex = 0.5,
        col = "#6B7280"
      )
    }
  }
)

title("Coxph results for other lipid-related and low-molecular weight metabolites",font = 3)

# Convert the metabolite correlation matrix to long format
r_long <- as.data.frame(as.table(pooled_metabolite_cor_abs)) %>%
  rename(from = Var1, to = Var2, r = Freq) %>%
  mutate(from = as.character(from),
         to   = as.character(to)) %>%
  filter(from < to)  # Retain only one copy of each metabolite pair and remove self-correlations

# Calculate the absolute value of correlation magnitude for each metabolite pair
r_long <- r_long %>%
  mutate(r_abs = abs(r))

# Create a lookup table linking each metabolite to its Circos sector
feat_map <- core_abs_part2 %>%
  transmute(
    Metabolite = as.character(Metabolite),
    Group      = as.character(Group),
    xpos       = xpos
  )

# Create the set of strong correlation links to be drawn in the Circos plot.
links <- r_long %>%
  mutate(from = as.character(from),
         to   = as.character(to)) %>%
  left_join(feat_map, by = c("from" = "Metabolite")) %>%
  rename(from_group = Group, from_xpos = xpos) %>%
  left_join(feat_map, by = c("to" = "Metabolite")) %>%
  rename(to_group = Group, to_xpos = xpos) %>%
  filter(!is.na(from_group), !is.na(to_group)) %>%
  filter(r_abs >= 0.8)   # Retain only strongly correlated metabolite pairs

# Adjust color transparency 
adj_col <- function(col, alpha = 0.35) grDevices::adjustcolor(col, alpha.f = alpha)

# Draw the correlation links 
for (i in seq_len(nrow(links))) {
  fsec <- links$from_group[i];  fpos <- links$from_xpos[i]
  tsec <- links$to_group[i];    tpos <- links$to_xpos[i]
  if (is.na(fpos) || is.na(tpos)) next
  col_i <- if (links$r[i] >= 0) adj_col("#E41A1C", 0.30) else adj_col("#08306B", 0.30) # color = direction
  lwd_i <- 0.1 + 2.5 * links$r_abs[i] # line width = strength
  circos.link(
    sector.index1 = fsec, point1 = fpos,
    sector.index2 = tsec, point2 = tpos,
    col = col_i, lwd = lwd_i,
    h.ratio = 0.65, rou = 0.6
  )
}

# Legend
legend("bottomright",
       legend = c("Significant HR > 1", "Significant HR < 1", "Non-significant HR","Correlation >= 0.8"),
       fill = c("#D55E00", "#0072B2", "#9CA3AF", "#E41A1C"),
       border = NA, bty = "n", cex = 0.85)

circos.clear()
dev.off()
