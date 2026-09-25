###################################################################
# Functions to Generate Concentration Response Curves from 
# Hierarchical Generalized Additive Models using the 
# Linear predictor matrix

###################################################################


# hgam total PM2.5 exposure response curve prediction function

get_er_hier_gam_simple_bam <- function(model, data, pm_grid, who_ref = 9) {
  
  # ---- Extract city levels directly from the fitted model ----
  model_cities <- levels(model$model$city)
  if (is.null(model_cities)) {
    stop("Could not find factor levels for 'city' in the provided model object.")
  }
  
  # ---- Check for non-finite values ----
  for (v in c("temp03", "dewpt03", "t", "pm01")) {
    if (v %in% names(data)) {
      n_bad <- sum(!is.finite(data[[v]]), na.rm = TRUE)
      if (n_bad > 0) cat("WARNING:", n_bad, "non-finite values in", v, "\n")
    }
  }
  
  # ---- Reference row: median for ALL continuous, modal for categoricals ----
  ref_row <- data %>%
    summarise(
      temp03            = median(temp03,            na.rm = TRUE),
      dewpt03           = median(dewpt03,              na.rm = TRUE),
      t                 = median(t,                 na.rm = TRUE),  
      season_group      = names(sort(table(season_group),      decreasing = TRUE))[1],
      city_season_group = names(sort(table(city_season_group), decreasing = TRUE))[1],
      dow               = dow[1],    # preserves class
      is_holiday        = 0L,
      pcap01            = 0L,
      city              = model_cities[1],  # SAFE: Uses valid level from model
      icd_version_flag  = "ICD9"
    )
  
  # ---- Align factors to model ----
  ref_row$city <- factor(ref_row$city, levels = model_cities)
  
  if (is.factor(model$model$season_group))
    ref_row$season_group <- factor(ref_row$season_group, levels = 
                                     levels(model$model$season_group))
  if (is.factor(model$model$city_season_group))
    ref_row$city_season_group <- factor(ref_row$city_season_group, 
                                        levels = levels(model$model$city_season_group))
  if (is.factor(model$model$dow))
    ref_row$dow <- factor(ref_row$dow, levels = levels(model$model$dow))
  
  # ---- Verify no NAs in ref_row ----
  na_counts <- sapply(ref_row, function(x) sum(is.na(x)))
  if (any(na_counts > 0)) {
    cat("WARNING: NAs in ref_row:\n"); print(na_counts[na_counts > 0])
  }
  
  # ---- Build prediction frame and lpmatrix ----
  pred_df <- crossing(pm01 = pm_grid, ref_row)
  Xp      <- predict(model, newdata = pred_df, type = "lpmatrix",
                     newdata.guaranteed = TRUE, discrete = FALSE)
  
  # ---- Inspect columns ----
  col_names <- colnames(Xp)
  
  # ---- Identify columns (FIXED REGEX TYPO HERE) ----
  pop_pm_cols    <- grep("^s\\(pm01\\)\\.",                       col_names) # Exact match
  city_pm_cols   <- grep("s\\(pm01,city\\)",                       col_names)
  city_conf_cols <- grep("s\\(temp03,city\\)|s\\(dewpt03,city\\)",    col_names)
  city_re_cols   <- grep("^s\\(city\\)",                           col_names)
  city_t_cols    <- grep("s\\(t,city_season_group\\)",             col_names)
  city_fe_cols   <- grep(":city|city:",                            col_names)
  zero_cols      <- unique(c(city_pm_cols, city_conf_cols, city_re_cols,
                             city_t_cols,  city_fe_cols))
  
  cat("\nPop PM cols:", length(pop_pm_cols),
      "| Zero cols:", length(zero_cols), "\n\n")
  
  if(length(pop_pm_cols) == 0) {
    stop("Error: Main effect smooth 's(pm01)' not found in model matrix columns.")
  }
  
  # ---- Zero city-specific columns for Global Population Curve ----
  Xp_pop              <- Xp
  Xp_pop[, zero_cols] <- 0
  
  # ---- Population PM submatrix ----
  Xp_pm   <- Xp_pop[, pop_pm_cols, drop = FALSE]
  coef_pm <- coef(model)[pop_pm_cols]
  vcov_pm <- vcov(model)[pop_pm_cols, pop_pm_cols]
  
  ref_idx <- which.min(abs(pm_grid - who_ref))
  D       <- sweep(Xp_pm, 2, Xp_pm[ref_idx, ], FUN = "-")
  log_rr  <- as.vector(D %*% coef_pm)
  se      <- sqrt(diag(D %*% vcov_pm %*% t(D)))
  
  # ---- Population ER curve ----
  er_curve <- tibble(
    pm     = pm_grid,
    log_rr = log_rr,
    se     = se,
    rr     = exp(log_rr),
    rr_lo  = exp(log_rr - 1.96 * se),
    rr_hi  = exp(log_rr + 1.96 * se)
  )
  
  # ---- City-specific curves (FIXED LOGICAL CURVE MATH HERE) ----
  city_curves <- map_dfr(model_cities, function(cn) {
    pred_city <- crossing(pm01 = pm_grid, ref_row) %>%
      mutate(city = factor(cn, levels = model_cities))
    
    Xp_city <- predict(model, newdata = pred_city, type = "lpmatrix",
                       newdata.guaranteed = TRUE, discrete = FALSE)
    
    # Keep city PM random smooth components; zero out unrelated confounders
    Xp_city[, unique(c(city_conf_cols, city_re_cols,
                       city_t_cols, city_fe_cols))] <- 0
    
    # TRUE CITY CURVE = Global Smooth Matrix + Random Smooth Matrix Deviation
    city_total_smooth_cols <- c(pop_pm_cols, city_pm_cols)
    
    Xp_city_pm <- Xp_city[, city_total_smooth_cols, drop = FALSE]
    coef_city  <- coef(model)[city_total_smooth_cols]
    
    D_city     <- sweep(Xp_city_pm, 2, Xp_city_pm[ref_idx, ], FUN = "-")
    log_rr_c   <- as.vector(D_city %*% coef_city)
    
    tibble(pm = pm_grid, city = cn,
           log_rr = log_rr_c, rr = exp(log_rr_c))
  })
  
  list(
    er_curve    = er_curve,
    city_curves = city_curves,
    pop_pm_cols = pop_pm_cols,
    zero_cols   = zero_cols,
    D           = D,
    coef_pm     = coef_pm,
    vcov_pm     = vcov_pm
  )
}



# Function for predicting PM2.5 C-R curves from PM PCAP interaction HGAM models

get_er_hier_gam_simple_bam_v2 <- function(model, data, pm_grid, who_ref = 9, pcap_val = "0") {
  
  # ---- Extract factor levels directly from the fitted model ----
  model_cities <- levels(model$model$city)
  pcap_levels  <- levels(model$model$pcap01)
  
  if (is.null(model_cities) || is.null(pcap_levels)) {
    stop("Could not find factor levels for 'city' or 'pcap01' in the model object.")
  }
  
  # ---- Check for non-finite values ----
  for (v in c("temp03", "rh03", "t", "pm01")) {
    if (v %in% names(data)) {
      n_bad <- sum(!is.finite(data[[v]]), na.rm = TRUE)
      if (n_bad > 0) cat("WARNING:", n_bad, "non-finite values in", v, "\n")
    }
  }
  
  # ---- Reference row: median for ALL continuous, modal for categoricals ----
  ref_row <- data %>%
    summarise(
      temp03            = median(temp03,            na.rm = TRUE),
      dewpt03           = median(dewpt03,              na.rm = TRUE),
      t                 = median(t,                 na.rm = TRUE),  
      season_group      = names(sort(table(season_group),      decreasing = TRUE))[1],
      city_season_group = names(sort(table(city_season_group), decreasing = TRUE))[1],
      dow               = dow[1],    
      is_holiday        = 0L,
      pcap01            = pcap_val,   # Stratification target
      city              = model_cities[1],  
      city_pcap         = paste(model_cities[1], pcap_val, sep = "."), # New interaction variable
      icd_version_flag  = "ICD9"
    )
  
  # ---- Align factors strictly to model structures ----
  ref_row$city      <- factor(ref_row$city,      levels = model_cities)
  ref_row$pcap01    <- factor(ref_row$pcap01,    levels = pcap_levels)
  ref_row$city_pcap <- factor(ref_row$city_pcap, levels = levels(model$model$city_pcap))
  
  if (is.factor(model$model$season_group))
    ref_row$season_group <- factor(ref_row$season_group, levels = levels(model$model$season_group))
  if (is.factor(model$model$city_season_group))
    ref_row$city_season_group <- factor(ref_row$city_season_group, levels = levels(model$model$city_season_group))
  if (is.factor(model$model$dow))
    ref_row$dow <- factor(ref_row$dow, levels = levels(model$model$dow))
  
  # ---- Build dummy prediction matrix to scrape column names ----
  pred_df <- crossing(pm01 = pm_grid, ref_row)
  Xp      <- predict(model, newdata = pred_df, type = "lpmatrix",
                     newdata.guaranteed = TRUE, discrete = FALSE)
  col_names <- colnames(Xp)
  
  # ---- Identify columns via updated regex matching ----
  # Captures the specific population curve for the chosen PCAP level
  pop_pm_cols    <- grep(paste0("^s\\(pm01\\):pcap01", pcap_val, "\\."), col_names) 
  # Captures the new interaction random smooths
  city_pm_cols   <- grep("s\\(pm01,city_pcap\\)",                       col_names)
  
  # Legacy tracking elements (kept for downstream object list structure)
  city_conf_cols <- grep("s\\(temp03,city\\)|s\\(rh03,city\\)",    col_names)
  city_t_cols    <- grep("s\\(t,city_season_group\\)",             col_names)
  city_fe_cols   <- grep("^city|:city|city:",                      col_names)
  zero_cols      <- unique(c(city_pm_cols, city_conf_cols, city_t_cols, city_fe_cols))
  
  cat("Targeting PCAP Level:", pcap_val, "\n")
  cat("Pop PM columns found:", length(pop_pm_cols), "| Total City columns found:", length(city_pm_cols), "\n\n")
  
  if(length(pop_pm_cols) == 0) {
    stop("Error: Population smooth columns for this specific pcap_val were not found.")
  }
  
  # ---- Population ER curve (Isolated Slicing) ----
  Xp_pm   <- Xp[, pop_pm_cols, drop = FALSE]
  coef_pm <- coef(model)[pop_pm_cols]
  vcov_pm <- vcov(model)[pop_pm_cols, pop_pm_cols]
  
  ref_idx <- which.min(abs(pm_grid - who_ref))
  D       <- sweep(Xp_pm, 2, Xp_pm[ref_idx, ], FUN = "-")
  log_rr  <- as.vector(D %*% coef_pm)
  se      <- sqrt(diag(D %*% vcov_pm %*% t(D)))
  
  er_curve <- tibble(
    pm     = pm_grid,
    log_rr = log_rr,
    se     = se,
    rr     = exp(log_rr),
    rr_lo  = exp(log_rr - 1.96 * se),
    rr_hi  = exp(log_rr + 1.96 * se)
  )
  
  # ---- City-specific curves for this PCAP level ----
  city_curves <- map_dfr(model_cities, function(cn) {
    pred_city <- crossing(pm01 = pm_grid, ref_row) %>%
      mutate(
        city      = factor(cn, levels = model_cities),
        pcap01    = factor(pcap_val, levels = pcap_levels),
        city_pcap = factor(paste(cn, pcap_val, sep = "."), levels = levels(model$model$city_pcap))
      )
    
    Xp_city <- predict(model, newdata = pred_city, type = "lpmatrix",
                       newdata.guaranteed = TRUE, discrete = FALSE)
    
    # Isolate total smooth columns needed for this specific city curve
    city_total_smooth_cols <- c(pop_pm_cols, city_pm_cols)
    Xp_city_pm             <- Xp_city[, city_total_smooth_cols, drop = FALSE]
    coef_city              <- coef(model)[city_total_smooth_cols]
    
    # Mathematically isolate the net curve centered at your reference value
    D_city     <- sweep(Xp_city_pm, 2, Xp_city_pm[ref_idx, ], FUN = "-")
    log_rr_c   <- as.vector(D_city %*% coef_city)
    
    tibble(pm = pm_grid, city = cn, pcap = pcap_val,
           log_rr = log_rr_c, rr = exp(log_rr_c))
  })
  
  list(
    er_curve    = er_curve,
    city_curves = city_curves,
    pop_pm_cols = pop_pm_cols,
    zero_cols   = zero_cols,
    D           = D,
    coef_pm     = coef_pm,
    vcov_pm     = vcov_pm
  )
}








































