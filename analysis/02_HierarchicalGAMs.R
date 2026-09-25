#Code to Replicate Analysis Part 2: Hierarchical Generalized Additive Models

#First is models with no interaction to estimate association with total PM2.5 adjusting for PCAP status
#Second are models with factor by smooth to investigate PM2.5 and PCAP interaction
#Third are models with ordered factor coding of the PCAP indicator to assess whether the curves differ significantly

#Libraries
library(gratia)
library(tidyverse)
library(splines)
library(dlnm)
library(mgcv)
library(patchwork)
library(metafor)
library(broom)
library(conflicted)
library(stringr)
library(ggthemes)





# HGAMs

# 1. Setup Threading
n_threads <- 2
cat("Starting BAM with", n_threads, "threads (discrete method)...\n")

allcity_test3$city<-factor(allcity_test3$city,levels = c("Sacramento","Modesto","Fresno","Visalia","Bakersfield","Las Vegas","Reno","Salt Lake City","Provo"))

cities<-dput(unique(allcity_test3$city))

allcity_test3 <- allcity_test3 %>% 
  mutate(across(c('pcap01', 'dow', 'city', 'season_group', 'city_season_group'), as.factor))


### trim pm

test99<-allcity_test3 %>% 
  filter(pm01<75)%>% 
  mutate(icd_version_flag=as.factor(ifelse(Date<"2015-10-01","ICD9","ICD10")))

n_threads<-2


test99_clean <- subset(test99, 
                       !is.na(resp_prin) & 
                         !is.na(pm01) & 
                         !is.na(city) & 
                         !is.na(temp03) & 
                         !is.na(rh03) & 
                         !is.na(t)
)
head(test99_clean)


## HGAM PM Function
# Function with outcome parameter
# Includes parametric main effects and interactions (Clean intercept control)
# And global population-level smooth
hgam_pm_fun <- function(df, outcome) {
  formula_str <- paste(outcome, 
                       "~0 + city * dow + city:is_holiday + s(pm01, k=8, bs='ps', m=2) + s(temp03, city, bs='fs') + city:pcap01+
    s(dewpt03, city, bs='fs') + s(t,city,bs='fs',k=5,m=1) +   city:season_group")
  
  bam(as.formula(formula_str), 
      family = quasipoisson, 
      data = df,method="fREML",discrete=T,nthreads=2)
}

#Run HGAM with PM adjusting for PCAP for all outcomes
hgam_pm_results <- map(outcomes, ~hgam_pm_fun(test99_clean, .x))

#Name List
names(hgam_pm_results)<-outcome_labs

#All model summaries
map(hgam_pm_results,~summary(.x))

# hgam exposure response curve prediction function

get_er_hier_gam_simple_bam <- function(model, data, pm_grid, who_ref = 9) {
  
  # ---- Extract city levels directly from the fitted model ----
  model_cities <- levels(model$model$city)
  if (is.null(model_cities)) {
    stop("Could not find factor levels for 'city' in the provided model object.")
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
      dewpt03              = median(dewpt03,              na.rm = TRUE),
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
    ref_row$season_group <- factor(ref_row$season_group, levels = levels(model$model$season_group))
  if (is.factor(model$model$city_season_group))
    ref_row$city_season_group <- factor(ref_row$city_season_group, levels = levels(model$model$city_season_group))
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



## Get C-R curves for all outcomes
er_all <- map(
  hgam_pm_results,
  ~get_er_hier_gam_simple_bam(
    model   = .x,
    data    = test99_clean,
    pm_grid = pm_grid,
    who_ref = 9
  )
)

#Combine into single plot
er_combined <- imap_dfr(er_all, ~.x$er_curve %>% mutate(Outcome = .y)) %>%
  mutate(Outcome = factor(Outcome, levels = outcome_labs))

ggplot(er_combined, aes(pm, rr)) +
  geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi),
              alpha = 0.2, fill = "steelblue") +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  labs(x = "PM2.5 (µg/m³)", y = "RR vs 9 µg/m³",
       title = "Pooled C-R Curves") +
  theme_bw() +
  coord_cartesian(ylim = c(0.95, 1.2)) +
  facet_wrap(~Outcome, ncol = 3)


### Get 19 vs 9 contrasts

contrast_fun<-function(model){
  pred<-crosspred("pm01",model=model,at=19,cen=9,from=0,to=75)
  RR<-cbind(pred$matRRfit,pred$matRRlow,pred$matRRhigh)
}


map(hgam_pm_results,~contrast_fun(.x))





## HGAM PM PCAP Function
# Function with outcome parameter
# Includes parametric main effects and interactions (Clean intercept control)
# And global population-level smooth
hgam_pm_pcap_fun <- function(df, outcome) {
  formula_str <- paste(outcome, 
                       "~0 + city * dow + city:is_holiday + s(pm01, k=8, bs='ps', m=2,by=pcap01) + s(temp03, city, bs='fs') + city:pcap01+
    s(dewpt03, city, bs='fs') + s(t,city,bs='fs',k=5,m=1) +   city:season_group")
  
  bam(as.formula(formula_str), 
      family = quasipoisson, 
      data = df,method="fREML",discrete=T,nthreads=2)
}

#Run HGAM for all outcomes
hgam_pm_pcap_results <- map(outcomes, ~hgam_pm_pcap_fun(test99_clean, .x))

#Name List
names(hgam_pm_pcap_results)<-outcome_labs

#All model summaries
map(hgam_pm_pcap_results,~summary(.x))


# function for predicting PM C-R curves from PM PCAP hgam models

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


### Get C-R Curves
er_all_pm_pcap_0 <- map(
  hgam_pm_pcap_results,
  ~get_er_hier_gam_simple_bam_v2(
    model   = .x,
    data    = test99_clean,
    pm_grid = pm_grid,
    who_ref = 9,
    pcap_val = "0"
  )
)

er_all_pm_pcap_1 <- map(
  hgam_pm_pcap_results,
  ~get_er_hier_gam_simple_bam_v2(
    model   = .x,
    data    = test99_clean,
    pm_grid = pm_grid,
    who_ref = 9,
    pcap_val = "1"
  )
)

#Combine into single plot
er_combined_pcap0 <- imap_dfr(er_all_pm_pcap_0, ~.x$er_curve %>% mutate(Outcome = .y)) %>%
  mutate(Outcome = factor(Outcome, levels = outcome_labs))
er_combined_pcap1 <- imap_dfr(er_all_pm_pcap_1, ~.x$er_curve %>% mutate(Outcome = .y)) %>%
  mutate(Outcome = factor(Outcome, levels = outcome_labs))


# 2. Combine the population curves into a single dataframe
pm_pcap_combined_er_curves <- bind_rows(
  er_combined_pcap0 %>% mutate(group = "0"),
  er_combined_pcap1     %>% mutate(group = "1")
)

# 3. Generate the plot

pm_pcap_combined_er_curves %>%
  filter(group %in% c("0", "1")) %>%
  mutate(group = if_else(group == "0", "Non-PCAP", "PCAP"),
         group = factor(group, levels = c("Non-PCAP", "PCAP"))) %>%
  ggplot(aes(pm, rr, color = group, fill = group, linetype = group)) +
  geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 1,       linetype = "dashed", color = "gray40") +
  #geom_vline(xintercept = who_ref, linetype = "dotted", color = "red") +
  scale_linetype_manual(values = c("Non-PCAP" = "solid", "PCAP" = "dashed")) +
  labs(x = "PM2.5 (µg/m³)", y = "RR vs 9 µg/m³",title="Pooled C-R Curves",
       color = NULL, fill = NULL, linetype = NULL) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=10))+
  xlim(0,77)+
  facet_wrap(~Outcome,ncol=3)





## Ordered Factor Curve Differences

test99_clean$pcap01_of<-factor(test99_clean$pcap01,ordered = T,levels=c("0","1"))
contrasts(test99_clean$pcap01_of) <- "contr.treatment"

hgam_pm_pcap_of_fun <- function(df, outcome) {
  formula_str <- paste(outcome, 
                       "~0 + city * dow + city:is_holiday +s(pm01, k=8, bs='ps') + s(pm01, k=8, bs='ps', m=2,by=pcap01_of) + s(temp03, city, bs='fs') + city:pcap01_of+
    s(dewpt03, city, bs='fs') + s(t,city,bs='fs',k=5,m=1) +   city:season_group")
  
  bam(as.formula(formula_str), 
      family = quasipoisson, 
      data = df,method="fREML",discrete=T,nthreads=2)
}

#Run HGAM for all outcomes
hgam_pm_pcap_of_results <- map(outcomes, ~hgam_pm_pcap_of_fun(test99_clean, .x))

#Name List
names(hgam_pm_pcap_of_results)<-outcome_labs

#All model summaries
pm_pcap_of_summaries<-map(hgam_pm_pcap_of_results,~summary(.x))


pm_pcap_of_summary_tables<-map(hgam_pm_pcap_of_results,~summary(.x)$s.table)
names(pm_pcap_of_summary_tables)<-outcome_labs

pm_pcap_of_summary_tables













