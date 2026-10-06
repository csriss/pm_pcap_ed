#Code to Replicate Analysis Part 1: City Specific GAMs


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

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)

#Load Analysis Dataset

allcity_test2<-readRDS("allcity_pm_cap_ed_26aug2026.rds")


#Convert other analysis Variables to factors
allcity_test3<-allcity_test2 %>% 
  mutate_at( c('pcap01','dow','season_group','city_season_group'), as.factor) 


allcity_test3$city<-factor(allcity_test3$city,levels = c("Sacramento","Modesto","Fresno","Visalia","Bakersfield","Las Vegas","Reno","Salt Lake City","Provo"))


#Define Main Outcomes

outcomes<-c("resp_prin","lad_prin","uad_prin","asthma_attack_prin","bronchitis_prin",
            "pneumonia_prin")

outcome_labs<-c("RD","LAD","UAD","Asthma","Bronchitis","Pneumonia")

outcome_df<-data.frame(outcome=outcomes,outcome_labs=outcome_labs)


# ---------------------------------------------------------------------------
# Multivariate pooling of PM2.5 + PCAP estimates across respiratory outcomes
#
# Four models, same pipeline:
#   main              : pm01 + pcap01              (binary PCAP; primary)
#   interaction       : pm01 * pcap01              (binary PCAP; slopes by status)
#   main_phase        : pm01 + pcap_phase_3l       (NoPCAP / Build / Persist)
#   interaction_phase : pm01 * pcap_phase_3l       (slopes by phase)
#
# Stage 1: city-specific quasi-Poisson GAMs
# Stage 2: multivariate random-effects meta-analysis (mixmeta), then contrasts
# ---------------------------------------------------------------------------

outcomes <- c("resp_prin", "lad_prin", "uad_prin", "asthma_attack_prin",
              "bronchitis_prin", "pneumonia_prin")

pm_increment <- 10   # slopes reported per 10 ug/m3
pm_center    <- 9   # common centering value for PM2.5 (ug/m3); makes the
# interaction models' PCAP / phase terms = contrasts at this PM

phase_levels <- c("NoPCAP", "Build", "Persist")   # first level = reference

# ---- Data prep ---------------------------------------------------------------
# pcap01 numeric 0/1 so the coefficient is named "pcap01" (not "pcap011").
# pcap_phase_3l a factor with NoPCAP as reference, so its coefficients are
# "pcap_phase_3lBuild" and "pcap_phase_3lPersist".
# pm01c = PM2.5 centered at pm_center. Centering does not change any slope,
# only the meaning of the PCAP / phase main effects in the interaction models.

prep_data <- function(d) {
  if ("pcap01" %in% names(d)) {
    d$pcap01 <- as.numeric(as.character(d$pcap01))
  }
  if ("pcap_phase_3l" %in% names(d)) {
    d$pcap_phase_3l <- factor(as.character(d$pcap_phase_3l), levels = phase_levels)
  }
  d$pm01c <- d$pm01 - pm_center
  d
}

covariates <- "s(temp03, bs = 'cr', k = 4) + s(dewpt03, bs = 'cr',k=4) + s(t, bs = 'cr', k = 5, m = 2) + season_group + dow + is_holiday"

covariates

# ---- Model specifications ----------------------------------------------------
# coef_names: coefficients pooled jointly (order = columns of contrasts)
# contrasts:  rows = quantities to report
# slope_rows: rows that are per-ug/m3 slopes (get scaled by pm_increment)
# tests:      joint Wald tests on the pooled coefficients (named coef sets)
# bscov:      between-city covariance; "unstr" falls back to "diag"

B <- "pcap_phase_3lBuild"
P <- "pcap_phase_3lPersist"

model_specs <- list(
  
  main = list(
    rhs        = "pm01c + pcap01",
    coef_names = c("pm01c", "pcap01"),
    contrasts  = rbind(
      PM   = c(1, 0),     # PM2.5 mass, adjusted for PCAP
      PCAP = c(0, 1)      # PCAP vs non-PCAP, adjusted for PM2.5 mass
    ),
    slope_rows = "PM",
    tests      = list(),
    bscov      = "unstr"  # 3 between-city params
  ),
  
  interaction = list(
    rhs        = "pm01c * pcap01",
    coef_names = c("pm01c", "pcap01", "pm01c:pcap01"),
    contrasts  = rbind(
      noPCAP      = c(1, 0, 0),   # PM slope, non-PCAP days
      PCAP        = c(1, 0, 1),   # PM slope, PCAP days
      interaction = c(0, 0, 1),   # difference in slopes
      pcap_shift  = c(0, 1, 0)    # PCAP contrast at PM = pm_center
    ),
    slope_rows = c("noPCAP", "PCAP", "interaction"),
    tests      = list(),
    bscov      = "diag"           # 3x3 unstructured = 6 params; too many for 9 cities
  ),
  
  main_phase = list(
    rhs        = "pm01c + pcap_phase_3l",
    coef_names = c("pm01c", B, P),
    contrasts  = rbind(
      PM              = c(1,  0, 0),   # PM2.5 mass, adjusted for phase
      Build           = c(0,  1, 0),   # Build vs NoPCAP, at same PM2.5 mass
      Persist         = c(0,  0, 1),   # Persist vs NoPCAP, at same PM2.5 mass
      Persist_v_Build = c(0, -1, 1)    # Persist vs Build, at same PM2.5 mass
    ),
    slope_rows = "PM",
    tests      = list(phase_any = c(B, P)),   # any difference across phases (2 df)
    bscov      = "diag"           # 3x3 unstructured = 6 params
  ),
  
  interaction_phase = list(
    rhs        = "pm01c * pcap_phase_3l",
    coef_names = c("pm01c", B, P, paste0("pm01c:", B), paste0("pm01c:", P)),
    contrasts  = rbind(
      slope_NoPCAP       = c(1, 0, 0,  0, 0),   # PM slope, NoPCAP days
      slope_Build        = c(1, 0, 0,  1, 0),   # PM slope, Build days
      slope_Persist      = c(1, 0, 0,  0, 1),   # PM slope, Persist days
      int_Build          = c(0, 0, 0,  1, 0),   # Build slope - NoPCAP slope
      int_Persist        = c(0, 0, 0,  0, 1),   # Persist slope - NoPCAP slope
      int_Persist_v_Build= c(0, 0, 0, -1, 1),   # Persist slope - Build slope
      shift_Build        = c(0, 1, 0,  0, 0),   # Build vs NoPCAP at PM = pm_center
      shift_Persist      = c(0, 0, 1,  0, 0)    # Persist vs NoPCAP at PM = pm_center
    ),
    slope_rows = c("slope_NoPCAP", "slope_Build", "slope_Persist",
                   "int_Build", "int_Persist", "int_Persist_v_Build"),
    tests      = list(
      slopes_differ = paste0("pm01c:", c(B, P)),  # PM slope differs by phase (2 df)
      shifts_differ = c(B, P)                     # phase contrast at pm_center (2 df)
    ),
    bscov      = "diag"           # 5x5 unstructured = 15 params; far too many
  )
)

# ---- Stage 1: one city x one outcome x one model ----------------------------

fit_city <- function(d, outcome, spec) {
  f <- as.formula(paste(outcome, "~", spec$rhs, "+", covariates))
  fit <- gam(f, family = quasipoisson(), data = d, method = "REML")
  
  cn <- spec$coef_names
  missing_terms <- setdiff(cn, names(coef(fit)))
  if (length(missing_terms) > 0) {
    stop("Coefficients not found: ", paste(missing_terms, collapse = ", "),
         " (a phase level with no days in this city?). Model has: ",
         paste(head(names(coef(fit)), 8), collapse = ", "))
  }
  if (anyNA(coef(fit)[cn])) {
    stop("NA coefficient(s): ", paste(cn[is.na(coef(fit)[cn])], collapse = ", "),
         " (term not identifiable in this city)")
  }
  
  list(
    b        = coef(fit)[cn],
    V        = vcov(fit, freq = TRUE)[cn, cn],   # includes quasi-Poisson scale
    n_events = sum(d[[outcome]], na.rm = TRUE)
  )
}

safe_fit_city <- safely(fit_city)

# ---- Stage 2: pool one outcome across cities --------------------------------

pool_outcome <- function(outcome, spec, city_list) {
  res  <- map(city_list, \(d) safe_fit_city(d, outcome, spec))
  errs <- keep(map(res, "error"), negate(is.null))
  if (length(errs) > 0) {
    msgs <- map_chr(errs, conditionMessage)
    message(outcome, ": dropped cities ->\n",
            paste0("  ", names(msgs), ": ", msgs, collapse = "\n"))
  }
  fits <- compact(map(res, "result"))
  if (length(fits) < 2) {
    stop(outcome, ": fewer than 2 cities fit successfully; see messages above.")
  }
  
  Y <- do.call(rbind, map(fits, "b"))
  S <- do.call(rbind, map(fits, \(x) x$V[lower.tri(x$V, diag = TRUE)]))
  rownames(Y) <- rownames(S) <- names(fits)
  
  # Within-city correlation of the PM2.5 coefficient with each other coefficient
  within_cor <- imap_dfr(fits, function(x, city) {
    R <- cov2cor(x$V)
    tibble(city = city, coef = colnames(R)[-1], cor_with_pm = R[1, -1])
  })
  
  # Fallback chain: REML with spec$bscov -> REML diag -> method of moments.
  # mixmeta only *warns* on non-convergence, so check $converged as well as
  # catching errors. REML can also fail outright ("singular matrix") when all
  # between-city variances are ~0 (common for sparse outcomes); the
  # non-iterative method-of-moments estimator always returns an answer.
  fit_mm <- function(method, bscov = "diag") {
    tryCatch(
      # default maxiter = 100 is too low when a between-city variance is near 0
      suppressWarnings(mixmeta(Y ~ 1, S = S, method = method, bscov = bscov,
                               control = list(maxiter = 1000))),
      error = function(e) e
    )
  }
  ok <- \(m) inherits(m, "mixmeta") && (m$method == "mm" || isTRUE(m$converged))
  why <- \(m) if (inherits(m, "error")) conditionMessage(m) else "did not converge"
  
  bscov_used <- spec$bscov
  method_used <- "reml"
  pool <- fit_mm("reml", spec$bscov)
  if (!ok(pool) && spec$bscov != "diag") {
    message(outcome, ": REML bscov = '", spec$bscov, "' ", why(pool), "; trying diag")
    bscov_used <- "diag"
    pool <- fit_mm("reml", "diag")
  }
  if (!ok(pool)) {
    message(outcome, ": REML diag ", why(pool), "; using method of moments")
    bscov_used  <- "unstr (mm)"
    method_used <- "mm"
    pool <- tryCatch(mixmeta(Y ~ 1, S = S, method = "mm"),   # mm takes no bscov
                     error = function(e) e)
  }
  if (!ok(pool)) stop(outcome, ": pooling failed (", why(pool), ")")
  
  list(outcome = outcome, pool = pool, Y = Y, S = S,
       within_cor = within_cor, n_cities = nrow(Y),
       n_events = map_dbl(fits, "n_events"),
       bscov = bscov_used, method = method_used)
}

# ---- Contrasts -> tidy results ----------------------------------------------

tidy_pool <- function(p, spec) {
  L   <- spec$contrasts
  b   <- coef(p$pool)
  V   <- vcov(p$pool)
  Psi <- p$pool$Psi
  s   <- summary(p$pool)
  
  scale <- ifelse(rownames(L) %in% spec$slope_rows, pm_increment, 1)
  
  est  <- drop(L %*% b)
  se   <- sqrt(diag(L %*% V %*% t(L)))
  tau2 <- diag(L %*% Psi %*% t(L))
  
  # summary(mixmeta) returns Q, df, p and I^2 as vectors named
  # c(".all", <coef 1>, <coef 2>, ...): element 1 is the multivariate test,
  # element j + 1 is the univariate test for coefficient j.
  # Per-term stats are only defined for contrasts that pick out one
  # coefficient; derived combinations (e.g. the PCAP-day slope) get NA.
  coef_idx <- apply(L, 1, \(r) if (sum(r != 0) == 1 && sum(r) == 1) which(r == 1) else NA_integer_)
  per_term <- \(x) ifelse(is.na(coef_idx), NA_real_, unname(x[coef_idx + 1]))
  
  i2     <- per_term(s$i2stat)
  Q_term <- per_term(s$qstat$Q)
  Q_p    <- per_term(s$qstat$pvalue)
  
  wald_all <- as.numeric(t(b) %*% solve(V) %*% b)
  
  tibble(
    outcome  = p$outcome,
    term     = rownames(L),
    est      = est * scale,
    se       = se * scale,
    tau2     = tau2 * scale^2,
    I2       = i2,
    p_value  = 2 * pnorm(-abs(est / se)),
    Q        = Q_term,                     # univariate Cochran Q, df = k - 1
    Q_p      = Q_p,
    Q_all    = unname(s$qstat$Q[1]),       # multivariate Q, df = (k - 1) x n_coef
    Q_all_p  = unname(s$qstat$pvalue[1]),
    joint_p  = pchisq(wald_all, df = length(b), lower.tail = FALSE),
    n_cities = p$n_cities,
    method   = p$method,                   # "reml", or "mm" if REML failed
    bscov    = p$bscov
  ) |>
    mutate(
      RR     = exp(est),
      lower  = exp(est - 1.96 * se),
      upper  = exp(est + 1.96 * se),
      RR_fmt = sprintf("%.3f (%.3f, %.3f)", RR, lower, upper)
    )
}

# ---- Joint Wald tests on subsets of pooled coefficients ---------------------

tidy_tests <- function(p, spec) {
  if (length(spec$tests) == 0) return(tibble())
  b <- coef(p$pool)
  V <- vcov(p$pool)
  names(b) <- colnames(V) <- rownames(V) <- spec$coef_names
  
  imap_dfr(spec$tests, function(cn, test) {
    chi2 <- as.numeric(t(b[cn]) %*% solve(V[cn, cn]) %*% b[cn])
    tibble(outcome = p$outcome, test = test, chi2 = chi2, df = length(cn),
           p_value = pchisq(chi2, df = length(cn), lower.tail = FALSE))
  })
}

# ---- City-level BLUPs on the contrast scale (for forest plots only) ---------

tidy_blups <- function(p, spec) {
  L     <- spec$contrasts
  scale <- ifelse(rownames(L) %in% spec$slope_rows, pm_increment, 1)
  bl    <- mixmeta::blup(p$pool, vcov = TRUE)   # explicit: metafor also has blup()
  
  map2_dfr(bl, rownames(p$Y), function(x, city) {
    tibble(
      outcome = p$outcome,
      city    = city,
      term    = rownames(L),
      est     = drop(L %*% x$blup) * scale,
      se      = sqrt(diag(L %*% x$vcov %*% t(L))) * scale
    )
  })
}

# ---- Run a model across all outcomes ----------------------------------------

run_model <- function(model, data, label = "full") {
  spec      <- model_specs[[model]]
  d         <- prep_data(data)
  city_list <- split(d, d$city)
  
  message("== ", label, " / ", model, " ==")
  pooled <- map(outcomes, possibly(\(o) pool_outcome(o, spec, city_list),
                                   otherwise = NULL, quiet = FALSE)) |>
    set_names(outcomes) |>
    compact()   # an outcome that can't be pooled in this dataset is skipped, not fatal
  
  tag <- \(x) if (nrow(x) == 0) x else mutate(x, dataset = label, model = model, .before = 1)
  
  list(
    pooled     = pooled,
    results    = map_dfr(pooled, \(p) tidy_pool(p, spec))  |> tag(),
    tests      = map_dfr(pooled, \(p) tidy_tests(p, spec)) |> tag(),
    blups      = map_dfr(pooled, \(p) tidy_blups(p, spec)) |> tag(),
    within_cor = imap_dfr(pooled, \(p, o) mutate(p$within_cor, outcome = o, .before = 1)) |> tag()
  )
}



# Run Models

dat<-allcity_test3

main_out     <- run_model("main",              dat)
int_out      <- run_model("interaction",       dat)















