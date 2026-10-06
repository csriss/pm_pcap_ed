# Code used to conduct random effects meta analysis of city-specific estimates
# from models with and without PM2.5 and PCAP interactions where PM2.5
# Concentration response is modeled linearly


#Conduct Meta Analysis
# ---- Forest plots ------------------------------------------------------------
# City rows = observed Stage 1 estimates (as in an rma forest plot), put on the
# contrast scale using each city's full vcov. Pooled row = multivariate estimate.

vech2mat <- function(v, k) {
  m <- matrix(0, k, k)
  m[lower.tri(m, diag = TRUE)] <- v
  m + t(m) - diag(diag(m))
}

# Observed city-level estimates for one contrast (e.g. "PCAP", "noPCAP")
city_contrast <- function(p, spec, term) {
  L     <- spec$contrasts[term, ]
  k     <- length(L)
  scale <- if (term %in% spec$slope_rows) pm_increment else 1
  
  tibble(
    city = rownames(p$Y),
    est  = drop(p$Y %*% L) * scale,
    se   = map_dbl(seq_len(nrow(p$S)),
                   \(i) sqrt(drop(t(L) %*% vech2mat(p$S[i, ], k) %*% L))) * scale
  )
}


# Multivariate Meta Analysis: metafor-style forest plot for one model x outcome x term
forest_mv <- function(out, model, outcome, term, ...) {
  spec   <- model_specs[[model]]
  p      <- out$pooled[[outcome]]
  cd     <- city_contrast(p, spec, term)
  pooled <- out$results |> filter(outcome == !!outcome, term == !!term)
  k      <- nrow(cd)
  labs<-out$results |> filter(outcome == !!outcome, term == !!term) %>% select(outcome_labs)
  forest(
    x       = cd$est,
    sei     = cd$se,
    slab    = cd$city,
    atransf = exp,
    main = labs,
    digits=3,
    refline = 0,
    ylim    = c(-1.5, k + 3),
    header  = c("City", "RR [95% CI]"),
    xlab = "",
    cex = 1.2,
    xlim=c(-.35,.45),
    # xlab    = if (term %in% spec$slope_rows)
    #   paste0("RR per ", pm_increment, " \u00b5g/m\u00b3") else "RR",
    ...
  )
  addpoly(
    x       = pooled$est,
    sei     = pooled$se,
    rows    = -1,
    atransf = exp,
    #mlab    = sprintf("Pooled (multivariate RE), \u03c4\u00b2 = %.4f", pooled$tau2)
    mlab    ="Pooled RR"
  )
  abline(h = 0)
}
















