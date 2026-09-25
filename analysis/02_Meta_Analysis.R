# Code used to conduct random effects meta analysis of city-specific estimates
# from models with and without PM2.5 and PCAP interactions where PM2.5
# Concentration response is modeled linearly


#Conduct Meta Analysis

meta_pm_gam_adj_pcap <- res_pm_gam_adj_pcap %>%
  group_by(outcome,term) %>%
  nest() %>%
  mutate(
    meta = map(data, ~rma(
      yi     = .x$estimate,
      sei    = .x$std.error,
      method = "REML",
      slab   = .x$city
    )),
    pooled_est = map_dbl(meta, ~as.numeric(.x$beta)),
    pooled_se  = map_dbl(meta, ~.x$se),
    I2         = map_dbl(meta, ~.x$I2),
    tau2       = map_dbl(meta, ~.x$tau2),
    Q          = map_dbl(meta, ~.x$QE),
    Q_pval     = map_dbl(meta, ~.x$QEp)
  ) %>%
  select(-data, -meta) %>%
  mutate(
    p_pooled            = round(2 * pnorm(-abs(pooled_est / pooled_se)),digits=4)
  ) %>% 
  mutate(
    conf.low              = pooled_est - 1.96 * pooled_se,
    conf.high             = pooled_est + 1.96 * pooled_se,
    RR                    = ifelse(term=="pm01",exp(10 * pooled_est),exp(pooled_est)),
    RR_conf.low           = ifelse(term=="pm01",exp(10 * pooled_est - 1.96 * (pooled_se * 10)),
                                   exp(pooled_est - 1.96 * (pooled_se))),
    RR_conf.high   = ifelse(term=="pm01",exp(10 * pooled_est + 1.96 * (pooled_se * 10)),
                            exp(pooled_est + 1.96 * (pooled_se)))
  )



#Create Forest Plots of PM2.5 Effect Estimates
res_pm_gam_adj_pcap2<-res_pm_gam_adj_pcap %>% 
  mutate(est10unit=10*estimate,
         se10unit=10*std.error) %>% 
  filter(term=="pm01")

par(mfrow=c(2,3))
for (i in 1:6){
  res_pm_gam_adj_pcap2 %>% 
    filter(outcome==outcomes[i]) %>% 
    rma(data=.,yi=est10unit,sei=se10unit,method = "REML") %>% 
    forest(slab=city,xlab = "Per 10 µg/m³ increase in PM2.5",mlab="Pooled RR",
           atransf=exp, digits = 3,header = "City",xlim=c(-.2,.3),
           main=outcome_labs[i],cex = 1.2)
}




# Meta Analysis with estimates for interaction term

meta_pm_pcap <- function(results_df, outcome_name) {
  
  results_df <- results_df %>% filter(outcome == outcome_name)
  
  map_dfr(c("noPCAP", "PCAP", "interaction"), function(trm) {
    dat <- results_df %>% filter(term == trm)
    
    fit <- rma(
      yi     = estimate,
      vi     = se^2,
      data   = dat,
      method = "REML"
    )
    
    tibble(
      outcome            = outcome_name,
      term               = trm,
      pooled_estimate    = fit$beta[1],
      pooled_se          = fit$se,
      ci_lo              = fit$ci.lb,
      ci_hi              = fit$ci.ub,
      p_pooled           = fit$pval,
      I2                 = fit$I2,
      Q                  = fit$QE,
      Q_pval             = fit$QEp,
      RR_10unit          = exp(10 * fit$beta[1]),
      RR_10unit_conf.low = exp(10 * fit$ci.lb),
      RR_10unit_conf.high= exp(10 * fit$ci.ub)
    )
  })
}


meta_results <- map_dfr(outcomes, ~meta_pm_pcap(all_results_pm_pcap_gam, .x))

meta_results<-left_join(meta_results,outcome_df,by="outcome")

meta_results$p_pooled<-round(meta_results$p_pooled,digits=4)

meta_results<-meta_results %>% 
  mutate(rr_fmt=paste(round(RR_10unit,digits=3),"(",round(RR_10unit_conf.low,digits=3),", ",round(RR_10unit_conf.high,digits=3),")"))

















