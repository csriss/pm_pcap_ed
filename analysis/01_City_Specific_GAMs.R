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


# Define GAM fitting function with linear PM term and PCAP included in model

# Function with outcome parameter
city_pm_gam_adj_pcap <- function(df, outcome) {
  formula_str <- paste(outcome, 
                       "~pm01 + pcap01 +  s(temp03,bs='cr',k=10) + s(dewpt03,bs='cr') +  s(t,bs='cr',k=5,m=2) + season_group + dow + is_holiday ")
  
  gam(as.formula(formula_str), 
      family = quasipoisson, 
      data = df,method="REML")
}

#Nest data by city
allcity_test3<-as_tibble(allcity_test3)

nested<-allcity_test3 %>% 
  nest(data = -city)


#Run model for all cities

res_pm_gam_adj_pcap <- map_dfr(outcomes, function(current_outcome) {
  nested %>%
    mutate(
      models  = map(data, ~city_pm_gam_adj_pcap(.x, current_outcome)),
      tidied  = map(models, ~tidy(.x, parametric = TRUE)),  # parametric terms only
      outcome = current_outcome
    ) %>%
    unnest(tidied) %>%
    filter(term %in% c("pm01","pcap011"))})# keep only the PM2.5 and pcap coefs

#Remove models and data used to fit models
res_pm_gam_adj_pcap<-res_pm_gam_adj_pcap%>% select(-c(models,data))


#Calculate Rate Ratio estimates
res_pm_gam_adj_pcap<-res_pm_gam_adj_pcap%>% 
  mutate(
    conf.low              = estimate - 1.96 * std.error,
    conf.high             = estimate + 1.96 * std.error,
    RR                    = ifelse(term=="pm01",exp(10 * estimate),exp(estimate)),
    RR_conf.low           = ifelse(term=="pm01",exp(10 * estimate - 1.96 * (std.error * 10)),
                                   exp(estimate - 1.96 * (std.error))),
    RR_conf.high   = ifelse(term=="pm01",exp(10 * estimate + 1.96 * (std.error * 10)),exp(estimate + 1.96 * (std.error)))
  )



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




#Define function for fitting PM*PCAP interaction models

city_pm_pcap_gam <- function(df, outcome) {
  formula_str <- paste(outcome, 
                       "~pm01*pcap01 + s(temp03,bs='cr',k=10) + s(dewpt03,bs='cr') + s(t,bs='cr',k=5,m=2) + season_group + dow + is_holiday ")
  
  gam(as.formula(formula_str), 
      family = quasipoisson, 
      data = df,method="REML")
}



#Run PM*PCAP Interaction Models for all cities

all_results_pm_pcap_gam <- map_dfr(outcomes, function(current_outcome) {
  nested %>%
    mutate(
      models = map(data, ~city_pm_pcap_gam(.x, current_outcome)),
      outcome = current_outcome,
      
      coefs = map(models, ~{
        td  <- tidy(.x, parametric = TRUE)
        vc  <- vcov(.x)
        
        est_main <- td$estimate[td$term == "pm01"]
        est_int  <- td$estimate[td$term == "pm01:pcap011"]
        se_main  <- td$std.error[td$term == "pm01"]
        se_int   <- td$std.error[td$term == "pm01:pcap011"]
        se_pcap  <- sqrt(vc["pm01","pm01"] + 
                           vc["pm01:pcap011","pm01:pcap011"] + 
                           2 * vc["pm01","pm01:pcap011"])
        
        tibble(
          term                 = c("noPCAP", "PCAP","interaction"),
          estimate             = c(est_main, est_main + est_int, est_int),
          se                   = c(se_main,  se_pcap, se_int),
          RR_10unit            = exp(10 * estimate),
          RR_10unit_conf.low   = exp(10 * (estimate - 1.96 * se)),
          RR_10unit_conf.high  = exp(10 * (estimate + 1.96 * se))
        )
      })
    ) %>%
    select(city, outcome, coefs) %>%
    unnest(coefs)
})


#Meta analysis and forest plots of interaction model estimates

for (i in 1:6) {
  
  # Plot 1: No PCAP
  all_results_pm_pcap_gam %>% 
    filter(outcome == outcomes[i], term == "noPCAP") %>% 
    rma(yi = estimate, sei = se, data = ., method = "REML") %>% 
    forest(slab = city, xlab = "RR PM No PCAP",
           atransf = exp_10_unit, digits = 3, header = outcome_labs[i])
  
  # Plot 2: PCAP
  all_results_pm_pcap_gam %>% 
    filter(outcome == outcomes[i], term == "PCAP") %>% 
    rma(yi = estimate, sei = se, data = ., method = "REML") %>% 
    forest(slab = city, xlab = "RR PM PCAP",
           atransf = exp_10_unit, digits = 3, header = outcome_labs[i])
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





