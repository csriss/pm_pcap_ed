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























