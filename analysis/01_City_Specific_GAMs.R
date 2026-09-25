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






