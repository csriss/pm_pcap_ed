##########################################
## Sensitivity Analysis: this code replicates the sensitivity analysis
# comparing pooled associations before and after the transition 
# from ICD 9 to ICD 10
########################################

## ICD Stratify



# Create ICD Version Indicator
allcity_test3<-allcity_test3 %>% 
  mutate(icd_version_flag=ifelse(Date>"2015-11-14","10","9"))

nested_icd9<-allcity_test3 %>% 
  filter(icd_version_flag=="9") %>% 
  nest(data = -city)


### ICD 9 PM No Interaction
pm_gam_icd9 <- map_dfr(outcomes, function(current_outcome) {
  nested_icd9 %>%
    mutate(
      models  = map(data, ~city_pm_gam_adj_pcap(.x, current_outcome)),
      tidied  = map(models, ~tidy(.x, parametric = TRUE)),  # parametric terms only
      outcome = current_outcome
    ) %>%
    unnest(tidied) %>%
    filter(term %in%  c("pm01","pcap011"))}) %>% select(-c(models,data))



# Run meta analysis
meta_pm_gam_icd9<-meta_pm_gam_fun(pm_gam_icd9)



## ICD 10 PM No Interaction

nested_icd10<-allcity_test3 %>% 
  filter(icd_version_flag=="10") %>% 
  nest(data = -city)

pm_gam_icd10 <- map_dfr(outcomes, function(current_outcome) {
  nested_icd10 %>%
    mutate(
      models  = map(data, ~city_pm_gam_adj_pcap(.x, current_outcome)),
      tidied  = map(models, ~tidy(.x, parametric = TRUE)),  # parametric terms only
      outcome = current_outcome
    ) %>%
    unnest(tidied) %>%
    filter(term %in%  c("pm01","pcap011"))}) %>% select(-c(models,data))



# Run Meta Analysis
meta_pm_gam_icd10<-meta_pm_gam_fun(pm_gam_icd10)


# Calculate percent change in risk estimates across periods

meta_pm_gam_icd9$icd_version<-"9"
meta_pm_gam_icd10$icd_version<-"10"


t9<-meta_pm_gam_icd9 %>% 
  filter(term=="pm01") %>% 
  mutate(est10unit=10*pooled_est,
         se10unit=10*pooled_se) %>% 
  select(outcome,term,est10unit,se10unit) %>% 
  mutate(icd_version="9")

t10<-meta_pm_gam_icd10 %>% 
  filter(term=="pm01") %>% 
  mutate(est10unit=10*pooled_est,
         se10unit=10*pooled_se) %>% 
  select(outcome,term,est10unit,se10unit) %>% 
  mutate(icd_version="10")

t910<-rbind(t9,t10)

t910wide<-t910 %>% 
  pivot_wider(names_from=icd_version,values_from = c("est10unit","se10unit"),id_cols=c(outcome,term)) %>% 
  mutate(diff_log_rr = est10unit_10 - est10unit_9,
         se_diff=sqrt(se10unit_9^2 + se10unit_10^2),
         lower_log = diff_log_rr - (1.96 * se_diff),
         upper_log = diff_log_rr + (1.96 * se_diff),
         pct_diff  = (exp(diff_log_rr) - 1) * 100,
         pct_lower = (exp(lower_log) - 1) * 100,
         pct_upper = (exp(upper_log) - 1) * 100
  )



## ICD 9 PM PCAP Interaction Models


pm_pcap_icd9 <- map_dfr(outcomes, function(current_outcome) {
  nested_icd9 %>%
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

# Run Meta Analysis
meta_pm_pcap_icd9 <- map_dfr(outcomes, ~meta_pm_pcap(pm_pcap_icd9, .x))

meta_pm_pcap_icd9<-left_join(meta_pm_pcap_icd9,outcome_df,by="outcome")



## ICD 10 PM PCAP Interaction Models

pm_pcap_icd10 <- map_dfr(outcomes, function(current_outcome) {
  nested_icd10 %>%
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

# Run Meta Analysis
meta_pm_pcap_icd10 <- map_dfr(outcomes, ~meta_pm_pcap(pm_pcap_icd10, .x))

meta_pm_pcap_icd10<-left_join(meta_pm_pcap_icd10,outcome_df,by="outcome")













































