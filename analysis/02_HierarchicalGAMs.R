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


##Get C-R curves for all outcomes
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













