#Create Multicity Exposure Dataset 

# Libraries ----
library(readxl)
library(tidyverse)
library(weathermetrics)


#Load PM Data
##Modesto

modesto_pm<-read_xlsx("Raw/PM Update Nov 2024/Modesto.xlsx",
                      sheet="Final_PM2.5") %>% 
  mutate(city="Modesto")


##Bakersfield
bakersfield_pm<-read_xlsx("Raw/PM Update Nov 2024/Bakersfield.xlsx",
                          sheet="Final_PM2.5") %>% 
  mutate(city="Bakersfield")

##Fresno
fresno_pm<-read_xlsx("Raw/PM Update Nov 2024/Fresno.xlsx",
                     sheet="Final_PM2.5") %>% 
  rename("SiteCode"="SiteCode1",
         "SiteName"="SiteName1") %>% 
  mutate(city="Fresno")

##Visalia
visalia_pm<-read_xlsx("Raw/PM Update Nov 2024/Visalia.xlsx",
                      sheet="Final_PM2.5") %>% 
  mutate(city="Visalia")

##Sacramento
sacramento_pm<-read_xlsx("Raw/PM Update Nov 2024/Sacramento.xlsx",
                         sheet="Final_PM2.5") %>% 
  mutate(city="Sacramento")

##Reno
reno_pm<-read_xlsx("Raw/PM Update Nov 2024/Reno.xlsx", sheet = "Final_PM2.5") %>% 
  rename("SiteCode"="SiteCode1",
         "SiteName"="SiteName1")%>% 
  mutate(city="Reno")

##Las Vegas
lv_pm<-read_xlsx("Raw/PM Update Nov 2024/Las Vegas.xlsx",sheet = "Final_PM2.5") %>%
  rename("SiteCode"="SiteCode1",
         "SiteName"="SiteName1")%>% 
  mutate(city="Las Vegas")


##Salt Lake City
slc_pm<-read_xlsx("Raw/PM Update Nov 2024/Hawthorne.xlsx",sheet = "Final_PM2.5") %>% 
  mutate(city="SLC")

##Provo
provo_pm<-read_xlsx("Raw/PM Update Nov 2024/Lindon.xlsx",sheet = "Final_PM2.5") %>% 
  mutate(city="Provo")


#Combine
combined_allcity_pm<-rbind(modesto_pm,bakersfield_pm,fresno_pm,sacramento_pm,visalia_pm,
                           reno_pm,lv_pm,slc_pm,provo_pm)


#Automatically Fill Missing Dates
combined_allcity_pm<-combined_allcity_pm %>% 
  # Ensure Date column is proper Date class
  mutate(Date = as.Date(Date)) %>% 
  group_by(city) %>% 
  # Fill in missing dates from min to max date for each city
  complete(Date = seq.Date(min(Date), max(Date), by = "day")) %>% 
  arrange(Date, .by_group = TRUE) %>% 
  ungroup()



#Load Met Data
##California
bakersfield_weather<-read.csv("Raw/epiweather raw/epiweather_updated/weather-Bakersfield.csv")
modesto_weather<-read.csv("Raw/epiweather raw/epiweather_updated/weather-Modesto.csv")
fresno_weather<-read.csv("Raw/epiweather raw/epiweather_updated/weather-Fresno.csv")
sacramento_weather<-read.csv("Raw/epiweather raw/epiweather_updated/weather-Sacramento.csv")
visalia_weather<-read.csv("Raw/epiweather raw/epiweather_updated/weather-Visalia.csv")


##Nevada
lasvegas_weather<-read.csv("Raw/epiweather raw/epiweather_updated/weather-LasVegas.csv")
reno_weather<-read.csv("Raw/epiweather raw/epiweather_updated/weather-Reno.csv")

##Utah
slc_weather<-read.csv("Raw/epiweather raw/epiweather/weather-SLC.csv")%>% 
  mutate(date=ymd(date))

provo_weather<-read.csv("Raw/epiweather raw/epiweather_updated/weather-Lindon.csv") %>% mutate(city="Provo")


#All City weather Combined
allcity_weather<-rbind(bakersfield_weather,modesto_weather,fresno_weather,sacramento_weather,
                     visalia_weather,lasvegas_weather,reno_weather,slc_weather,
                     provo_weather,denver_weather) %>% 
  mutate(date=ymd(date)) %>% 
  select(c("date", "city", "max_dewpt_temp", "mean_dewpt_temp", "min_dewpt_temp", 
           "max_temp", "mean_temp", "min_temp"))


#Automatically Fill Missing Dates
allcity_weather <- allcity_weather %>% 
  # Ensure Date column is proper Date class
  mutate(Date = as.Date(Date)) %>% 
  group_by(city) %>% 
  # Fill in missing dates from min to max date for each city
  complete(Date = seq.Date(min(Date), max(Date), by = "day")) %>% 
  arrange(Date, .by_group = TRUE) %>% 
  ungroup()


#Calculate Relative Humidity
allcity_weather$RH <- dewpoint.to.humidity(t = allcity_weather$mean_temp,
                                           dp = allcity_weather$mean_dewpt_temp,
                                           temperature.metric = "fahrenheit")



#Join PM and Weather Data and Create Lags

allcity_pm_met<-left_join(allcity_weather,combined_allcity_pm,by=c("Date","city")) %>% 
  select(-c(SiteCode,SiteName))

#Set Max Lag
maxlag <- 6

# Create a named list of lag functions: lag1, lag2, ..., lag6
lag_fns <- setNames(
  lapply(1:maxlag, function(i) function(x) lag(x, n = i)),
  paste0("lag", 1:maxlag)
)

# Create lag variables
allcity_pm_met <- allcity_pm_met %>% 
  group_by(city) %>% 
  # 1. Guarantee chronological ordering
  arrange(Date, .by_group = TRUE) %>% 
  # 2. Create all 24 lags in 1 single pass
  mutate(
    across(
      .cols = c(PM2.5 = PM2.5, mean_temp = mean_temp, RH = RH, mean_dewpt_temp = mean_dewpt_temp),
      .fns = lag_fns,
      .names = "{.col}_{.fn}"
    )
  ) %>% 
  # 3. Clean up grouping
  ungroup()




#Load CAP Classification Dataset
cap_vhd_new<-read.csv("Raw/CAP_VHD_Health_final_Sep2025(1).csv")

##Create event length variables

cap_vhd_new <- cap_vhd_new %>% 
  # 1. Clean City Names
  mutate(
    city = case_when(
      Station == "BAK" ~ "Bakersfield",
      Station == "DNR" ~ "Denver",
      Station == "REV" ~ "Reno",
      Station == "BOI" ~ "Boise",
      Station == "VEF" ~ "Las Vegas",
      Station == "OGD" ~ "Ogden",
      Station == "PVO" ~ "Provo",
      Station == "SAC" ~ "Sacramento",
      Station == "FNO" ~ "Fresno",
      Station == "VIS" ~ "Visalia",
      Station == "MOD" ~ "Modesto",
      Station == "MFR" ~ "Medford",
      Station == "SLC" ~ "SLC",
      TRUE ~ NA_character_
    )
  ) %>% 
  # Filter out unwanted cities and unmapped stations
  filter(!city %in% c("Boise", "Medford") & !is.na(city)) %>% 
  
  # 2. Group & Order by Date
  group_by(city) %>% 
  arrange(Date, .by_group = TRUE) %>% 
  
  # 3. Handle 3-State Logic (TRUE, FALSE, Maybe) & Binary Flag
  mutate(
    new_cap = case_when(
      ERA.Adj.CAP == "TRUE" ~ TRUE,
      ERA.Adj.CAP == "Maybe" & (
        lag(ERA.Adj.CAP, default = "") == "TRUE" | 
          lead(ERA.Adj.CAP, default = "") == "TRUE"
      ) ~ TRUE,
      TRUE ~ FALSE
    ),
    cap_code = as.integer(new_cap)
  ) %>% 
  
  # 4. Identify Event Starts & Assign Group IDs
  mutate(
    # Start event if cap_code is 1 and previous row was not 1 (safely defaults border NA to 0)
    event_start = if_else(
      cap_code == 1 & lag(cap_code, default = 0) == 0, 1, 0
    ),
    event_group = cumsum(event_start),
    event_group = if_else(cap_code == 1, event_group, NA_real_)
  ) %>% 
  
  # 5. Calculate Event Length and Duration Metrics
  group_by(city, event_group) %>% 
  mutate(
    event_length = if_else(!is.na(event_group), n(), NA_integer_),
    days_into_event = if_else(!is.na(event_group), row_number(), NA_integer_)
  ) %>% 
  
  # 6. Final Cleanup & Date Filtering
  ungroup() %>% 
  filter(Date <= "2020-02-15")


# Create PCAP Indicator Variable
cap_vhd_new <- cap_vhd_new %>% 
  group_by(city) %>% 
  mutate(
    pcap01 = as.integer(
      new_cap & lag(new_cap, default = FALSE)
    )
  ) %>% ungroup()


#Join Updated CAP data with pm 
pm_met_updated_cap<-left_join(cap_vhd_new,allcity_pm_met,by=c("Date","city"))


#Filter Data to match start dat of health data
pm_met_cap_filtered <- pm_met_updated_cap %>% 
  mutate(Date = as.Date(Date)) %>% 
  filter(
    (city == "Bakersfield"    & Date >= as.Date("2005-01-01")) |
      (city == "Fresno"         & Date >= as.Date("2005-01-01")) |
      (city == "Las Vegas"      & Date >= as.Date("2009-01-01")) |
      (city == "Modesto"        & Date >= as.Date("2010-11-15")) |
      (city == "Provo"          & Date >= as.Date("2001-01-01")) |
      (city == "Reno"           & Date >= as.Date("2010-12-15")) |
      (city == "Sacramento"     & Date >= as.Date("2010-11-15")) |
      (city == "Salt Lake City" & Date >= as.Date("2001-01-01")) |
      (city == "Visalia"        & Date >= as.Date("2005-01-01"))
  )


saveRDS(pm_met_cap_filtered,"pm_met_updated_cap_filtered_Aug2026.rds")








