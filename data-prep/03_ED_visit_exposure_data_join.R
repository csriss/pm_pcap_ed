#Join ED visit daily count datasets and Multicity Exposure Dataset 

#Libraries
library(tidyverse)

#Utah

ut_winter_daily<-readRDS("ut_winter_daily_new_15apr2026.rds")


#California

ca_winter_daily<-readRDS("ca_winter_daily_new_15apr2026.rds")

ca_ut_winter_daily<-rbind(ca_winter_daily,ut_winter_daily)



# Nevada
nv_winter_daily<-readRDS("nv_winter_daily_new_15apr2026.rds")


#Combined ED Data from Each STate

allstates_winter_daily<-rbind(ca_ut_winter_daily,nv_winter_daily)


#Create year variable
allstates_winter_daily$year<-year(allstates_winter_daily$Date)

#Group cities by state
allstates_winter_daily<-allstates_winter_daily %>% 
  mutate(state=case_when(
    city %in% c("Visalia","Fresno","Bakersfield","Modesto","Sacramento")~"CA",
    city %in% c("Salt Lake City","Provo")~"UT",
    city %in% c("Reno","Las Vegas")~"NV"))



#Create day of winter (t) variable as 1-93 across study period days

t_seq<-data.frame(Date=seq.Date(from=as.Date("2000-11-15"),to=as.Date("2001-02-15"),by="days"))



t_seq<-t_seq %>% 
  mutate(monthday=substr(Date,6,10),
         t=1:93) %>% 
  select(-c(Date))

#Create Monthday variable  to match t_seq
allstates_winter_daily<-allstates_winter_daily %>% 
  mutate(monthday=substr(Date,6,10))


#Join with ED visit Data
allstates_winter_daily2<-left_join(allstates_winter_daily,t_seq,by="monthday")



#Create Season Group and Winter Year Variables
allstates_winter_daily2<-allstates_winter_daily2 %>% 
  arrange(Date) %>% 
  group_by(Date,city) %>% 
  mutate(dow=wday(Date),
         doy=yday(Date),
         month=month(Date),
         year=year(Date)) %>% 
  ungroup() %>% 
  mutate(
    # Define winter year (Nov-Dec belong to current year's winter, Jan-Feb to previous year's winter)
    winter_year = case_when(
      month %in% c(11, 12) ~ year,
      month %in% c(1, 2) ~ year - 1,
      TRUE ~ NA_real_
    ),
    # Create season group
    season_group = case_when(
      month %in% c(11, 12, 1, 2) & winter_year >= 2000 & winter_year <= 2020 ~ 
        paste0("Winter_", sprintf("%02d", winter_year %% 100), "_", sprintf("%02d", (winter_year + 1) %% 100)),
      TRUE ~ NA_character_),
    city_season_group=paste0(city,season_group)
  ) 



#Load Exposure Dataset
pm_met_updated_cap<-readRDS("pm_met_updated_cap_filtered_Aug2026.rds")


#Join ED visit and Exposure Dataset
allcity_pm_cap_ed_test<-left_join(allstates_winter_daily2,
                                  pm_met_updated_cap,
                                  by=c("Date","city"))

#Save
saveRDS(allcity_pm_cap_ed_test, "allcity_pm_cap_ed_test_26aug2026.rds")


#Create moving average variables of PM2.5, temperature, humidity, and dewpoint
allcity_pm_cap_ed_test2<-allcity_pm_cap_ed_test %>% 
  mutate( 
    temp03 = rowMeans(
      across(c(mean_temp, mean_temp_lag1,mean_temp_lag2,mean_temp_lag3), ~ .x),
      na.rm = TRUE
    ),
    rh03=rowMeans(
      across(c(RH, RH_lag1,RH_lag2,RH_lag3), ~ .x),
      na.rm = TRUE
    ),
    dewpt03=rowMeans(
      across(c(mean_dewpt_temp, mean_dewpt_temp_lag1,mean_dewpt_temp_lag2,mean_dewpt_temp_lag3), ~ .x),
      na.rm = TRUE
    ),
    pm01=rowMeans(
      across(c(PM2.5, PM2.5_lag1), ~ .x),
      na.rm = TRUE
    )) %>% 
  relocate(c(pm01,temp03,rh03,dewpt03),.after = PM2.5)


#Create Holiday Indicator Variable
# Load necessary libraries
library(dplyr)
library(lubridate)

# Function to create a dataset with holiday indicators
create_holiday_indicator <- function(start_date = "2013-01-01", end_date = "2020-12-31") {
  # Create a sequence of dates
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "day")
  
  # Create a data frame with dates
  date_df <- data.frame(date = dates, stringsAsFactors = FALSE)
  
  # Extract year for holiday calculations
  date_df$year <- year(date_df$date)
  date_df$month <- month(date_df$date)
  date_df$day <- day(date_df$date)
  
  # Initialize holiday indicator
  date_df$is_holiday <- 0
  
  # Function to calculate Thanksgiving (4th Thursday in November)
  get_thanksgiving <- function(year) {
    # Get the first day of November
    nov_first <- as.Date(paste0(year, "-11-01"))
    # Get the day of the week (1-7, where 1 is Sunday)
    first_day_of_week <- wday(nov_first)
    # Calculate days until first Thursday (if Sunday=1, then Thursday=5)
    days_until_first_thurs <- (5 - first_day_of_week) %% 7
    # First Thursday
    first_thurs <- nov_first + days_until_first_thurs
    # Fourth Thursday is 3 weeks later
    fourth_thurs <- first_thurs + 21
    return(fourth_thurs)
  }
  
  # Function to calculate MLK Day (3rd Monday in January)
  get_mlk_day <- function(year) {
    # Get the first day of January
    jan_first <- as.Date(paste0(year, "-01-01"))
    # Get the day of the week (1-7, where 1 is Sunday)
    first_day_of_week <- wday(jan_first)
    # Calculate days until first Monday (if Sunday=1, then Monday=2)
    days_until_first_mon <- (2 - first_day_of_week) %% 7
    # First Monday
    first_mon <- jan_first + days_until_first_mon
    # Third Monday is 2 weeks later
    third_mon <- first_mon + 14
    return(third_mon)
  }
  
  # Function to calculate Presidents' Day (3rd Monday in February)
  get_presidents_day <- function(year) {
    # Get the first day of February
    feb_first <- as.Date(paste0(year, "-02-01"))
    # Get the day of the week (1-7, where 1 is Sunday)
    first_day_of_week <- wday(feb_first)
    # Calculate days until first Monday
    days_until_first_mon <- (2 - first_day_of_week) %% 7
    # First Monday
    first_mon <- feb_first + days_until_first_mon
    # Third Monday is 2 weeks later
    third_mon <- first_mon + 14
    return(third_mon)
  }
  
  # Loop through each year and mark holidays
  for (yr in unique(date_df$year)) {
    # Fixed date holidays
    date_df$is_holiday[date_df$month == 12 & date_df$day == 24 & date_df$year == yr] <- 1  # Christmas Eve
    date_df$is_holiday[date_df$month == 12 & date_df$day == 25 & date_df$year == yr] <- 1  # Christmas Day
    date_df$is_holiday[date_df$month == 12 & date_df$day == 31 & date_df$year == yr] <- 1  # New Year's Eve
    date_df$is_holiday[date_df$month == 1 & date_df$day == 1 & date_df$year == yr] <- 1   # New Year's Day
    date_df$is_holiday[date_df$month == 2 & date_df$day == 14 & date_df$year == yr] <- 1  # Valentine's Day
    
    # Variable date holidays
    thanksgiving <- get_thanksgiving(yr)
    black_friday <- thanksgiving + 1
    mlk_day <- get_mlk_day(yr)
    presidents_day <- get_presidents_day(yr)
    
    date_df$is_holiday[date_df$date == thanksgiving] <- 1    # Thanksgiving
    date_df$is_holiday[date_df$date == black_friday] <- 1    # Black Friday
    date_df$is_holiday[date_df$date == mlk_day] <- 1         # MLK Day
    
    # Only mark Presidents' Day if it falls on or before Feb 15
    if (presidents_day <= as.Date(paste0(yr, "-02-15"))) {
      date_df$is_holiday[date_df$date == presidents_day] <- 1
    }
  }
  
  # Filter for just the dates between 11/15 and 2/15 each year
  filtered_df <- date_df %>%
    mutate(month_day = paste0(month, "-", day)) %>%
    filter((month == 11 & day >= 15) | 
             month == 12 | 
             month == 1 | 
             (month == 2 & day <= 15))
  
  # Return both the full dataset and filtered dataset
  return(list(
    full_date_df = date_df %>% select(date, is_holiday),
    winter_holidays = filtered_df %>% select(date, is_holiday)
  ))
}

# Example usage:
holiday_data <- create_holiday_indicator("2001-01-01", "2020-12-31")

# Get the full date range with holiday indicators
full_holidays <- holiday_data$full_date_df


# Join Analysis data with holiday indicator dataset
allcity_pm_cap_ed_test2<-left_join(allcity_pm_cap_ed_test2,full_holidays,by=c("Date"="date"))


saveRDS(allcity_pm_cap_ed_test2,"allcity_pm_cap_ed_26aug2026.rds")










