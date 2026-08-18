### Quick script to look at data during RTO and calculate 3-day means

# Load packages------------------
library(here)
library(CDECRetrieve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(readr)
library(sharpshootR)
library(scales)
library(dbscan)
library(zoo)
library(janitor)


# Download Freeport turbidity and also flow

start_date <- as.POSIXct("2025-10-01 00:00:00", tz = "UTC")
end_date   <- as.POSIXct("2026-01-07 23:45:00", tz = "UTC")


turb_data <- CDECquery(id = "FPT", sensor = 221, interval = "D", start = start_date, end = end_date)

q_data <- CDECquery(id = "FPT", sensor = 20, interval = "D", start = start_date, end = end_date)
                  

turb_data <- turb_data %>% 
  mutate(three_d_mean = rollmean(value, k=3, fill=NA, align='right'))

q_data <- q_data %>% 
  mutate(three_d_mean = rollmean(value, k=3, fill=NA, align='right'))
