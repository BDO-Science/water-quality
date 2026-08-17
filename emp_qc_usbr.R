# WY 2025 Data QC - Reclamation EMP stations
# Code originally created by Catarina Pien
# Modified for 2025 by Lilly McCormick (lmccormick@usbr.gov)
# This code cleans up USBR's EMP data from CDEC data for EMP reporting 
# Lots of very specific identifying of points rather than a general script

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

# Notes -------------------------
# C2 – Sacramento River at Collinsville (CLL): EC (100 - event), WT (25 - hourly)
  # also has bottom EC
# C4 – SJR at San Andreas Landing (SAL): EC (100 - event), WT (25 - hourly)
# C5 – Contra Costa Canal at Pumping Plant #1 (CNT): EC (100 - event), WT (25 - hourly)
# C8 – Old River near Middle River (UNI): EC (100 - event), WT (25 - hourly)
# C14 – Sacramento River at Port Chicago (PCT): EC (sensor 100), WT (25 - event)
  # also has bottom EC
# *C19 – Cache Slough at City of Vallejo Intake (CCS): discontinued 1/1/2015
# D22 – Sacramento River at Emmaton (EMM): EC (100 - event), WT (25 - hourly)
  # Also has bottom EC
# DMC1 – Delta-Mendota Canal at Jones Pumping Plant (DMC): EC (100 - event), WT (25 - hourly)

#####  UPDATES for WY 2025
# Got access to HAR system for QAQC Reclamation EC data.
# Other data (temp) will still need to be downloaded from CDEC

# Download data --------------------------
# removed CCS since no data
stations <- c("CLL", "SAL", "CNT", "UNI", "PCT", "EMM", "DMC")
start_date <- as.POSIXct("2024-10-01 00:00:00", tz = "UTC")
end_date   <- as.POSIXct("2025-09-30 23:45:00", tz = "UTC")
latlon <- data.frame(
  station = stations,
  lat = c(38.073981,38.103301,37.994950,37.822101, 38.056788, 38.084272, 37.781567),
  lon = c(-121.850123,-121.591315,-121.702809, -121.374990, -121.999971, -121.738924, -121.590239)
)


## Water Temperature --------------------------

wt_data <- lapply(stations,
                  function(x) {
                    CDECquery(id = x, sensor = 25, interval = "H", start = start_date,end = end_date)
                  })
wt_raw <- bind_rows(wt_data)%>%
  mutate(datetime = format(as.POSIXct(datetime), "%Y-%m-%d %H:%M:%S"),
         date = date(datetime),
         year = year(date))

saveRDS(wt_raw, "data/data_raw/wt_raw_2025.rds")


## Water Temperature -------------------
wt_raw0 <- readRDS(here("data/data_raw/wt_raw_2025.rds"))%>%
  mutate(datetime = ymd_hms(datetime)) %>%
  filter(date > ymd("2024-09-30")) %>%
  rename(station = station_id)

# Make hourly dataset
hour_time <- tibble(
  datetime = seq(from = start_date,
                 to   = end_date,
                 by   = "1 hour")
)

# Combine with stations
sta_hour <- crossing(stations, hour_time) %>%
  rename(station = stations)

# Join with WT data
wt_raw <- sta_hour %>%
  left_join(wt_raw0, by = c("station", "datetime")) %>%
  mutate(station_d1641 = case_when(station == "CLL" ~ "C2",
                                   station== "SAL" ~ "C4",
                                   station == "CNT" ~ "C5",
                                   station == "UNI" ~ "C8",
                                   station == "PCT" ~ "C14",
                                   station == "CCS" ~ "C19",
                                   station == "EMM" ~ "D22",
                                   station == "DMC" ~ "DMC1")) %>%
  mutate(date = date(datetime),
         month = factor(month(date))) %>%
  select(station,
         station_d1641,
         wt = value,
         datetime, date, month)

### QA/QC --------------
#### Attempt 1 ---------------
# tried a few rules here but too inconsistent in either flagging too many or too few points. 

# all in one dataframe: 
wt_qc <- wt_raw %>%
  # belowzero 
  mutate(range = if_else(wt <=32 | wt>=90, 1L, 0L))%>%
  # duplicated
  group_by(station) %>%
  mutate(
    wt_lag = lag(wt, 1),
    same = case_when(
      is.na(wt) & is.na(wt_lag) ~ 0L,   # both NA: treat as a "change" (breaks the run)
      is.na(wt) | is.na(wt_lag) ~ 0L,   # one NA, one real value: also breaks the run
      wt == wt_lag ~ 1L,
      TRUE ~ 0L
    ),
    issame = cumsum(same == 0L)
  ) %>%
  ungroup() %>%
  group_by(station, issame) %>%
  mutate(flag = sum(same) + 1) %>%
  ungroup() %>%
  mutate(repeated = ifelse(flag > 36, 1L, 0L)) %>% 
  group_by(station) %>%
  mutate(rmean = zoo::rollapply(wt, 36, mean, partial =TRUE),
         rsd = zoo::rollapply(wt, 36, sd, partial = TRUE)) %>%
  ungroup() %>% 
  mutate(outlier = ifelse(wt < rmean-3*rsd | wt > rmean+3*rsd, 1L, 0L)) %>%
  ungroup() %>%
  # spike test
  group_by(station) %>%
  mutate(spike = if_else(wt > 1.5 * lag(wt, 1, default = first(wt)) & 
                           wt > 1.5 * lead(wt, 1, default = first(wt)), 1L, 0L)) %>%
  ungroup() %>%
  # summarize
  mutate(flagged = if_else(range == 1L |  repeated == 1L | outlier == 1L | spike == 1L, "flagged", "not flagged")) 

# check what this looks like
ggplot(wt_qc)+ #%>% filter(!is.na(datetime))) + 
  geom_point(aes(datetime, wt, color = flagged, shape = flagged))+
  facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "wt (F)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom")


### View data
wt_raw$station <- as.factor(wt_raw$station)

ggplot(wt_qc %>% filter(wt < 100))+
  geom_point(aes(datetime, wt))+
  facet_wrap(~station)

cll <- wt_qc %>% 
  filter(station == "CLL") #%>% 
  #filter(wt< 100)

ggplot(cll, aes(datetime, wt))+
  geom_point(aes(datetime, wt, color = flagged, shape = flagged))+
  #facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "wt (F)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom")



cnt <- wt_qc %>% 
  filter(station == "CNT") #%>% 
#filter(wt< 100)

ggplot(cnt, aes(datetime, wt))+
  geom_point(aes(datetime, wt, color = flagged, shape = flagged))+
  #facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "wt (F)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom")


dmc <- wt_qc %>% 
  filter(station == "DMC") #%>% 
#filter(wt< 100)
library(plotly)

ggplot(dmc, aes(datetime, wt))+
  geom_point(aes(datetime, wt, color = flagged, shape = flagged))+
  #facet_wrap(~station, scales = "free_y", ncol = 2) +
  #labs(y = "wt (F)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_y_continuous(name = "wt (F)", n.breaks = 16)+
  scale_color_manual(values = c("magenta", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom")



emm <- wt_qc %>% 
  filter(station == "EMM") %>% 
  filter(wt< 100)

ggplot(emm, aes(datetime, wt))+
  geom_point(aes(datetime, wt, color = flagged, shape = flagged))+
  #facet_wrap(~station, scales = "free_y", ncol = 2) +
  #labs(y = "wt (F)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_y_continuous(name = "wt (F)", n.breaks = 16)+
  scale_color_manual(values = c("magenta", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom")

pct <- wt_qc %>% 
  filter(station == "PCT") %>% 
  filter(wt< 100)

ggplot(pct, aes(datetime, wt))+
  geom_point(aes(datetime, wt, color = flagged, shape = flagged))+
  #facet_wrap(~station, scales = "free_y", ncol = 2) +
  #labs(y = "wt (F)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_y_continuous(name = "wt (F)", n.breaks = 16)+
  scale_color_manual(values = c("magenta", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom")


sal <- wt_qc %>% 
  filter(station == "SAL") #%>% 
  #filter(wt< 100)

ggplot(sal, aes(datetime, wt))+
  geom_point(aes(datetime, wt, color = flagged, shape = flagged))+
  #facet_wrap(~station, scales = "free_y", ncol = 2) +
  #labs(y = "wt (F)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_y_continuous(name = "wt (F)", n.breaks = 16)+
  scale_color_manual(values = c("magenta", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom")

uni <- wt_qc %>% 
  filter(station == "UNI") #%>% 
#filter(wt< 100)

ggplot(uni, aes(datetime, wt))+
  geom_point(aes(datetime, wt, color = flagged, shape = flagged))+
  #facet_wrap(~station, scales = "free_y", ncol = 2) +
  #labs(y = "wt (F)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_y_continuous(name = "wt (F)", n.breaks = 20)+
  scale_color_manual(values = c("magenta", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom")



#### Attempt 2 (this is what I used) -----------
# manual removal of visual outliers
wt_qc[61237:61320, 17] <- "not flagged"


wt_qc2 <- wt_qc %>%
  mutate(vis_flag = case_when(station == "CLL" & datetime >= "2024-11-29 12:00:00 UTC" & datetime <="2024-12-04 10:00:00 UTC" & wt > 100 ~1L,
                              station == "CNT" & date >=ymd("2025-03-04") & date <= ymd("2025-04-09") ~1L,
                              station == "DMC" & date >=ymd("2024-10-01") & date <=ymd("2024-10-06") & wt > 76 ~1L,
                              station == "DMC" & date >=ymd("2024-10-01") & date <=ymd("2024-10-06") & wt <= 73  ~1L,
                              station == "DMC" & date >=ymd("2024-12-25") & date <= ymd("2025-01-02") & wt <= 51 ~1L,
                              station == "DMC" & date >=ymd("2025-01-20") & date <= ymd("2025-02-01") & wt <= 48 ~1L,
                              station == "DMC" & date >=ymd("2025-01-20") & date <= ymd("2025-02-01") & wt >=52 ~1L,
                              station == "DMC" & datetime >= "2025-04-17 16:00:00 UTC" & datetime <= "2025-05-01 09:00:00 UTC" ~1L,
                              station == "DMC" & date >=ymd("2025-08-25") & wt < 70 ~1L,
                              station == "EMM" & wt > 100 ~1L,
                              station == "EMM" & datetime >= "2025-05-15 11:00:00 UTC" & datetime <= "2025-05-16 11:00:00 UTC" & wt <65 ~1L,
                              station == "PCT" & wt > 100 ~1L,
                              station == "PCT" & date >=ymd("2025-09-01") & wt >= 73 ~1L,
                              station == "PCT" & date >=ymd("2025-03-15") & date <= ymd("2025-03-25") & wt >= 58 ~1L,
                              station == "UNI" & datetime == "2024-10-01 16:00:00 UTC" ~1L,
                              station == "UNI" & datetime == "2024-10-03 16:00:00 UTC" ~1L,
                              station == "UNI" & datetime == "2024-10-04 16:00:00 UTC" ~1L,
                              station == "UNI" & datetime == "2024-10-05 16:00:00 UTC" ~1L,
                              station == "UNI" & datetime == "2024-11-02 14:00:00 UTC" ~1L,
                              station == "UNI" & datetime >= "2024-10-01 00:00:00 UTC" & datetime <="2024-10-06 15:00:00 UTC" & wt > 75 ~1L,
                              station == "UNI" & datetime >= "2024-11-22 05:00:00 UTC" & datetime <="2024-12-18 07:00:00 UTC" & wt < 50 ~1L,
                              station == "UNI" & datetime >= "2024-11-22 05:00:00 UTC" & datetime <="2024-12-18 07:00:00 UTC" & wt >55.6 ~1L, 
                              station == "UNI" & datetime >= "2025-01-18 06:00:00 UTC" & datetime <="2025-02-03 05:00:00 UTC" & wt <= 45 ~1L,
                              station == "UNI" & datetime >= "2025-01-18 06:00:00 UTC" & datetime <="2025-02-03 05:00:00 UTC" & wt >= 56.5 ~1L,
                              TRUE~0L)) %>%
  mutate(flagged = if_else(range == 1L |  vis_flag == 1L, "flagged", "not flagged")) 

# check with this looks like
ggplot(wt_qc2 %>% filter(!is.na(datetime))) + 
  geom_point(aes(datetime,wt, color = flagged, shape = flagged), size = 0.75)+
  facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "wt (F)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom")

# more useful to get rid of high vals first
ggplot(wt_qc2 %>% filter(range == 0L)) + 
  geom_point(aes(datetime,wt, color = flagged, shape = flagged), size = 0.75)+
  facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "wt (F)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom")

#### Cleaned ----------------
wt_check = wt_qc %>%
  filter(flagged == "flagged")

wt_clean <- wt_qc2 %>%
  filter(flagged == "not flagged")

ggplot(wt_clean %>% filter(!is.na(datetime))) + 
  geom_point(aes(datetime,wt))+
  facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "wt (uS/cm)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  theme_bw() + 
  theme(legend.position = "bottom")

#### Individual stations-------------------
(cnt_wt <- ggplot(wt_qc2 %>% filter(station == "CNT"))+
  geom_point(aes(datetime, wt, color = flagged), size = 1) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))) +
  theme_bw()
plotly::ggplotly(cnt_wt)

pct_wt <- ggplot(wt_qc2 %>% filter(station == "PCT"))+
  geom_point(aes(datetime, wt, color = flagged), size = 1) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))
plotly::ggplotly(pct_wt)

emm_wt <- ggplot(wt_qc2 %>% filter(station == "EMM"))+
  geom_point(aes(datetime, wt, color = flagged), size = 1) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))
plotly::ggplotly(emm_wt)

sal_wt <- ggplot(wt_qc2 %>% filter(station == "SAL"))+
  geom_point(aes(datetime, wt, color = flagged), size = 1) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))
plotly::ggplotly(sal_wt)

(cll_wt <- ggplot(wt_qc2 %>% filter(station == "CLL"))+
    geom_point(aes(datetime, wt, color = flagged), size = 1) +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
    scale_color_manual(values = c("magenta", "black"))) 
plotly::ggplotly(cll_wt)

(dmc_wt <- ggplot(wt_qc2 %>% filter(station == "DMC"))+
    geom_point(aes(datetime, wt, color = flagged), size = 1) +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
    scale_color_manual(values = c("magenta", "black")))
plotly::ggplotly(dmc_wt)

(uni_wt <- ggplot(wt_qc2 %>% filter(station == "UNI"))+
    geom_point(aes(datetime, wt, color = flagged), size = 1) +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
    scale_color_manual(values = c("magenta", "black")))
plotly::ggplotly(uni_wt)


### Write --------------------------
saveRDS(wt_qc2, here("data/data_clean/wt_WY2025_flagged.rds"))
saveRDS(wt_clean, here("data/data_clean/wt_WY2025_clean.rds"))





#######################################################################

## EC DATA

########################################################################


## Electrical Conductivity -----------------------

# Load EC data and rearrange
# Data are in daily format
ec_data <- read_csv("data/data_clean/EC_USBR_2025_Daily_20260814.csv")

ec_data <- ec_data %>% 
  mutate(date = date(OBSERV_DATE),
         year = year(date))

ec_data_long <- pivot_longer(ec_data, cols = 2:8, names_to = "station_id")

saveRDS(ec_data_long, "data/data_raw/ec_raw_long_2025.rds")

# Clean data ---------------------------------

## Electrical Conductivity --------------------------
ec_raw0 <- readRDS(here("data/data_raw/ec_raw_long_2025.rds"))%>%
  mutate(date = ymd(OBSERV_DATE)) %>%
  #filter(date > ymd("2024-09-30")) #%>%
  rename(station = station_id)

# Complete data set with missing timestamps

# Build full 15-min sequence
full_time <- tibble(
  date = seq(from = start_date,
                 to   = end_date,
                 by   = "1 day")
)

# Combine station and time combos
sta_time <- crossing(stations, full_time) %>%
  rename(station = stations)

# Join with EC data
ec_raw <- sta_time %>%
  left_join(ec_raw0, by = c("station", "date")) %>%
  mutate(station_d1641 = case_when(station == "CLL" ~ "C2",
                                   station== "SAL" ~ "C4",
                                   station == "CNT" ~ "C5",
                                   station == "UNI" ~ "C8",
                                   station == "PCT" ~ "C14",
                                   station == "CCS" ~ "C19",
                                   station == "EMM" ~ "D22",
                                   station == "DMC" ~ "DMC1")) %>%
  mutate(date = date(date),
         month = factor(month(date))) %>%
  select(station,
         station_d1641,
         spc = value,
         date, month)

### Examine range of values
ggplot(ec_raw %>% filter(!is.na(month))) + 
  geom_boxplot(aes(month, spc, fill = month))+
  facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "EC (uS/cm)")+
  theme_bw() + 
  theme(legend.position = "bottom")

# see the outliers well here
ggplot(ec_raw %>% filter(!is.na(date))) + 
  geom_point(aes(date, spc, fill = station))+
  facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "EC (uS/cm)")+
  theme_bw() + 
  theme(legend.position = "bottom")

# Look at variability of data 
ec_sd <- ec_raw %>%
  filter(!is.na(ec)) %>%
  group_by(station, month, date) %>%
  summarize(sd = sd(ec, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(station, month) %>%
  summarize(mean_sd = mean(sd, na.rm = TRUE))

ggplot(ec_sd) + 
  geom_point(aes(month, mean_sd, color = station))+
  geom_line(aes(month, mean_sd, color = station, group = station)) +
  facet_wrap(~station, scales= "free_y")+
  theme_bw()






### QA/QC --------------
#### Attempt 1 ---------------
# tried a few rules here but too inconsistent in either flagging too many or too few points. 

# all in one dataframe: 
spc_qc <- ec_raw %>%
  # belowzero 
  mutate(range = if_else(spc <=5 | spc>25000, 1L, 0L))%>%
  # duplicated
  group_by(station) %>%
  mutate(
    spc_lag = lag(spc, 1),
    same = case_when(
      is.na(spc) & is.na(spc_lag) ~ 0L,   # both NA: treat as a "change" (breaks the run)
      is.na(spc) | is.na(spc_lag) ~ 0L,   # one NA, one real value: also breaks the run
      spc == spc_lag ~ 1L,
      TRUE ~ 0L
    ),
    issame = cumsum(same == 0L)
  ) %>%
  ungroup() %>%
  group_by(station, issame) %>%
  mutate(flag = sum(same) + 1) %>%
  ungroup() %>%
  mutate(repeated = ifelse(flag > 36, 1L, 0L)) %>% 
  group_by(station) %>%
  mutate(rmean = zoo::rollapply(spc, 36, mean, partial =TRUE),
         rsd = zoo::rollapply(spc, 36, sd, partial = TRUE)) %>%
  ungroup() %>% 
  mutate(outlier = ifelse(spc < rmean-3*rsd | spc > rmean+3*rsd, 1L, 0L)) %>%
  ungroup() %>%
  # spike test
  group_by(station) %>%
  mutate(spike = if_else(spc > 1.5 * lag(spc, 1, default = first(spc)) & 
                           spc > 1.5 * lead(spc, 1, default = first(spc)), 1L, 0L)) %>%
  ungroup() %>%
  # summarize
  mutate(flagged = if_else(range == 1L |  repeated == 1L | outlier == 1L | spike == 1L, "flagged", "not flagged")) 

# check what this looks like
ggplot(spc_qc)+ #%>% filter(!is.na(datetime))) + 
  geom_point(aes(date, spc, color = flagged, shape = flagged))+
  facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "spc (F)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom")


#### Attempt 2 (this is what I used) -----------
# manual removal of visual outliers
spc_qc2 <- spc_qc %>%
  mutate(vis_flag = case_when(station == "CNT" & spc < 50 ~ 1L,
                              station == "EMM" & date== ymd("2024-12-14") & spc > 1000 ~ 1L,
                              station == "EMM" & spc < 50 ~ 1L,
                              station == "CLL" & spc < 50 ~ 1L,
                              station == "CLL" & date== ymd("2024-12-14") & spc > 1000 ~ 1L,
                              station == "DMC" & date== ymd("2024-12-26") ~ 1L,
                              station == "DMC" & spc < 50 ~ 1L,
                              station == "PCT" & spc < 50 ~ 1L,
                              station == "SAL" & spc < 50 ~ 1L,
                              station == "UNI" & spc < 50 ~ 1L,
                              station == "UNI" & date >= ymd("2025-06-17") & date <= ymd("2025-06-19") & spc > 400 ~ 1L,
                              station == "UNI" & date >= ymd("2025-06-21") & date <= ymd("2025-06-24") & spc < 200 ~ 1L,
                              TRUE~0L)) %>%
  mutate(flagged = if_else(range == 1L | vis_flag == 1L, "flagged", "not flagged")) 

# look at this version of flags
ggplot(spc_qc2 %>% filter(!is.na(date))) + 
  geom_point(aes(date, spc, color = flagged, shape = flagged), size = 0.75)+
  facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "SpC (uS/cm)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom")

#### Cleaned ----------------
spc_check = spc_qc2 %>%
  filter(flagged == "flagged")

spc_clean <- spc_qc2 %>%
  filter(flagged == "not flagged")

ggplot(spc_clean %>% filter(!is.na(date))) + 
  geom_point(aes(date, spc))+
  facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "EC (uS/cm)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  theme_bw() + 
  theme(legend.position = "bottom")

#### Individual stations-------------------
cnt_spc <- ggplot(spc_qc %>% filter(station == "CNT"))+
  geom_point(aes(date, spc, color = flagged), size = 1) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))
plotly::ggplotly(cnt_spc)

emm_spc <- ggplot(spc_qc %>% filter(station == "EMM"))+
  geom_point(aes(date, spc, color = flagged), size = 1) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))
plotly::ggplotly(emm_spc)

cll_spc <- ggplot(spc_qc %>% filter(station == "CLL"))+
    geom_point(aes(date, spc, color = flagged), size = 1) +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
    scale_color_manual(values = c("magenta", "black"))
plotly::ggplotly(cll_spc)

dmc_spc <- ggplot(spc_qc %>% filter(station == "DMC"))+
  geom_point(aes(date, spc, color = flagged), size = 1) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))
plotly::ggplotly(dmc_spc)

pct_spc <- ggplot(spc_qc %>% filter(station == "PCT"))+
  geom_point(aes(date, spc, color = flagged), size = 1) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))
plotly::ggplotly(pct_spc)

sal_spc <- ggplot(spc_qc %>% filter(station == "SAL"))+
  geom_point(aes(date, spc, color = flagged), size = 1) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))
plotly::ggplotly(sal_spc)

uni_spc <- ggplot(spc_qc %>% filter(station == "UNI"))+
  geom_point(aes(date, spc, color = flagged), size = 1) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))
plotly::ggplotly(uni_spc)


### Write --------------------------
saveRDS(spc_qc2, here("data/data_clean/spc_WY2025_flagged.rds"))
saveRDS(spc_clean, here("data/data_clean/spc_WY2025_clean.rds"))






## Combine WT and EC data ------------------------


### Change WT data to daily mean
dailydata_wt <- wt_clean %>%
  mutate(wt_c = (wt - 32)/1.8) %>% 
  group_by(Site = station, Station = station_d1641, Date =date, Month = month) %>%
  summarize(Count = sum(!is.na(wt_c)),
            Value = round(mean(wt_c, na.rm = TRUE),1)) %>%
  ungroup() %>% 
  mutate(Analyte = "WaterTemperature")


### reorg SPC data

dailydata_spc <- spc_clean %>% 
  mutate(Count = 24, Analyte = "SpC") %>% 
  select(Site = station, Station = station_d1641, Date =date, Month = month,
         Count, Value= spc, Analyte)
  
### Combine df

data_all <- rbind.data.frame(dailydata_spc, dailydata_wt) %>% 
  mutate(Survey= "USBR") %>% 
  select(Survey, Date, Count, Value, Site, Station, Month, Analyte) # change to correct order


## Save as csv
write.csv(data_all, file = "data/data_clean/USBR_data_d1641_2025.csv", row.names = FALSE)


# # function to convert EC to SPC
# ec_to_spc <- function(EC, temp, alpha = 0.019) {
#   SpC <- EC / (1 + alpha * (temp - 25))
#   return(SpC)
# }
