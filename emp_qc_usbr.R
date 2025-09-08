# WY 2024 Data QC - Reclamation EMP stations
# Catarina Pien (cpien@usbr.gov)
# This code cleans up CDEC data for EMP reporting 
# Lots of very specific QC rather than a general script

# Load packages------------------
library(here)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(readr)
library(sharpshootR)
library(scales)
library(dbscan)
library(zoo)

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

# Download data --------------------------
## Electrical Conductivity -----------------------
stations <- c("CLL", "SAL", "CNT", "UNI", "PCT", "EMM", "DMC")
start_date <- as.POSIXct("2023-10-01 00:00:00", tz = "UTC")
end_date   <- as.POSIXct("2024-09-30 23:45:00", tz = "UTC")

ec_data <- lapply(stations,
                   function(x) {
                     CDECquery(id = x,sensor = 100, interval = "E", start = start_date,end = end_date)
                     })

ec_raw <- bind_rows(ec_data) %>%
  mutate(datetime = format(as.POSIXct(datetime), "%Y-%m-%d %H:%M:%S", tz = "US/Pacific"),
         date = date(datetime))

saveRDS(ec_raw, "data/data_raw/ec_raw_2024.rds")

## Water Temperature --------------------------

wt_data <- lapply(stations,
                   function(x) {
                     CDECquery(id = x,sensor = 25, interval = "H", start = start_date,end = end_date)
                     })
wt_raw <- bind_rows(wt_data)%>%
  mutate(datetime = format(as.POSIXct(datetime), "%Y-%m-%d %H:%M:%S"),
         date = date(datetime),
         year = year(date))

saveRDS(wt_raw, "data/data_raw/wt_raw_2024.rds")

# Clean data ---------------------------------
## Electrical Conductivity --------------------------
ec_raw0 <- readRDS(here("data/data_raw/ec_raw_2024.rds"))%>%
  mutate(datetime = ymd_hms(datetime)) %>%
  filter(date > ymd("2023-09-30")) %>%
  rename(station = station_id)

# Complete data set 

# Build full 15-min sequence
full_time <- tibble(
  datetime = seq(from = start_date,
                 to   = end_date,
                 by   = "15 min")
)

sta_time <- crossing(stations, full_time) %>%
  rename(station = stations)

# Join with your data
ec_raw <- sta_time %>%
  left_join(ec_raw0, by = c("station", "datetime")) %>%
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
         ec = value,
         datetime, date, month)

### Examine range of values
ggplot(ec_raw %>% filter(!is.na(month))) + 
  geom_boxplot(aes(month,ec, fill = month))+
  facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "EC (uS/cm)")+
  theme_bw() + 
  theme(legend.position = "bottom")

# see the outliers well here
ggplot(ec_raw %>% filter(!is.na(datetime))) + 
  geom_point(aes(datetime,ec, fill = station))+
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

library(dbscan)
# all in one dataframe: 
ec_qc <- ec_raw %>%
  # belowzero 
  mutate(belowzero = if_else(station != "CLL" & ec <=0, 1L, if_else(ec<0, 1L, 0L)))%>%
  # missing data 
  group_by(station, date) %>%
  mutate(missingvals = sum(is.na(ec)),
         missing_flag = if_else(missingvals/96>=0.7, 1L, 0L)) %>%
  ungroup() %>%
  # duplicated
  group_by(station) %>%
  mutate(same = ifelse(ec == lag(ec, 1, default = 0), 1L, 0L),
         issame = cumsum(same == 0L)) %>%
  ungroup() %>%
  group_by(station, issame) %>%
  mutate(flag = sum(same) + 1) %>%
  ungroup() %>%
  mutate(repeated = ifelse(flag>36, 1L, 0L)) %>%
  # # MAD 
  # group_by(station) %>%
  # mutate(
  #   rmed = zoo::rollapply(ec, 336, median, partial = TRUE, na.rm = TRUE),
  #   mad_val = mad(ec, constant = 1, na.rm = TRUE),
  #   mad_outlier = if_else(abs(ec - rmed) > 3 * mad_val, 1L, 0L) # 5×MAD rule
  # ) %>%
  # ungroup()%>%
  # local outlier factor
  # group_by(station) %>%
  # mutate(
  #   lof = lof(as.matrix(ec), minPts = 1440),
  #   lof_outlier = if_else(lof > 6, 1L, 0L)
  # ) %>%
  # ungroup()%>%
  # outlier 
  group_by(station) %>%
  mutate(rmean = zoo::rollapply(ec, 1440, mean, partial =TRUE),
         rsd = zoo::rollapply(ec, 1440, sd, partial = TRUE)) %>%
  ungroup() %>% 
  mutate(outlier = ifelse(ec < rmean-4*rsd | ec > rmean+4*rsd, 1L, 0L)) %>%
  ungroup() %>%
  
  # filter(complete.cases(ec)) %>%
  # mutate(repeated = if_else((lag(ec,1) == ec & lag(ec,2) == ec & lag(ec,3) == ec)| 
           # (lead(ec,1) == ec & lead(ec,2) == ec & lead(ec,3) == ec), 1L, 0L)) %>%
  # ungroup()%>%
  # spike test
  group_by(station) %>%
  mutate(spike = if_else(ec > 1.5 * lag(ec, 1, default = first(ec)) & 
                          ec > 1.5 * lead(ec, 1, default = first(ec)), 1L, 0L)) %>%
  ungroup() %>%
  # summarize
  mutate(flagged = if_else(belowzero == 1L | missing_flag == 1L |repeated == 1L | outlier == 1L | spike == 1L, "flagged", "not flagged")) 

# look at flagged outliers
ggplot(ec_qc %>% filter(!is.na(datetime))) + 
  geom_point(aes(datetime,ec, color = flagged, shape = flagged))+
  facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "EC (uS/cm)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom")

#### Attempt 2 (this is what I used) -----------
# manual removal of visual outliers
ec_qc2 <- ec_qc %>%
  mutate(vis_flag = case_when(station == "CLL" & ec <= 150 ~ 1L,
                          station == "CLL" & month %in% c(11,12) & ec < 120 ~1L,
                          station == "CLL" & date== ymd("2023-11-29") & ec == 1074 ~ 1L,
                          station == "CNT" & ec>1000 ~1L,
                          station == "CNT" & date >= ymd("2023-10-17") & date<=ymd("2023-10-20") & ec <270 ~ 1L,
                          station == "CNT" & date >=ymd("2023-11-20") & date <= ymd("2023-12-19") ~ 1L,
                          station == "DMC" & ec > 2000~1L,
                          station == "CNT" & ec <160 ~ 1L,
                          station == "SAL" & ec > 950 ~1L,
                          station == "PCT" & ec < 10 ~ 1L,
                          station == "DMC" & ec < 100 ~ 1L,
                          station == "DMC" & date == ymd("2024-04-18") & ec < 200 ~ 1L,
                          TRUE~0L)) %>%
  mutate(flagged = if_else(belowzero == 1L | missing_flag == 1L | vis_flag == 1L, "flagged", "not flagged")) 

# look at this version of flags
ggplot(ec_qc2 %>% filter(!is.na(datetime))) + 
  geom_point(aes(datetime,ec, color = flagged, shape = flagged), size = 0.75)+
  facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "EC (uS/cm)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom")

#### Cleaned ----------------
ec_check = ec_qc2 %>%
  filter(flagged == "flagged")

ec_clean <- ec_qc2 %>%
  filter(flagged == "not flagged")

ggplot(ec_clean %>% filter(!is.na(datetime))) + 
  geom_point(aes(datetime,ec))+
  facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "EC (uS/cm)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  theme_bw() + 
  theme(legend.position = "bottom")

#### Individual stations-------------------
cnt_ec <- ggplot(ec_qc2 %>% filter(station == "CNT"))+
  geom_point(aes(datetime, ec, color = flagged), size = 1) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))
plotly::ggplotly(cnt_ec)

emm_ec <- ggplot(ec_qc2 %>% filter(station == "EMM"))+
  geom_point(aes(datetime, ec, color = flagged), size = 1) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))
plotly::ggplotly(emm_ec)

(cll_ec <- ggplot(ec_qc2 %>% filter(station == "CLL"))+
  geom_point(aes(datetime, ec, color = flagged), size = 1) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black")))

plotly::ggplotly(cll_ec)

### Summary of data removed --------------

ec_qc2 %>%
  group_by(station) %>%
  summarize(removed = sum(flagged == "flagged"),
            values = n(),
            pct_removed = removed/values*100)


### Write --------------------------
saveRDS(ec_qc2, here("data/data_clean/ec_WY2024_flagged.rds"))
saveRDS(ec_clean, here("data/data_clean/ec_WY2024_clean.rds"))


## Water Temperature -------------------
wt_raw0 <- readRDS(here("data/data_raw/wt_raw_2024.rds"))%>%
  mutate(datetime = ymd_hms(datetime)) %>%
  filter(date > ymd("2023-09-30")) %>%
  rename(station = station_id)

hour_time <- tibble(
  datetime = seq(from = start_date,
                 to   = end_date,
                 by   = "1 hour")
)

sta_hour <- crossing(stations, hour_time) %>%
  rename(station = stations)

# Join with your data
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
  # missing data 
  group_by(station, date) %>%
  mutate(missingvals = sum(is.na(wt)),
         missing_flag = if_else(missingvals/24>=0.7, 1L, 0L)) %>%
  ungroup() %>%
  # duplicated
  group_by(station) %>%
  mutate(same = ifelse(wt == lag(wt, 1, default = 0), 1L, 0L),
         issame = cumsum(same == 0L)) %>%
  ungroup() %>%
  group_by(station, issame) %>%
  mutate(flag = sum(same) + 1) %>%
  ungroup() %>%
  mutate(repeated = ifelse(flag>36, 1L, 0L)) %>%
  # # MAD 
  # group_by(station) %>%
  # mutate(
  #   rmed = zoo::rollapply(wt, 336, median, partial = TRUE, na.rm = TRUE),
  #   mad_val = mad(wt, constant = 1, na.rm = TRUE),
  #   mad_outlier = if_else(abs(wt - rmed) > 3 * mad_val, 1L, 0L) # 5×MAD rule
  # ) %>%
  # ungroup()%>%
  # local outlier factor
  # group_by(station) %>%
  # mutate(
  #   lof = lof(as.matrix(wt), minPts = 1440),
  #   lof_outlier = if_else(lof > 6, 1L, 0L)
  # ) %>%
  # ungroup()%>%
  # outlier 
  group_by(station) %>%
  mutate(rmean = zoo::rollapply(wt, 360, mean, partial =TRUE),
         rsd = zoo::rollapply(wt, 360, sd, partial = TRUE)) %>%
  ungroup() %>% 
  mutate(outlier = ifelse(wt < rmean-3*rsd | wt > rmean+3*rsd, 1L, 0L)) %>%
  ungroup() %>%
  # filter(complete.cases(wt)) %>%
  # mutate(repeated = if_else((lag(wt,1) == wt & lag(wt,2) == wt & lag(wt,3) == wt)| 
  # (lead(wt,1) == wt & lead(wt,2) == wt & lead(wt,3) == wt), 1L, 0L)) %>%
  # ungroup()%>%
  # spike test
  group_by(station) %>%
  mutate(spike = if_else(wt > 1.5 * lag(wt, 1, default = first(wt)) & 
                           wt > 1.5 * lead(wt, 1, default = first(wt)), 1L, 0L)) %>%
  ungroup() %>%
  # summarize
  mutate(flagged = if_else(range == 1L | missing_flag == 1L| repeated == 1L | outlier == 1L | spike == 1L, "flagged", "not flagged")) 

# check what this looks like
ggplot(wt_qc %>% filter(!is.na(datetime))) + 
  geom_point(aes(datetime,wt, color = flagged, shape = flagged))+
  facet_wrap(~station, scales = "free_y", ncol = 2) +
  labs(y = "wt (uS/cm)")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("magenta", "black"))+
  theme_bw() + 
  theme(legend.position = "bottom")

#### Attempt 2 (this is what I used) -----------
# manual removal of visual outliers
wt_qc2 <- wt_qc %>%
  mutate(vis_flag = case_when(station == "SAL" & date== ymd("2023-12-20")& wt == 58.4 ~1L,
                              station == "EMM" & wt < 48~1L,
                              station == "EMM" & date >= ymd("2024-08-14") & date <=ymd("2024-08-15") & wt > 64~1L,
                              station == "EMM" & date >= ymd("2024-09-16") & date <=ymd("2024-09-17") & wt > 64~1L,
                              station == "PCT" & wt > 74~1L,
                              station == "DMC" & date >=ymd("2024-01-26") & date <= ymd("2024-03-03") ~1L,
                              station == "DMC" & date ==ymd("2024-04-17") & wt==76.9 ~1L,
                              station == "DMC" & date ==ymd("2024-04-17") & wt==76.9 ~1L,
                              station == "DMC" & date >=ymd("2024-05-14") & date <= ymd("2024-05-20") & wt >70~1L,
                              station == "DMC" & date >=ymd("2024-09-27") & date <= ymd("2024-09-29") & ((wt>74)|(wt<70)) ~1L,
                              station == "CLL" & month == 3 & wt>60~1L,
                              station == "CLL"  & wt<49.5~1L,
                              station == "CLL"  &date >=ymd("2024-03-22") & date <= ymd("2024-03-27") &wt<55~1L,
                              station == "CNT" &date >=ymd("2023-11-20") & date <= ymd("2023-12-19") ~ 1L,
                              TRUE~0L)) %>%
  mutate(flagged = if_else(range == 1L | missing_flag == 1L | vis_flag == 1L, "flagged", "not flagged")) 

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
(cnt_wt <- ggplot(wt_qc2 %>% filter(station == "CNT", range == 0L))+
  geom_point(aes(datetime, wt, color = "flagged"), size = 1) +
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

(cll_wt <- ggplot(wt_qc2 %>% filter(station == "CLL", range == 0L))+
    geom_point(aes(datetime, wt, color = flagged), size = 1) +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
    scale_color_manual(values = c("magenta", "black"))) 
plotly::ggplotly(cll_wt)

(dmc_wt <- ggplot(wt_qc2 %>% filter(station == "DMC"))+
    geom_point(aes(datetime, wt, color = flagged), size = 1) +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
    scale_color_manual(values = c("magenta", "black")))
plotly::ggplotly(dmc_wt)

plotly::ggplotly(cll_wt)

### Summary of data removed --------------
wt_qc2 %>%
  group_by(station) %>%
  summarize(removed = sum(flagged == "flagged"),
            values = n(),
            pct_removed = removed/values*100)


### Write --------------------------
saveRDS(wt_qc2, here("data/data_clean/wt_WY2024_flagged.rds"))
saveRDS(wt_clean, here("data/data_clean/wt_WY2024_clean.rds"))

## Combine data ------------------------

# function to convert EC to SPC
ec_to_spc <- function(EC, temp, alpha = 0.019) {
  SpC <- EC / (1 + alpha * (temp - 25))
  return(SpC)
}

alldata <- left_join(ec_clean %>% select(station, station_d1641, month, datetime, date, ec),
          wt_clean%>% select(station, station_d1641, datetime, date, wt)) %>%
  mutate(spc = ec_to_spc(ec, wt)) 

# Note: once we convert ec to spc we now 
# just have hourly values for spc since we only have hourly water temp
# data. Could fill downward to calculate more values, but probably not necessary 
# since this all gets averaged to daily anyways. 

dailydata <- alldata %>%
  group_by(station, station_d1641, date, month) %>%
  summarize(wt_n = sum(!is.na(wt)),
            wt_mean = round(mean(wt, na.rm = TRUE),1),
         spc_mean = round(mean(spc, na.rm = TRUE),1),
         spc_n = sum(!is.na(spc))) %>%
  ungroup()

dailydata_formatted <- dailydata %>%
  pivot_longer(cols = c(wt_mean, wt_n, spc_mean, spc_n),
               names_to = c("Analyte", "ValueType"),
               names_sep = "_") %>%
  pivot_wider(names_from = "ValueType",
              values_from = "value") %>%
  select(Date = date, Count = n, Value=mean, Site = station, Station = station_d1641, Month=month,
         Analyte) %>%
  mutate(Analyte = case_when(Analyte == "wt" ~ "WaterTemperature",
                             Analyte == "spc" ~ "SpC")) %>%
  arrange(Analyte, Site, Date)

# Look at each station
plot_vals <- function(CDEC_sta) {
  ggplot(dailydata_formatted %>% filter(Site == CDEC_sta))+
    geom_point(aes(Date, Value, color = Analyte)) +
    facet_wrap(~Analyte, scales = "free_y", nrow = 2) + 
    theme_bw()
}

plot_vals("CLL")


# Write dataset 

saveRDS(dailydata_formatted, here("data/data_clean/usbr_contdata_wy2024.rds"))
