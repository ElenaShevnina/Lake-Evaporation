# The supplement to Shevnina, Potes, Vihma and Naakka 2025
# Authorship: eshevnina@gmail.com
# phone: +358449185446
# April, 2025

library(stringr)
library(dplyr)
library(tidyr) 
library(ggplot2)
library(lubridate)

############## IRGASON DATA ##############################################################################################
### the fluxes and evaporation after Potes et al., 2017  ###########################################################
# Date Time, is the time stamp (TS) with the format “yy/mm/dd HH:MM:SS” [1]
# npoints, is number of measurements [2]
# u_starr (m/s), taur Kg(m s^2) are u_star and momentum flux [3], [4]
# Hcr (W/m^2) is a sensible heat flux humidity [5]
# LE_wplr (W/m^2), is a latent heat flux //  J m-2 s-1 [6]
# Evap (L/m^2) is an evaporation by Potes M. ## equal to mm of evaporation [7]
# Fc_wplr (mg/(m^2 s)), Fc_wpl_umol (umol/(m^2 s)), is carbon flux with WPL correction, mg/(m^2 s) and umol/(m^2 s) [8], [9]
# CO2_conc (mg/m^3), CO2_PPM (PPM), CO2_sig_str (arb), are CO2 concentration (two unit) and the Irgason’s signal strength [10], [11]
# H2O_conc (g/m^3), H2O_sig_str (arb), are H2O concentration and the Irgason’s signal strength 
# Temp_amb (C), Amb_Press (kPa), are ambient air temperature and atmospheric pressure
# wind_speed (m/s), wind dir_sonic (degree), are wind speed and direction
# obukhov (m), zeta (arb), rou (m), sigmaw_ustar (arb), are Obuhkov length, stability parameter, roughness
# Xmax (m), XR90 (m), are 90% of footprint (m) and maximum of footprint (m)

# To estimate Evaporation:
# L = 2.5x10 6 J kg water–1 = latent heat of vaporization; //Stull, 2015 p. 108
# m2 = 100 00 cm2
# Kg = 1000 g /ρ =1000 cm3 (ρ, density of water = 1 g/cm3 ) 
#  ET is Latent heat energy divided latent heat of vaporization of water: ET = LE/L = J m-2 s-1 /2.5x10 6 J kg-1 = J m-2 s-1 / 2.5x10 6 J kg-1 = J kg /2.5x10 6 J m2 s 
# Substituting for kg of water and m2 to cm2:  = J 1000 cm3 /2.5 MJ 10000 cm2 s = J cm3 /2.5x10 6 J 10 cm2 s = cm /2.5x10 6 s = mm/2.5x10 6 s = 0.0000004 mm/s
# In 30 minutes, = 1800 * 0.0000004 mm/s =0.00072 mm / 30 min 

###########################################################################################################
####################      Functions #############

#################### to read the raw data
to_read_Irgason<- function(Irgason)  ## the argument is the path to file with the Irgason measurements
{
 # variable_unit <- scan(Irgason, nlines = 18, sep=",", what = character())
  dat1 <- read.csv(Irgason, skip = 18, header = FALSE, fill=TRUE)
  dat1 <- dat1[,1:24]
  colnames(dat1) <- c("TS","npoints","u_star","taur","Hcr","LE_wplr", "Evap","Fc_wplr","Fc_wpl_umol","CO2_conc","CO2_PPM","CO2_sig_str",
                      "H2O_conc","H2O_sig_str","Temp_amb","Amb_Press","wind_speed","wind_dir_sonic", "obukhov","zeta","rou","sigmaw_ustar","Xmax","XR90")
  
  dat1$Timestamp_UTC <- as.POSIXct(strptime(substr(dat1$TS,1,15),tz="UTC",format="%y/%m/%d %H:%M") )   # to create proper UTC Timestamp, values corresponded to the end of time interval
  return(dat1)
}

# to read the Hobo data
# in function! the argument is the pahh to file with the HOBO mearuremets
to_read_hobo<- function(Indat4)
{
  dat4 <- read.table(Indat4,sep=",",skip=2,header=F)                             #### HOBO RAW DATA
  dat4 <- dat4[,2:5]
  colnam <- c("Timestamp","Abs_Press","TW","WL")
  colnames(dat4) <- colnam
  
  # to correct the date which is given with am (ap.) and pm (ip.)
  timestamp1 <- substr(dat4$Timestamp,1,8)
  timestamp2 <- substr(dat4$Timestamp,14,24)
  timestamp3<-str_replace(timestamp2, "ip","PM")
  timestamp3<-str_replace(timestamp3, "ap","AM")
  dat4$Timestamp <- paste(timestamp1,timestamp3,sep=" ")
  rm(timestamp1,timestamp2,timestamp3)
  dat4$Timestamp_EET <- as.POSIXlt(dat4$Timestamp,tz="Europe/Helsinki",format="%m.%d.%y %I.%M.%S %p")
  dat4$Timestamp_UTC <- dat4$Timestamp_EET-2*60*60      #original data given in EET, needs to be converted into UTC
  #dat4$Timestamp_UTC <- dat4$Timestamp_EET
  
  # whole observation period
  dat4$one_min<-strptime(dat4$Timestamp_UTC, format="%Y-%m-%d %H:%M:%S")
  hobo24hour<- aggregate(list(TW_mean = dat4$TW),list(Timestamp= cut(dat4$one_min, "24 hour" )), mean)
  min<-aggregate(list(TW_min = dat4$TW),list(Timestamp= cut(dat4$one_min, "24 hour" )), min)
  max<-aggregate(list(TW_max = dat4$TW),list(Timestamp= cut(dat4$one_min, "24 hour" )), max)
  hobo24hour$newTS<-as.Date(hobo24hour$Timestamp, "%Y-%m-%d")
  hobo24hour$TW_min<-min$TW_min
  hobo24hour$TW_max<-max$TW_max
  
  return(hobo24hour)
}


#################### to filter the fluxes 
### Flag conc and fluxes: to flag signal strengths lower than 0.7 + wind direction from the land + individually outliers
# dat1 is the raw data from the EC, out of to_read_irgason() function
# wd1 and wd2 are the min/max wind directions covering the footprint
# res1 is the filtered data, and res2 is list of statistics st1,2,3,4
to_filter_footprint <- function(dat1, wd1, wd2) 
{
  all <- nrow(dat1[!is.na(dat1$H2O_conc),])   
  
  filt <- nrow(dat1[(!is.na(dat1$H2O_sig_str) & dat1$H2O_sig_str >= 0.7),])       # st1: to filter out by the signal strengths (lower than 0.7)
  st1 <- (all-filt)*100/all # 0.11 %
  
  filt <- nrow(dat1[(dat1$wind_dir_sonic >= wd1 & dat1$wind_dir_sonic <= wd2),])   # st2 to filter out the winds from the land 
  st2 <- (all-filt)*100/all  # 17.8 %
  
  # to filter out by all above, and to get the percentage of data for the  Evaporation (E) >[st3]
  filt <- nrow(dat1[!is.na(dat1$Evap) & dat1$H2O_sig_str >= 0.7 & (dat1$wind_dir_sonic >= wd1 & dat1$wind_dir_sonic <= wd2),])
  st3<- 100 -((all-filt)*100/all)
  
  
  # to filter all
  dat1$H2O_conc_filter <- ifelse(dat1$H2O_sig_str >= 0.7 &
                                   (dat1$wind_dir_sonic >= wd1 & dat1$wind_dir_sonic <= wd2),
                                 dat1$H2O_conc,
                                 NA)
  dat1$LE_filter <- ifelse(dat1$H2O_sig_str >= 0.7 &
                             (dat1$wind_dir_sonic >=wd1 & dat1$wind_dir_sonic <= wd2),
                           dat1$LE_wplr,
                           NA)
  
  dat1$Evap_filter <- ifelse(dat1$H2O_sig_str >= 0.7 &
                               (dat1$wind_dir_sonic >= wd1 & dat1$wind_dir_sonic <= wd2),
                             dat1$Evap,
                             NA)
  
  dat1$H_filter <- ifelse((dat1$wind_dir_sonic >= wd1 & dat1$wind_dir_sonic <= wd2),
                          dat1$Hcr,
                          NA)
  dat1$tau_filter <- ifelse((dat1$wind_dir_sonic >= wd1 & dat1$wind_dir_sonic <= wd2),
                            dat1$taur,
                            NA)
  
  dat1$u_star_filter <- ifelse((dat1$wind_dir_sonic >= wd1 & dat1$wind_dir_sonic <= wd2),
                               dat1$u_star,
                               NA)
  
  #Outputs
  res1<- dat1
  numbers<-c(st1,st2, st3)
  cat(numbers, sep = " ")
  return (res1)
  
} #

####################to calculate daily min/mean/max temperature of air and water
# ECH is dataframe merged EC and LSWT 
tmp_stat_AirWater <- function(ECH) 
{
  t24hour <- aggregate(list(TW_mean = ECH$TW),list(Timestamp= cut(ECH$Timestamp_UTC, "24 hour" )), mean, na.rm=TRUE)
  t24hour$newTS <- as.Date(t24hour$Timestamp, "%Y-%m-%d")
  
  LWSTmin <- aggregate(list(TW_min = ECH$TW),list(Timestamp= cut(ECH$Timestamp_UTC, "24 hour" )), min, na.rm=TRUE)
  LWSTmax <- aggregate(list(TW_max = ECH$TW),list(Timestamp= cut(ECH$Timestamp_UTC, "24 hour" )), max, na.rm=TRUE)
  
  ATmin <- aggregate(list(AT_min = ECH$Temp_amb),list(Timestamp= cut(ECH$Timestamp_UTC, "24 hour" )), min, na.rm=TRUE)
  ATmax <- aggregate(list(AT_max = ECH$Temp_amb),list(Timestamp= cut(ECH$Timestamp_UTC, "24 hour" )), max, na.rm=TRUE)
  ATmean <- aggregate(list(AT_mean = ECH$Temp_amb),list(Timestamp= cut(ECH$Timestamp_UTC, "24 hour" )), mean, na.rm=TRUE)
  
  t24hour$TW_min <- LWSTmin$TW_min
  t24hour$TW_max <- LWSTmax$TW_max
  
  t24hour$AT_min <- ATmin$AT_min
  t24hour$AT_max <- ATmax$AT_max
  t24hour$AT_mean <- ATmean$AT_mean
  
  return(t24hour)
}


tmp_stat_HcLe_flux <- function(ECH) 
{
  t24hour <- aggregate(list(Hcr_mean = ECH$Hcr),list(Timestamp= cut(ECH$Timestamp_UTC, "24 hour" )), mean, na.rm=TRUE)
  t24hour$newTS <- as.Date(t24hour$Timestamp, "%Y-%m-%d")
  
  Hcr_min <- aggregate(list(Hcr_min = ECH$Hcr),list(Timestamp= cut(ECH$Timestamp_UTC, "24 hour" )), min, na.rm=TRUE)
  Hcr_max <- aggregate(list(Hcr_max = ECH$Hcr),list(Timestamp= cut(ECH$Timestamp_UTC, "24 hour" )), max, na.rm=TRUE)
  
  Le_min <- aggregate(list(Le_min = ECH$LE_wplr),list(Timestamp= cut(ECH$Timestamp_UTC, "24 hour" )), min, na.rm=TRUE)
  Le_max <- aggregate(list(Le_max = ECH$LE_wplr),list(Timestamp= cut(ECH$Timestamp_UTC, "24 hour" )), max, na.rm=TRUE)
  Le_mean <- aggregate(list(Le_mean = ECH$LE_wplr),list(Timestamp= cut(ECH$Timestamp_UTC, "24 hour" )), mean, na.rm=TRUE)
  
  t24hour$Hc_min <- Hcr_min$Hcr_min
  t24hour$Hc_max <- Hcr_max$Hcr_max
  
  t24hour$Le_min <- Le_min$Le_min
  t24hour$Le_max <- Le_max$Le_max
  t24hour$Le_mean <- Le_mean$Le_mean
  
  return(t24hour)
}

############################################################################################################################
### to calculate the relative humidity from the H2O concentration
# Tair is air temperature, C
# P is atmospheric pressure, kPa
# wv is H2O concentration,  g/m3
to_calulate_RH <- function (Tair, wv) {
  R <- 287.05                                                                     ###the universal gas constant for a dry air, J kg-1 K-1
  SatPre <- 0.6113 * exp(17.27*Tair/(Tair + 237.3))                             # Tetens formula in Stull, 2007, p. 89
  SpecHumid <- R*wv*0.001*(Tair + 273.15) / 622
  RH <- 100*SpecHumid/SatPre
  return(RH)
}
################################################################################################################################
######################### EOF: functions #######################################################################################



######################## DATA 30 min: EC (Ts, WS, RH) + HOBO / iButton (LSWT) ##############################
#############################################################################################################
### EC system on shore LAKE ZUB
Irgason<-"/home/shevnina/Manuscripts/2022/Antarctica_evaporation/Gryosphere2021/new_supplements/20180101_20180207_EC_FLUX.txt"  
Irg_ZB<-to_read_Irgason(Irgason) 
Irg_ZB$wind_dir_sonic <- (Irg_ZB$wind_dir_sonic+43)%%360                                     #correct for the  wind direction: added by Daniela Franz

######################### filtering the fluxes as suggested by Daniela Franz ###############################################
### Flag conc and fluxes: to flag signal strengths lower than 0.7 + wind direction from the land + individually outliers
### and to count the percentage of the filtered data
dat1 <- Irg_ZB
all <- nrow(dat1[!is.na(dat1$H2O_conc),])    # to account to the percentage of data filtered due to the signal strengths (lower than 0.7)
filt <- nrow(dat1[(!is.na(dat1$H2O_sig_str) & dat1$H2O_sig_str >= 0.7),])
(all-filt)*100/all # 0.11 %

filt <- nrow(dat1[(dat1$wind_dir_sonic >= 105 & dat1$wind_dir_sonic <= 240),])            # to wind direction from the land
(all-filt)*100/all  # 17.8 %

# to account the percentage of the data after the filtering
# total after the Latent sensible flux and Evaporation
filt <- nrow(dat1[!is.na(dat1$Evap) & dat1$H2O_sig_str >= 0.7 & (dat1$wind_dir_sonic >= 105 & dat1$wind_dir_sonic <= 240),])
100-((all-filt)*100/all)

dat1$Evap_filter <- ifelse(dat1$H2O_sig_str >= 0.7 &
                             (dat1$wind_dir_sonic >= 105 & dat1$wind_dir_sonic <= 240),
                           dat1$Evap,
                           NA)

####################################################################################################################################
####  daily EC evaporation calculated with NA in the values: used further for the errors estimations:
double1<-dat1$Evap_filter  # by mean
double1[is.na(double1)]<-mean(dat1$Evap_filter, na.rm=T)
dat1$Evap_filter<-double1
rm(double1)

evap_ECbyMean<-aggregate(list(Evap=dat1$Evap_filter),list(Timestamp= cut(dat1$Timestamp,"1440 min")), sum, na.rm=T) ### daily sum
evap_ECbyMean$newTS<-as.Date(evap_ECbyMean$Timestamp, "%Y-%m-%d")
# sum(evap_ECbyMean$Evap) [1] 114.133, the daily evaporation measured // 115.8282

 Irg_ZB <- dat1

# to aggeragate to 24 hours
evap_ECbyMean<-aggregate(list(Evap=dat1$Evap_filter),list(Timestamp= cut(dat1$Timestamp,"1440 min")), sum, na.rm=T) ### daily sum
evap_ECbyMean$newTS<-as.Date(evap_ECbyMean$Timestamp, "%Y-%m-%d")
evap_daily_ZB <- evap_ECbyMean

# write to the expernal file ! indirect methods!
write.csv(evap_daily_ZB, "/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/result_tab/evap_daily_ZB.csv")
rm(dat1, evap_ECbyMean)

########################################################################
# AWS MAITRI: [AT, WS, RH], Dhote et al., 2021
name_aws_data <- "/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/meteo_sites/MeteoMaitri.csv"
aws_data <- read.csv(name_aws_data, header = TRUE, sep = ",")
aws_data$date_new <- as.Date(aws_data$Date,format="%m/%d/%Y")         


# LWST, Lake ZUB
setwd("/home/shevnina/Manuscripts/2022/Antarctica_evaporation/Gryosphere2021/new_supplements/")
Indat4 <- "wl_parwati.txt"                                                         #### ZUB
dat4 <- read.table(Indat4, sep=",", skip=2, header=F)                              #### HOBO RAW DATA
dat4 <- dat4[,2:5]
colnam <- c("Timestamp","Abs_Press","TW","WL")
colnames(dat4) <- colnam
# to correct the date which is given with am (ap.) and pm (ip.)
timestamp1 <- substr(dat4$Timestamp, 1, 8)
timestamp2 <- substr(dat4$Timestamp, 14, 24)
timestamp3<-str_replace(timestamp2, "ip","PM")
timestamp3<-str_replace(timestamp3, "ap","AM")
dat4$Timestamp <- paste(timestamp1, timestamp3,sep=" ")
rm(timestamp1,timestamp2, timestamp3)
dat4$Timestamp_EET <- as.POSIXlt(dat4$Timestamp,tz="Europe/Helsinki",format="%m.%d.%y %I.%M.%S %p")
dat4$Timestamp_UTC <- dat4$Timestamp_EET-2*60*60      #original data given in EET, needs to be converted into UTC
Hobo_ZB <- dat4   
rm(dat4)

######## relative humidity + saturation pressure of air and water 
Irg_ZB$RH <- to_calulate_RH(Irg_ZB$Temp_amb, Irg_ZB$H2O_conc)                           # Hoeltgebaum et al. (2020)
####################################################### Irgason data + Hobo data, 30 min
ECH_ZB <- merge(Irg_ZB, Hobo_ZB, by = "Timestamp_UTC")
# the LSWT by HOBO temperature sensor<-ECH_GL$TW 


#### the Supplement
write.csv(ECH_ZB,"/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/result_tab/ECH_ZB.csv")

# to write the output for the indirect methods
BA_ZB <- ECH_ZB[,c("Timestamp_UTC","RH","Temp_amb","wind_speed","Amb_Press","Evap","Evap_filter","TW")]
write.csv(BA_ZB,"/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/result_tab/BA_ZB.csv")
####################################################################################################################

##############################################################
## Sub-daily cycle
#### hourly aggregations
ECH_ZB_60m <- ECH_ZB %>%
  group_by(time = floor_date(Timestamp_UTC, "60 minutes")) %>%
  summarize(
    AT=  mean(Temp_amb, na.rm = TRUE),
    WS = mean(wind_speed, na.rm = TRUE),
    WD = mean(wind_dir_sonic, na.rm = TRUE),
    E60 = sum(Evap, na.rm = TRUE),                                 ## non-filtered data EvF = sum(Evap_filter, na.rm = TRUE),                          ## filtered data
    LSWT = mean(TW, na.rm = TRUE),
    RH = mean(RH, na.rm = TRUE))

ECH_ZB_60m$dWT_AT <- ECH_ZB_60m$LSWT - ECH_ZB_60m$AT
ECH_ZB_60m$es <-  0.6113 * exp((17.27 * ECH_ZB_60m$LSWT) / (237.3 + ECH_ZB_60m$LSWT))   # saturation vapor pressure, hPa
ECH_ZB_60m$ea <- 0.01*ECH_ZB_60m$RH*0.6113* exp((17.27*ECH_ZB_60m$AT)/(ECH_ZB_60m$AT + 237.3))       # ea in kPa
ECH_ZB_60m$es_ea <- ECH_ZB_60m$es - ECH_ZB_60m$ea


# !non-filtered data !
df <-  ECH_ZB_60m
df <- df %>%
  mutate(TimeOfDay = format(time, "%H:%M:%S"))

# Convert the data to long format for easy plotting of multiple variables
df_long <- df %>%
  pivot_longer(cols = c(AT, LSWT, es_ea, WS, E60), names_to = "Variable", values_to = "Value")

##################################################### FIGURES, fig. 8 and 9 in Shevnina et al., 2025
###############################################
################### Fig. 8 a
# Boxplot of intradaily cycle for multiple variables by Copilot
p_AT_hour <- ggplot(df, aes(x = TimeOfDay, y = AT)) +
  geom_boxplot(fill = "coral", color = "black") +
  labs(x = "Time of Day", y = "AT", title = "") +
  scale_y_continuous(breaks = seq(-10, 5, by = 5), limits = c(-10,5)) +                  # Customize y-axis tick intervals
  theme_minimal() +
  theme(
  axis.text.x = element_text(angle = 90, hjust = 1, size = 16),        # Increase font size for x-axis text
  axis.text.y = element_text(size = 16),                               # Increase font size for y-axis text
  axis.title.x = element_text(size = 18),                              # Increase font size for x-axis label
  axis.title.y = element_text(size = 18),                              # Increase font size for y-axis label
  legend.position = "none"
)
#print(p_AT_hour)
ggsave("/home/shevnina/Manuscripts/2025/Ant/april/InDayAT_ZB.png", plot = p_AT_hour, width = 10, height = 8, dpi = 300)

p_LSWT_hour <- ggplot(df, aes(x = TimeOfDay, y = LSWT)) +
  geom_boxplot(fill = "cyan", color = "black") +
  labs(x = "Time of Day", y = "LWST", title = "") +
  scale_y_continuous(breaks = seq(0, 10, by = 5), limits = c(0,10)) +                  # Customize y-axis tick intervals
  theme_minimal() +
  theme(
  axis.text.x = element_text(angle = 90, hjust = 1, size = 16),        # Increase font size for x-axis text
  axis.text.y = element_text(size = 16),                               # Increase font size for y-axis text
  axis.title.x = element_text(size = 18),                              # Increase font size for x-axis label
  axis.title.y = element_text(size = 18),                              # Increase font size for y-axis label
  legend.position = "none"
)
#print(p_LSWT_hour)
ggsave("/home/shevnina/Manuscripts/2025/Ant/april/InDayLSWT_ZB.png", plot = p_LSWT_hour, width = 10, height = 8, dpi = 300)

############### Fig. 9 a
# Boxplot of intradaily cycle for multiple variables by Copilot
p_SPD_hour <- ggplot(df, aes(x = TimeOfDay, y = es_ea)) +
  geom_boxplot(fill = "lightsalmon", color = "black") +
  labs(x = "Time of Day", y = "Saturation deficit", title = "") +
  scale_y_continuous(breaks = seq(-0.25, 0.8, by = 0.25), limits = c(-0.25,0.8)) +                  # Customize y-axis tick intervals
  theme_minimal() +
  theme(
  axis.text.x = element_text(angle = 90, hjust = 1, size = 16),        # Increase font size for x-axis text
  axis.text.y = element_text(size = 16),                               # Increase font size for y-axis text
  axis.title.x = element_text(size = 18),                              # Increase font size for x-axis label
  axis.title.y = element_text(size = 18),                              # Increase font size for y-axis label
  legend.position = "none"
)
#print(p_SPD_hour)
ggsave("/home/shevnina/Manuscripts/2025/Ant/april/InDaySPD_ZB.png", plot = p_SPD_hour, width = 10, height = 8, dpi = 300)


p_WS_hour <- ggplot(df, aes(x = TimeOfDay, y = WS)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(x = "Time of Day", y = expression(WS~(m~s^-1)), title = "") +
  scale_y_continuous(breaks = seq(0, 15, by = 5), limits = c(0,15)) +
  theme_minimal() +
  theme(
  axis.text.x = element_text(angle = 90, hjust = 1, size = 16),  # Increase font size for x-axis text
  axis.text.y = element_text(size = 16),                        # Increase font size for y-axis text
  axis.title.x = element_text(size = 18),                       # Increase font size for x-axis label
  axis.title.y = element_text(size = 18),                       # Increase font size for y-axis label
  legend.position = "none"
)
#print(p_WS_hour)
ggsave("/home/shevnina/Manuscripts/2025/Ant/april/InDayWS_ZB.png", plot = p_WS_hour, width = 8, height = 6, dpi = 300)


# Boxplot of intradaily cycle for multiple variables by Copilot
p_eva_hour <- ggplot(df, aes(x = TimeOfDay, y = E60)) +
  geom_boxplot(fill = "grey", color = "black") +
  labs(x = "Time of Day", y = expression(Evaporation~(L~m^-2)), title = "") +
 scale_y_continuous(breaks = seq(-0.05, 0.15, by = 0.05), limits = c(-0.05,0.15)) +
  theme_minimal() +
  theme(
  axis.text.x = element_text(angle = 90, hjust = 1, size = 16),  # Increase font size for x-axis text
  axis.text.y = element_text(size = 16),                        # Increase font size for y-axis text
  axis.title.x = element_text(size = 18),                       # Increase font size for x-axis label
  axis.title.y = element_text(size = 18),                       # Increase font size for y-axis label
  legend.position = "none"
)
#print(p_eva_hour)
ggsave("/home/shevnina/Manuscripts/2025/Ant/april/InDayEVA_ZB.png", plot = p_eva_hour, width = 8, height = 6, dpi = 300)
############################################
### EOF of sub-daily cycle
#######################################################

########################################################################### FIGURES, Fig. 3
AWTemp_ZB <-tmp_stat_AirWater(ECH_ZB)
# min/mean/max of LWST, Lake ZUB + IRgason Air Temp [Dec, 2017 - Feb,2018]
pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/result_fig/LSWTbyHobo_TSplot_ZB.pdf")
par(family = "serif", cex = 1.0)
par(xpd = TRUE)
#par(mar = par()$mar + c(0, 0, 0, 7))
plot(AWTemp_ZB$newTS, AWTemp_ZB$TW_mean, type="l", xlim=c(as.Date("2017-12-01"),as.Date("2018-02-28")), ylab= "T, C", xlab="Date", ylim=c(-10,10),col="blue",lwd=2)
lines(AWTemp_ZB$newTS, AWTemp_ZB$TW_max, col="blue", lwd=0.5, lty=2)
lines(AWTemp_ZB$newTS, AWTemp_ZB$TW_min, col="blue",lwd=0.5, lty=2)

#par(new=TRUE)
lines(AWTemp_ZB$newTS, AWTemp_ZB$AT_mean,  col="red",lwd=2)
lines(AWTemp_ZB$newTS, AWTemp_ZB$AT_max, col="red", lwd=0.5, lty=2)
lines(AWTemp_ZB$newTS, AWTemp_ZB$AT_min, col="red",lwd=0.5, lty=2)
lines(aws_data$date_new, aws_data$AT, col = "red", lwd = 0.5)

abline(v=as.Date("2018-01-01"), col="black",lty=2, lwd=1)                # camera + IRGASON ON
abline(v=as.Date("2018-02-07"), col="black",lty=2, lwd=1)                # Irgason OFF

abline(h=0, col="black", lwd=0.5)
grid()
#legend("topright", inset = c(-0.4, 0), legend=c("HOBO", "iBUT", "PanoIMG", "IRG[ON/OFF]"), col=c("blue", "green" ,"red", "black"), lwd=1)
dev.off()


################################################### FIGURES, Fig. 4
HZ_daily_stat <-tmp_stat_HcLe_flux(ECH_ZB)
pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/result_fig/Hcr_Le_TSplot_ZB.pdf")
par(family = "serif", cex = 1.0)
par(xpd = TRUE)
par(mar = par()$mar + c(0, 0, 0, 7))

plot(HZ_daily_stat$newTS, HZ_daily_stat$Hcr_mean, type="l", xlim=c(as.Date("2017-12-01"),as.Date("2018-02-28")), ylab= "H, Wm-2", xlab="Date", ylim=c(-110,290),col="orange",lwd=1)
lines(HZ_daily_stat$newTS, HZ_daily_stat$Hc_max, col="orange", lwd=0.5, lty=2)
lines(HZ_daily_stat$newTS, HZ_daily_stat$Hc_min, col="orange",lwd=0.5, lty=2)

par(new=TRUE)
lines(HZ_daily_stat$newTS, HZ_daily_stat$Le_mean, type="l", xlim=c(as.Date("2017-12-01"),as.Date("2018-02-28")), col="green", lwd=1)
lines(HZ_daily_stat$newTS, HZ_daily_stat$Le_max, col="green", lwd=0.5, lty=2)
lines(HZ_daily_stat$newTS, HZ_daily_stat$Le_min, col="green",lwd=0.5, lty=2)

abline(v=as.Date("2018-01-01"), col="black",lty=2, lwd=1)                # camera + IRGASON ON
abline(v=as.Date("2018-02-07"), col="black",lty=2, lwd=1)                # Irgason OFF

abline(h=0, col="black", lwd=0.5)
grid()
dev.off()


