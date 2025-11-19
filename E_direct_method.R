# The supplement to Shevnina, Potes, Vihma and Naakka 2025
# Authorship: eshevnina@gmail.com
# phone: +358449185446
# October, 2025

library(stringr)
library(dplyr)
library(tidyr) 
library(ggplot2)
library(lubridate)
library(patchwork)

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

#################### by Copilot
tmp_stat_MultiVars <- function(ECH) {
  t24hour <- aggregate(list(
    Amb_Press_mean = ECH$Amb_Press,
    wind_speed_mean = ECH$wind_speed,
    Temp_amb_mean = ECH$Temp_amb,
    RH_mean = ECH$RH
  ), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), mean, na.rm = TRUE)
  t24hour$newTS <- as.Date(t24hour$Timestamp, "%Y-%m-%d")

  Amb_Press_min <- aggregate(list(Amb_Press_min = ECH$Amb_Press), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), min, na.rm = TRUE)
  Amb_Press_max <- aggregate(list(Amb_Press_max = ECH$Amb_Press), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), max, na.rm = TRUE)

  wind_speed_min <- aggregate(list(wind_speed_min = ECH$wind_speed), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), min, na.rm = TRUE)
  wind_speed_max <- aggregate(list(wind_speed_max = ECH$wind_speed), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), max, na.rm = TRUE)

  Temp_amb_min <- aggregate(list(Temp_amb_min = ECH$Temp_amb), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), min, na.rm = TRUE)
  Temp_amb_max <- aggregate(list(Temp_amb_max = ECH$Temp_amb), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), max, na.rm = TRUE)

  RH_min <- aggregate(list(RH_min = ECH$RH), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), min, na.rm = TRUE)
  RH_max <- aggregate(list(RH_max = ECH$RH), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), max, na.rm = TRUE)

  t24hour$Amb_Press_min <- Amb_Press_min$Amb_Press_min
  t24hour$Amb_Press_max <- Amb_Press_max$Amb_Press_max

  t24hour$wind_speed_min <- wind_speed_min$wind_speed_min
  t24hour$wind_speed_max <- wind_speed_max$wind_speed_max

  t24hour$Temp_amb_min <- Temp_amb_min$Temp_amb_min
  t24hour$Temp_amb_max <- Temp_amb_max$Temp_amb_max

  t24hour$RH_min <- RH_min$RH_min
  t24hour$RH_max <- RH_max$RH_max

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

# to calculaye the saturation vapor pressure deficit SVD
# T is air temperature, C
# LWST is H2O concentration  C
to_calulate_SVD <- function (T, LWST) {
  SVD <- 0.6113 * exp(17.27*LWST/(LWST + 237.3)) - 0.6113 * exp(17.27*T/(T + 237.3))   # saturation vapor pressure deficit, hPa
  return(SVD)
}


# by Copilot
#' Aggregate ECH data to hourly means/sums
hourly_aggregate_ECH <- function(ECH) {
  ECH60m <- ECH %>%
    group_by(time = floor_date(Timestamp_UTC, "60 minutes")) %>%
    summarize(
      AT = mean(Temp_amb, na.rm = TRUE),
      WS = mean(wind_speed, na.rm = TRUE),
      WD = mean(wind_dir_sonic, na.rm = TRUE),
    #  E60 = sum(Evap, na.rm = TRUE),                   # mm per hour is 1 L m-2 per hour
      E60 = sum(Evap[Evap > 0], na.rm = TRUE),      # sum only positive values
      LSWT = mean(TW, na.rm = TRUE),
      RH = mean(RH, na.rm = TRUE),
      SVD = mean(SVD, na.rm = TRUE)
    )
  return(ECH60m)
}

tmp_stat_MultiVars <- function(ECH) {
  t24hour <- aggregate(list(
    Amb_Press_mean = ECH$Amb_Press,
    wind_speed_mean = ECH$wind_speed,
    Temp_amb_mean = ECH$Temp_amb,
    RH_mean = ECH$RH
  ), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), mean, na.rm = TRUE)
  t24hour$newTS <- as.Date(t24hour$Timestamp, "%Y-%m-%d")

  Amb_Press_min <- aggregate(list(Amb_Press_min = ECH$Amb_Press), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), min, na.rm = TRUE)
  Amb_Press_max <- aggregate(list(Amb_Press_max = ECH$Amb_Press), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), max, na.rm = TRUE)

  wind_speed_min <- aggregate(list(wind_speed_min = ECH$wind_speed), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), min, na.rm = TRUE)
  wind_speed_max <- aggregate(list(wind_speed_max = ECH$wind_speed), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), max, na.rm = TRUE)

  Temp_amb_min <- aggregate(list(Temp_amb_min = ECH$Temp_amb), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), min, na.rm = TRUE)
  Temp_amb_max <- aggregate(list(Temp_amb_max = ECH$Temp_amb), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), max, na.rm = TRUE)

  RH_min <- aggregate(list(RH_min = ECH$RH), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), min, na.rm = TRUE)
  RH_max <- aggregate(list(RH_max = ECH$RH), list(Timestamp = cut(ECH$Timestamp_UTC, "24 hour")), max, na.rm = TRUE)

  t24hour$Amb_Press_min <- Amb_Press_min$Amb_Press_min
  t24hour$Amb_Press_max <- Amb_Press_max$Amb_Press_max

  t24hour$wind_speed_min <- wind_speed_min$wind_speed_min
  t24hour$wind_speed_max <- wind_speed_max$wind_speed_max

  t24hour$Temp_amb_min <- Temp_amb_min$Temp_amb_min
  t24hour$Temp_amb_max <- Temp_amb_max$Temp_amb_max

  t24hour$RH_min <- RH_min$RH_min
  t24hour$RH_max <- RH_max$RH_max

  return(t24hour)
}
################################################################################################################################
######################### EOF: functions #######################################################################################



######################## DATA 30 min: EC (Ts, WS, RH) + HOBO / iButton (LSWT) ##############################
#############################################################################################################
### EC system on shore LAKE ZUB
#Irgason<-"/home/shevnina/Manuscripts/2022/Antarctica_evaporation/Gryosphere2021/new_supplements/20180101_20180207_EC_FLUX.txt"
Irgason<-"/home/shevnina/Manuscripts/2025/Ant/datasets/20171230_20180208_EC_FLUX.txt"
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

#write.csv(evap_daily_ZB, "/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/result_tab/evap_daily_ZB.csv")
write.csv(evap_daily_ZB, "/home/shevnina/Manuscripts/2025/Ant/september/result_tab/evap_daily_ZB.csv")

rm(dat1, evap_ECbyMean, all, filt)

########################################################################
# AWS MAITRI: [AT, WS, RH], Dhote et al., 2021
name_aws_data <- "/home/shevnina/Manuscripts/2025/Ant/datasets/MeteoMaitri.csv"
aws_data <- read.csv(name_aws_data, header = TRUE, sep = ",")
aws_data$date_new <- as.Date(aws_data$Date,format="%m/%d/%Y")         


# LWST, Lake ZUB
setwd("/home/shevnina/Manuscripts/2025/Ant/datasets/")
#Indat4 <- "wl_parwati.txt"
Indat4 <- "HOBO20172018.csv"  #### ZUB
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
rm(dat4, Indat4)

######## relative humidity + saturation pressure of air and water 
Irg_ZB$RH <- to_calulate_RH(Irg_ZB$Temp_amb, Irg_ZB$H2O_conc)
####################################################### Irgason data + Hobo data, 30 min
ECH_ZB <- merge(Irg_ZB, Hobo_ZB, by = "Timestamp_UTC")
ECH_ZB$SVD <- to_calulate_SVD(ECH_ZB$Temp_amb, ECH_ZB$TW)

#### the Supplement
write.csv(ECH_ZB,"/home/shevnina/Manuscripts/2025/Ant/september/result_tab/ECH_ZB.csv")
# to write the output for the indirect methods
BA_ZB <- ECH_ZB[,c("Timestamp_UTC","RH","Temp_amb","wind_speed","Amb_Press","Evap","Evap_filter","TW")]
write.csv(BA_ZB,"/home/shevnina/Manuscripts/2025/Ant/september/result_tab/BA_ZB.csv")
####################################################################################################################
# EOF of the observations on Lake ZUB
####################################################################################################################

##############################################################
# Sub-daily cycle
#### hourly aggregations
ECH_ZB_60m <- hourly_aggregate_ECH(ECH_ZB)
ECH_GL_60m <- hourly_aggregate_ECH(ECH_LG)

# Figure 9
# Boxplot of intradaily cycle for multiple variables by Copilot
# Diurnal cycle of E60, WS, SVD
# Combine and prepare data
df_GL <- ECH_GL_60m %>% mutate(Hour = format(time, "%H"), Site = "LG")
df_ZB <- ECH_ZB_60m %>% mutate(Hour = format(time, "%H"), Site = "ZB")
df_both <- bind_rows(df_GL, df_ZB)

# E60 boxplot
p_E60 <- ggplot(df_both, aes(x = Hour, y = E60, fill = Site)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  labs(x = NULL, y = expression(E~(mm~h^-1)), title = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "top")

# WS boxplot
p_WS <- ggplot(df_both, aes(x = Hour, y = WS, fill = Site)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  labs(x = NULL, y = expression(WS~(m~s^-1)), title = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "none")

# SVD boxplot
p_SVD <- ggplot(df_both, aes(x = Hour, y = SVD, fill = Site)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  labs(x = "Hour", y = "SVD, kPa", title = "") +
  theme_minimal() +
  theme(legend.position = "none")

# Combine vertically
combined_plot <- p_E60 / p_WS / p_SVD + plot_layout(heights = c(1, 1, 1))
ggsave('/home/shevnina/Manuscripts/2025/Ant/september/result_fig/InDay_E60_WS_SVD_LGZB.png', plot = combined_plot, width = 8, height = 16, dpi = 300)
### EOF of sub-daily cycle
#######################################################

########################################################################### FIGURES, Fig. 3
AWTemp_ZB <-tmp_stat_AirWater(ECH_ZB)
# min/mean/max of LWST, Lake ZUB + IRgason Air Temp [Dec, 2017 - Feb,2018]
pdf("/home/shevnina/Manuscripts/2025/Ant/september/result_fig/LSWTbyHobo_TSplot_ZB.pdf")
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

#abline(v=as.Date("2018-01-01"), col="black",lty=2, lwd=1)                # camera + IRGASON ON
#abline(v=as.Date("2018-02-07"), col="black",lty=2, lwd=1)                # Irgason OFF

abline(h=0, col="black", lwd=0.5)
grid()
#legend("topright", inset = c(-0.4, 0), legend=c("HOBO", "iBUT", "PanoIMG", "IRG[ON/OFF]"), col=c("blue", "green" ,"red", "black"), lwd=1)
dev.off()
#########


#################################
#################################################################################################################################
### EC system on shore LAKE GLUBOKOE 2019 - 2020
Irgason<-"/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/eva_wl_net2020/20191207_20200108_EC_FLUX.txt"  # Lake Glubokoe
Irg_LG<-to_read_Irgason(Irgason)
cor_angle <- 36                      # 180 = 144 + 36, where 144 is the sonic direction
Irg_LG$wind_dir_sonic <- (Irg_LG$wind_dir_sonic + cor_angle)%%360                                     #correct for the  wind direction
#################################################################################################################

######################## filtering the fluxes as suggested by Daniela Franz ###############################################
### Flag conc and fluxes: to flag signal strengths lower than 0.7 + wind direction from the land + individually outliers
### and to count the percentage of the filtered data
dat1 <- Irg_LG
all <- nrow(dat1[!is.na(dat1$H2O_conc),])    # to account to the percentage of data filtered due to the signal strengths (lower than 0.7)
filt <- nrow(dat1[(!is.na(dat1$H2O_sig_str) & dat1$H2O_sig_str >= 0.7),])
(all-filt)*100/all # 0 %

filt <- nrow(dat1[(dat1$wind_dir_sonic >= 90 & dat1$wind_dir_sonic <= 225),])            # to wind direction from the land
(all-filt)*100/all  # 11.3 %

# to account the percentage of the data after the filtering
# total after the Latent sensible flux and Evaporation
filt <- nrow(dat1[!is.na(dat1$Evap) & dat1$H2O_sig_str >= 0.7 & (dat1$wind_dir_sonic >= 90 & dat1$wind_dir_sonic <= 225),])
100-((all-filt)*100/all) # 87.5 %

dat1$Evap_filter <- ifelse(dat1$H2O_sig_str >= 0.7 &
                             (dat1$wind_dir_sonic >= 90 & dat1$wind_dir_sonic <= 225),
                           dat1$Evap,
                           NA)

####################################################################################################################################

#####################################
####  daily EC evaporation calculated with NA in the values: used further for the errors estimations:
double1<-dat1$Evap_filter  # by mean
double1[is.na(double1)]<-mean(dat1$Evap_filter, na.rm=T)
dat1$Evap_filter<-double1
rm(double1)

Irg_LG <- dat1
##################################################################################################################
# to aggeragate to 24 hours
evap_ECbyMean<-aggregate(list(Evap=dat1$Evap_filter),list(Timestamp= cut(dat1$Timestamp,"1440 min")), sum, na.rm=T)     ### daily sum
evap_ECbyMean$newTS<-as.Date(evap_ECbyMean$Timestamp, "%Y-%m-%d")
evap_daily_LG <- evap_ECbyMean

# write to the expernal file ! indirect methods!
write.csv(evap_daily_LG, "/home/shevnina/Manuscripts/2025/Ant/september/result_tab/evap_daily_LG.csv") # to write the EC output for the indirect methods
rm(dat1, evap_ECbyMean, cor_angle,filt,all)


##################################################################################################################
## LSWT, Lake Glubokoe
# Function to read and process LSWT data for Lake Glubokoe [by Copilot]
read_process_hobo <- function(filepath, start_date = "2019-12-08", end_date = "2020-02-12") {
  dat <- read.table(filepath, sep = ",", skip = 2, header = FALSE)[, 2:5]
  colnames(dat) <- c("Timestamp", "Abs_Press", "TW", "WL")

  # Correct the date format
  dat$Timestamp <- paste(
    substr(dat$Timestamp, 1, 8),
    str_replace_all(substr(dat$Timestamp, 14, 24), c("ip" = "PM", "ap" = "AM")),
    sep = " "
  )
  dat$Timestamp_EET <- as.POSIXlt(dat$Timestamp, tz = "Europe/Helsinki", format = "%m.%d.%y %I.%M.%S %p")
  dat$Timestamp_UTC <- dat$Timestamp_EET - 2 * 60 * 60  # Convert EET to UTC
  dat$one_min <- strptime(dat$Timestamp_UTC, format = "%Y-%m-%d %H:%M:%S")

  # Daily summary
  hobo24hour <- dat %>%
    group_by(Timestamp = cut(one_min, "24 hour")) %>%
    summarize(
      TW_mean = mean(TW, na.rm = TRUE),
      TW_min = min(TW, na.rm = TRUE),
      TW_max = max(TW, na.rm = TRUE)
    )
  hobo24hour$newTS <- as.Date(hobo24hour$Timestamp, "%Y-%m-%d")
  hobo_day <- subset(hobo24hour, newTS >= start_date & newTS <= end_date)

  list(hobo_day = hobo_day, hobo_data = dat)
}

# Usage
hobo_result <- read_process_hobo("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/eva_wl_net2020/LAKSHMI_season2020.csv")
hobo_Glubokoe <- hobo_result$hobo_day
Hobo_LG <- hobo_result$hobo_data
############################## EOF Copilot

#### to calculate the relative humidity
Irg_LG$RH <- to_calulate_RH(Irg_LG$Temp_amb, Irg_LG$H2O_conc)                           # Hoeltgebaum et al. (2020)
Irg_LG$RH[Irg_LG$RH > 100] <- NA

#### to merge the LSWT and EC
ECH_LG <- merge(Irg_LG, Hobo_LG, by = "Timestamp_UTC")

# the LSWT by HOBO temperature sensor<-ECH_GL$TW
ECH_LG$SVD <- to_calulate_SVD(ECH_LG$Temp_amb, ECH_LG$TW)

# the Supplements
write.csv(ECH_LG,"/home/shevnina/Manuscripts/2025/Ant/september/result_tab/ECH_GL.csv")
# to write the EC output for the indirect methods
BA_LG <- ECH_LG[,c("Timestamp_UTC","RH","Temp_amb","wind_speed","Amb_Press","Evap","Evap_filter","TW")]
write.csv(BA_LG,"/home/shevnina/Manuscripts/2025/Ant/september/result_tab/BA_LG.csv")
####################################################################################################################
# EOF of the EC system on shore LAKE GLUBOKOE 2019 - 2020

##########################################################################################################################
#### Figures
### Figure 3 (revied text) by Copilot + manual correction

# to calculate daily min/mean/max of multiple variables
t24_Lg <- tmp_stat_MultiVars(ECH_LG)
t24_Lg$d20191201 <- as.numeric(t24_Lg$newTS - as.Date("2019-12-01"))

t24_Zb <- tmp_stat_MultiVars(ECH_ZB)
t24_Zb$d20171201 <- as.numeric(t24_Zb$newTS - as.Date("2017-12-01"))
#new Fig. 3
# Plot daily min/mean/max for Temp_amb, RH, and wind_speed for both sites
pdf('/home/shevnina/Manuscripts/2025/Ant/september/result_fig/t24_LgZb_vars_daysince1Dec.pdf')
par(family = "serif", cex = 1.0)
par(mfrow = c(3, 1), mar = c(4, 6, 2, 1))

# Temp_amb
plot(t24_Lg$d20191201, t24_Lg$Temp_amb_mean, type = "l", col = "red", lwd = 0.8,
     ylab = expression(T~(degree~C)), xlab = "Days since 01.12", main = "", ylim = c(-10, 5), xlim = c(1, 95), bty = "n", cex.lab = 1.5)
lines(t24_Lg$d20191201, t24_Lg$Temp_amb_min, col = "red", lty = 2, lwd = 0.5)
lines(t24_Lg$d20191201, t24_Lg$Temp_amb_max, col = "red", lty = 2, lwd = 0.5)
lines(t24_Zb$d20171201, t24_Zb$Temp_amb_mean, col = "pink", lwd = 0.8)
lines(t24_Zb$d20171201, t24_Zb$Temp_amb_min, col = "pink", lty = 2, lwd = 0.5)
lines(t24_Zb$d20171201, t24_Zb$Temp_amb_max, col = "pink", lty = 2, lwd = 0.5)
abline(h = 0, col = "black", lwd = 0.2, lty = 1)
grid()

# Place legend to the right of the first subplot
legend("right", legend = c("LG", "ZB"), col = c("red", "pink"), lwd = 2, bty = "n")

# RH
plot(t24_Lg$d20191201, t24_Lg$RH_mean, type = "l", col = "green", lwd = 0.8,
     ylab = "RH (%)", xlab = "Days since 01.12", main = "", ylim = c(0, 100), xlim = c(1, 95), bty = "n", cex.lab = 1.5)
lines(t24_Lg$d20191201, t24_Lg$RH_min, col = "green", lty = 2, lwd = 0.5)
lines(t24_Lg$d20191201, t24_Lg$RH_max, col = "green", lty = 2, lwd = 0.5)
lines(t24_Zb$d20171201, t24_Zb$RH_mean, col = "darkgreen", lwd = 0.8)
lines(t24_Zb$d20171201, t24_Zb$RH_min, col = "darkgreen", lty = 2, lwd = 0.5)
lines(t24_Zb$d20171201, t24_Zb$RH_max, col = "darkgreen", lty = 2, lwd = 0.5)
abline(h = 0, col = "black", lwd = 0.2, lty = 2)
grid()
legend("right", legend = c("LG", "ZB"), col = c("green", "darkgreen"), lwd = 2, bty = "n", inset = 0.01)


# wind_speed
plot(t24_Lg$d20191201, t24_Lg$wind_speed_mean, type = "l", col = "slateblue2", lwd = 0.8,
     ylab = expression(WS~(ms^{-1})), xlab = "Days since 01.12", main = "", ylim = c(0, 20), xlim = c(1, 95), bty = "n", cex.lab = 1.5)
lines(t24_Lg$d20191201, t24_Lg$wind_speed_min, col = "slateblue2", lty = 2, lwd = 0.5)
lines(t24_Lg$d20191201, t24_Lg$wind_speed_max, col = "slateblue2", lty = 2, lwd = 0.5)
lines(t24_Zb$d20171201, t24_Zb$wind_speed_mean, col = "blue", lwd = 0.8)
lines(t24_Zb$d20171201, t24_Zb$wind_speed_min, col = "blue", lty = 2, lwd = 0.5)
lines(t24_Zb$d20171201, t24_Zb$wind_speed_max, col = "blue", lty = 2, lwd = 0.5)
abline(h = 0, col = "black", lwd = 0.2, lty = 1)
grid()
legend("right", legend = c("LG", "ZB"), col = c("slateblue2", "blue"), lwd = 2, bty = "n")

dev.off()
#############################################################################################

############### FIGURES
## intradaily cycle
## fig. 8 (revised text) Copilot + manual correction
## Prepare data for Temp_amb
df_LG_AT <- ECH_LG %>%
  mutate(TimeOfDay = format(floor_date(Timestamp_UTC, "60 minutes"), "%H:%M:%S"),
         Site = "LG") %>%
  select(TimeOfDay, AT = Temp_amb, Site)
df_ZB_AT <- ECH_ZB %>%
  mutate(TimeOfDay = format(floor_date(Timestamp_UTC, "60 minutes"), "%H:%M:%S"),
         Site = "ZB") %>%
  select(TimeOfDay, AT = Temp_amb, Site)
df_both_AT <- bind_rows(df_LG_AT, df_ZB_AT)

# Prepare data for TW
df_LG_TW <- ECH_LG %>%
  mutate(TimeOfDay = format(floor_date(Timestamp_UTC, "60 minutes"), "%H:%M:%S"),
         Site = "LG") %>%
  select(TimeOfDay, TW, Site)
df_ZB_TW <- ECH_ZB %>%
  mutate(TimeOfDay = format(floor_date(Timestamp_UTC, "60 minutes"), "%H:%M:%S"),
         Site = "ZB") %>%
  select(TimeOfDay, TW, Site)
df_both_TW <- bind_rows(df_LG_TW, df_ZB_TW)

# Plot for Temp_amb
p_AT_hour <- ggplot(df_both_AT, aes(x = TimeOfDay, y = AT, fill = Site)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  labs(x = "Time of Day", y = "T, (\u00B0C)", title = "")+
  scale_y_continuous(breaks = seq(-5, 5, by = 5), limits = c(-5, 5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    legend.position = "top"
  )

# Plot for TW
p_TW_hour <- ggplot(df_both_TW, aes(x = TimeOfDay, y = TW, fill = Site)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  labs(x = "Time of Day", y = "LSWT, (\u00B0C)", title = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = "none"
  )

# Combine plots vertically
combined_plot <- p_AT_hour / p_TW_hour + plot_layout(heights = c(1, 1))
ggsave('/home/shevnina/Manuscripts/2025/Ant/september/result_fig/Intradaily_AT_TW.png', plot = combined_plot, width = 12, height = 12, dpi = 300)
#########################################

# Prepare wind_speed data for LG and ZB
df_LG_WS <- ECH_LG %>%
  mutate(TimeOfDay = format(floor_date(Timestamp_UTC, "60 minutes"), "%H:%M:%S"),
         Site = "LG") %>%
  select(TimeOfDay, wind_speed, Site)
df_ZB_WS <- ECH_ZB %>%
  mutate(TimeOfDay = format(floor_date(Timestamp_UTC, "60 minutes"), "%H:%M:%S"),
         Site = "ZB") %>%
  select(TimeOfDay, wind_speed, Site)
df_both_WS <- bind_rows(df_LG_WS, df_ZB_WS)

# Prepare Evap data for LG and ZB
df_LG_Evap <- ECH_LG %>%
  mutate(TimeOfDay = format(floor_date(Timestamp_UTC, "60 minutes"), "%H:%M:%S"),
         Site = "LG") %>%
  select(TimeOfDay, Evap, Site)
df_ZB_Evap <- ECH_ZB %>%
  mutate(TimeOfDay = format(floor_date(Timestamp_UTC, "60 minutes"), "%H:%M:%S"),
         Site = "ZB") %>%
  select(TimeOfDay, Evap, Site)
df_both_Evap <- bind_rows(df_LG_Evap, df_ZB_Evap)

# Boxplot for wind_speed
p_WS_hour <- ggplot(df_both_WS, aes(x = TimeOfDay, y = wind_speed, fill = Site)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  labs(x = "Time of Day", y = "Wind speed (m/s)", title = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    legend.position = "top"
  )

# Boxplot for Evaporation
p_Evap_hour <- ggplot(df_both_Evap, aes(x = TimeOfDay, y = Evap, fill = Site)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  labs(x = "Time of Day", y = "Evaporation (mm h-1)", title = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = "none"
  )

# Combine plots vertically
library(patchwork)
combined_plot <- p_WS_hour / p_Evap_hour + plot_layout(heights = c(1, 1))
ggsave('/home/shevnina/Manuscripts/2025/Ant/september/result_fig/Intradaily_WS_Evap.png', plot = combined_plot, width = 12, height = 12, dpi = 300)

# Boxplot of intradaily cycle for multiple variables at LG and ZB (Copilot + manual)
# Prepare data for LG
df_LG <- ECH_LG %>%
  mutate(TimeOfDay = format(floor_date(Timestamp_UTC, "60 minutes"), "%H:%M:%S"),
         Site = "LG") %>%
  select(TimeOfDay, AT = Temp_amb, Site)

# Prepare data for ZB
df_ZB <- ECH_ZB %>%
  mutate(TimeOfDay = format(floor_date(Timestamp_UTC, "60 minutes"), "%H:%M:%S"),
         Site = "ZB") %>%
  select(TimeOfDay, AT = Temp_amb, Site)

# Combine
df_both <- bind_rows(df_LG, df_ZB)

# Plot boxplots for AT by site
p_AT_hour <- ggplot(df_both, aes(x = TimeOfDay, y = AT, fill = Site)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  labs(x = "Time of Day", y = "AT", title = "") +
  scale_y_continuous(breaks = seq(-5, 5, by = 5), limits = c(-5, 5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )
ggsave('/home/shevnina/Manuscripts/2025/Ant/september/result_fig/Intradaily_AT.png', plot = p_AT_hour, width = 10, height = 8, dpi = 300)

#### LSWT
# Prepare data for LG
df_LG_TW <- ECH_LG %>%
  mutate(TimeOfDay = format(floor_date(Timestamp_UTC, "60 minutes"), "%H:%M:%S"),
         Site = "LG") %>%
  select(TimeOfDay, TW, Site)

# Prepare data for ZB (if available)
df_ZB_TW <- ECH_ZB %>%
  mutate(TimeOfDay = format(floor_date(Timestamp_UTC, "60 minutes"), "%H:%M:%S"),
         Site = "ZB") %>%
  select(TimeOfDay, TW, Site)

# Combine
df_both_TW <- bind_rows(df_LG_TW, df_ZB_TW)

# Plot boxplots for TW by site
p_TW_hour <- ggplot(df_both_TW, aes(x = TimeOfDay, y = TW, fill = Site)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  labs(x = "Time of Day", y = "TW", title = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

ggsave('/home/shevnina/Manuscripts/2025/Ant/september/result_fig/Intradaily_TW.png', plot = p_TW_hour, width = 10, height = 8, dpi = 300)
###########################################################################
