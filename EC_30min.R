# The supplement to Shevnina, Potes, Vihma and Naakka 2025: Assessing evaporation... 
# Authorship: eshevnina@gmail.com/  mpotes@uevora.pt 
# phone: +358449185446
# December, 2024

library(stringr)
#library(openair)
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
to_read_Irgason<- function(Irgason)  ## the argument is the path to file with the HOBO measurements
{
  variable_unit <- scan(Irgason, nlines = 18, sep=",", what = character())
  dat1 <- read.csv(Irgason, skip = 18, header = FALSE, fill=TRUE)
  dat1 <- dat1[,1:24]
  colnames(dat1) <- c("TS","npoints","u_star","taur","Hcr","LE_wplr", "Evap","Fc_wplr","Fc_wpl_umol","CO2_conc","CO2_PPM","CO2_sig_str",
                      "H2O_conc","H2O_sig_str","Temp_amb","Amb_Press","wind_speed","wind_dir_sonic", "obukhov","zeta","rou","sigmaw_ustar","Xmax","XR90")
  
  dat1$Timestamp_UTC <- as.POSIXct(strptime(substr(dat1$TS,1,15),tz="UTC",format="%y/%m/%d %H:%M") )   # to create proper UTC Timestamp, values corresponded to the end of time interval
  return(dat1)
}


to_read_hobo<- function(Indat4)  ## in function! the argument is the phth to file with the HOBO mearuremets
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
### to calculate the relative humidity from the H2O concentration after Hoeltgebaum et al. (2020)
# wv is H2O concentration,  g/m3 
to_calulate_RH <- function (Tair, wv) {
  R <- 287.05                                                                     ###the universal gas constant for a dry air, J kg-1 K-1
  SatPre <- 0.6113 * exp(17.27*Tair/(Tair + 237.3))                             # Tetens formula in Stull, 2007, p. 89
  SpecHumid <- R*wv*0.001*(Tair + 273.15)/622                    
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
write.csv(ECH_ZB,"/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/ECH_ZB.csv")



## Sub-daily cycle
#### hourly aggregations
ECH_ZB_60m <- ECH_ZB %>%
  group_by(time = floor_date(Timestamp_UTC, "60 minutes")) %>%
  summarize(
    AT=  mean(Temp_amb),
    WS = mean(wind_speed),
    WD = mean(wind_dir_sonic),
    E60 = sum(Evap),
    LSWT = mean(TW))

# !non-filtered data ! 
df <-  ECH_ZB_60m
df <- df %>%
  mutate(TimeOfDay = format(time, "%H:%M:%S"))

# Convert the data to long format for easy plotting of multiple variables
df_long <- df %>%
  pivot_longer(cols = c(AT, WS, LSWT, WD), names_to = "Variable", values_to = "Value")

# Boxplot of intradaily cycle for multiple variables
p_multi <- ggplot(df_long, aes(x = TimeOfDay, y = Value)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  facet_wrap(~Variable, scales = "free_y") +
  labs(x = "Time of Day", y = "Value", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels
ggsave("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/InDay_ZB.png", plot = p_multi, width = 8, height = 6, dpi = 300)

p_E_hour <- ggplot(df, aes(x = TimeOfDay, y = E60)) +
  geom_boxplot(fill = "lightblue", color = "black") +
#  ylim(-0.05, 0.35) +
  labs(x = "Time of Day", y = "Evaporation rate (mm h-1)", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
ggsave("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/InDayEvap_ZB.png", plot = p_E_hour, width = 8, height = 6, dpi = 300)

#####################################################

p_WD <- ggplot(df, aes(x = TimeOfDay, y = WD)) +
  geom_boxplot(fill = "green", color = "black") +
  labs(x = "Time of Day", y = "Wind direction, degree", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
print(p_WD)

p_LSWT <- ggplot(df, aes(x = TimeOfDay, y = LSWT)) +
  geom_boxplot(fill = "blue", color = "black") +
  labs(x = "Time of Day", y = "LSWT (°C)", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
print(p_LSWT)

p_AT <- ggplot(df, aes(x = TimeOfDay, y = AT)) +
  geom_boxplot(fill = "pink", color = "black") +
  labs(x = "Time of Day", y = "Air Temperature (°C)", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
print(p_AT)

p_WS <- ggplot(df, aes(x = TimeOfDay, y = WS)) +
  geom_boxplot(fill = "lightgreen", color = "black") +
  labs(x = "Time of Day", y = "Wind speed (m s-1)", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
print(p_WS)

### EOF of sub-daily cycle
#######################################################

## Figure 3 a
AWTemp_ZB <-tmp_stat_AirWater(ECH_ZB)
# min/mean/max of LWST, Lake ZUB + IRgason Air Temp [Dec, 2017 - Feb,2018]
pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/LSWTbyHobo_TSplot_ZB.pdf")
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
#lines(aws_data$date_new, aws_data$AT, col = "pink", lwd = 1)

abline(v=as.Date("2018-01-01"), col="black",lty=2, lwd=1)                # camera + IRGASON ON
abline(v=as.Date("2018-02-07"), col="black",lty=2, lwd=1)                # Irgason OFF

abline(h=0, col="black", lwd=0.5) 
grid()
#legend("topright", inset = c(-0.4, 0), legend=c("HOBO", "iBUT", "PanoIMG", "IRG[ON/OFF]"), col=c("blue", "green" ,"red", "black"), lwd=1)
dev.off()


###################################################
HZ_daily_stat <-tmp_stat_HcLe_flux(ECH_ZB)
pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/Hcr_Le_TSplot_ZB.pdf")
par(family = "serif", cex = 1.0)
par(xpd = TRUE)
#par(mar = par()$mar + c(0, 0, 0, 7))

plot(HZ_daily_stat$newTS, HZ_daily_stat$Hcr_mean, type="l", xlim=c(as.Date("2017-12-01"),as.Date("2018-02-28")), ylab= "H, Wm-2", xlab="Date", ylim=c(-110,290),col="orange",lwd=1)
lines(HZ_daily_stat$newTS, HZ_daily_stat$Hc_max, col="orange", lwd=0.5, lty=2)
lines(HZ_daily_stat$newTS, HZ_daily_stat$Hc_min, col="orange",lwd=0.5, lty=2)
abline(h=0, col="black", lwd=0.5) 
grid()
dev.off()

pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/LWSTvsEvap.pdf")
 plot(inBA_LG$TW, inBA_LG$Evap, col="black", pch=20, type="p", ylim=c(-0.05,0.2), xlim=c(0, 20), ylab="E, mm/30min", xlab="LWST, C")
 points(inBA_ZB$TW, inBA_ZB$Evap, col="red", pch=1, cex = 0.5)
 abline(0,1, col="grey", lty=2) 
 #abline(1,0, col="grey", lty=2) 
 grid()
dev.off()
#######################################################################################################################


#################################################################################################################################
### EC system on shore LAKE GLUBOKOE
Irgason<-"/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/eva_wl_net2020/20191207_20200108_ANTARTIDA_FLUX.txt"  # Lake Glubokoe
Irg_LG<-to_read_Irgason(Irgason) 
cor_angle <- 36   # 180 = 144 + 36, where 144 is the sonic direction
Irg_LG$wind_dir_sonic <- (Irg_LG$wind_dir_sonic + cor_angle)%%360                                     #correct for the  wind direction
#################################################################################################################
## LSWT 
######### to read and calculate the daily LSWT from HOBO sensors: season 2019-2020
setwd("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/eva_wl_net2020/") 
GL<-to_read_hobo("LAKSHMI_season2020.csv")                                      # Glubokoe, daily 
#PM<-to_read_hobo("DURGA_season2020.csv")
#VR<-to_read_hobo("PARVATI_season2020.csv")

hobo_Glubokoe<-subset(GL, GL$newTS>="2019-12-08" & GL$newTS<="2020-02-17")      #sub setting for the period 8 Dec,2019 - 17.02.2020 Hobo operation
#hobo_Verhnee<-subset(VR, VR$newTS>="2019-12-08" & VR$newTS<="2020-02-17")
#hobo_Pomornik<-subset(PM, PM$newTS>="2019-12-08" & PM$newTS<="2020-02-17")

# LWST, Lake Glubokoe
Indat4 <- "/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/eva_wl_net2020/LAKSHMI_season2020.csv"        
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
Hobo_LG <- dat4   
rm(dat4)
rm(Indat4)
#######################################################################################################################
#### to calculate the relative humidity
Irg_LG$RH <- to_calulate_RH(Irg_LG$Temp_amb, Irg_LG$H2O_conc)                           # Hoeltgebaum et al. (2020)
#Irg_LG$RH[Irg_LG$RH > 100] <- NA

##summary(Irg_LG$RH)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  31.18   49.60   56.73   57.27   64.32   89.85      13 

#### to merge the LSWT and EC
ECH_LG <- merge(Irg_LG, Hobo_LG, by = "Timestamp_UTC")
# the Supplement 1B
write.csv(ECH_LG,"/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/ECH_GL.csv")


#############################################################################################
## Sub-daily cycle
#### hourly aggregation
ECH_GL_60m <- ECH_GL %>%
  group_by(time = floor_date(Timestamp_UTC, "60 minutes")) %>%
  summarize(
    AT=  mean(Temp_amb),
    WS = mean(wind_speed),
    WD = mean(wind_dir_sonic),
    E60 = sum(Evap),
    LSWT = mean(TW))

# !non-filtered data
#####################################################
df<-ECH_GL_60m 
# Extract time of day
df <- df %>%
  mutate(TimeOfDay = format(time, "%H:%M:%S"))

# Convert the data to long format for easy plotting of multiple variables
df_long <- df %>%
  pivot_longer(cols = c(AT, WS, LSWT, WD), names_to = "Variable", values_to = "Value")

# Boxplot of intradaily cycle for multiple variables
p_multi <- ggplot(df_long, aes(x = TimeOfDay, y = Value)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  facet_wrap(~Variable, scales = "free_y") +
  labs(x = "Time of Day", y = "Value", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels
ggsave("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/InDay_GL.png", plot = p_multi, width = 8, height = 6, dpi = 300)

p_E_hour <- ggplot(df, aes(x = TimeOfDay, y = E60)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(x = "Time of Day", y = "Evaporation rate (mm h-1)", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
ggsave("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/InDayEvap_GL.png", plot = p_E_hour, width = 8, height = 6, dpi = 300)

#+++++++
p_WD <- ggplot(df, aes(x = TimeOfDay, y = WD)) +
  geom_boxplot(fill = "green", color = "black") +
  labs(x = "Time of Day", y = "Wind direction, degree", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
print(p_WD)


p_LSWT <- ggplot(df, aes(x = TimeOfDay, y = LSWT)) +
  geom_boxplot(fill = "blue", color = "black") +
  labs(x = "Time of Day", y = "LSWT (°C)", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
print(p_LSWT)

p_AT <- ggplot(df, aes(x = TimeOfDay, y = AT)) +
  geom_boxplot(fill = "pink", color = "black") +
  labs(x = "Time of Day", y = "Air Temperature (°C)", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
print(p_AT)

p_E_hour <- ggplot(df, aes(x = TimeOfDay, y = E60)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(x = "Time of Day", y = "Evaporation rate (mm h-1)", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
print(p_E_hour)

p_WS <- ggplot(df, aes(x = TimeOfDay, y = WS)) +
  geom_boxplot(fill = "lightgreen", color = "black") +
  labs(x = "Time of Day", y = "Wind speed (m s-1)", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
print(p_WS)
### EOF of sub-daily cycle
#######################################################


#####################################
HL_daily_stat <-tmp_stat_HcLe_flux(ECH_LG)
pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/Hcr_LE_TSplot_ZB.pdf")
par(family = "serif", cex = 1.0)
par(xpd = TRUE)
#par(mar = par()$mar + c(0, 0, 0, 7))

plot(HL_daily_stat$newTS, HL_daily_stat$Hcr_mean, type="l", xlim=c(as.Date("2019-12-01"),as.Date("2020-02-28")), ylab= "Flux30m, Wm-3", xlab="Date", ylim=c(-110,290),col = "orange",lwd=1)
lines(HL_daily_stat$newTS, HL_daily_stat$Hc_max, col="orange", lwd=0.5, lty=2)
lines(HL_daily_stat$newTS, HL_daily_stat$Hc_min, col="orange",lwd=0.5, lty=2)

#par(new=TRUE)
#lines(HL_daily_stat$newTS, HL_daily_stat$Le_mean, type="l", xlim=c(as.Date("2019-12-01"),as.Date("2020-02-28")), col="green", lwd=1)
#lines(HL_daily_stat$newTS, HL_daily_stat$Le_max, col="red", lwd=0.5, lty=2)
#lines(HL_daily_stat$newTS, HL_daily_stat$Le_min, col="red",lwd=0.5, lty=2)

abline(v=as.Date("2019-12-07"), col="black",lty=2, lwd=1)                # camera + IRGASON ON
#abline(v=as.Date("2019-12-08"), col="red",lty=2, lwd=0.5)                # camera + IRGASON ON
abline(v=as.Date("2020-01-03"), col="red",lty=2, lwd=0.5)                # camera, panoramic photo
abline(v=as.Date("2020-02-25"), col="red",lty=2, lwd=0.5)                # camera, panoramic photo
abline(v=as.Date("2020-01-08"), col="black",lty=2, lwd=1)                # Irgason OFF

abline(h=0, col="black", lwd=0.5) 
grid()
#legend("topright", inset = c(-0.4, 0), legend=c("HOBO", "iBUT", "PanoIMG", "IRG[ON/OFF]"), col=c("blue", "green" ,"red", "black"), lwd=1)
dev.off()


# Figure 3 b
AWTemp_GL <-tmp_stat_AirWater(ECH_LG)
# min/mean/max of LWST, Lake Glubokoe + IRgason Air Temp [Dec, 2019 - Feb,2020]
pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/LSWTbyHobo_TSplot_GL.pdf")
par(family = "serif", cex = 1.0)
par(xpd = TRUE)
#par(mar = par()$mar + c(0, 0, 0, 7))

plot(AWTemp_GL$newTS, AWTemp_GL$TW_mean, type="l", xlim=c(as.Date("2019-12-01"),as.Date("2020-02-28")), ylab= "T, C", xlab="Date", ylim=c(-10,10),col="blue",lwd=2)
lines(AWTemp_GL$newTS, AWTemp_GL$TW_max, col="blue", lwd=0.5, lty=2)
lines(AWTemp_GL$newTS, AWTemp_GL$TW_min, col="blue",lwd=0.5, lty=2)

#lines(Hobo_LG$newTS, Hobo_LG$TW_mean, col="blue",lwd=2)
#lines(Hobo_LG$newTS, Hobo_LG$TW_max, col="blue", lwd=0.5, lty=2)
#lines(Hobo_LG$newTS, Hobo_LG$TW_min, col="blue",lwd=0.5, lty=2)
lines(hobo_Glubokoe$newTS, hobo_Glubokoe$TW_mean, col="blue", lwd=2)
lines(hobo_Glubokoe$newTS, hobo_Glubokoe$TW_max, col="blue", lwd=0.5, lty=2)
lines(hobo_Glubokoe$newTS, hobo_Glubokoe$TW_min, col="blue",lwd=0.5, lty=2)

#par(new=TRUE)
lines(AWTemp_GL$newTS, AWTemp_GL$AT_mean, type="l", xlim=c(as.Date("2019-12-01"),as.Date("2020-02-28")), col="red",lwd=2)
lines(AWTemp_GL$newTS, AWTemp_GL$AT_max, col="red", lwd=0.5, lty=2)
lines(AWTemp_GL$newTS, AWTemp_GL$AT_min, col="red",lwd=0.5, lty=2)
#lines(aws_data$newTS, aws_data$AT, col = "pink", lwd = 1)

abline(v=as.Date("2019-12-07"), col="black",lty=2, lwd=1)                # camera + IRGASON ON
#abline(v=as.Date("2019-12-08"), col="red",lty=2, lwd=0.5)                # camera + IRGASON ON
abline(v=as.Date("2020-01-03"), col="red",lty=2, lwd=0.5)                # camera, panoramic photo
abline(v=as.Date("2020-02-25"), col="red",lty=2, lwd=0.5)                # camera, panoramic photo
abline(v=as.Date("2020-01-08"), col="black",lty=2, lwd=1)                # Irgason OFF

abline(h=0, col="black", lwd=0.5) 
grid()
#legend("topright", inset = c(-0.4, 0), legend=c("HOBO", "iBUT", "PanoIMG", "IRG[ON/OFF]"), col=c("blue", "green" ,"red", "black"), lwd=1)
dev.off()

#lines(hobo_Verhnee$newTS, hobo_Verhnee$TW_mean, col="red", lwd=1)
#lines(hobo_Pomornik$newTS, hobo_Pomornik$TW_mean, col="green",lwd=1)

#################################################################################################################################################
# Novo AWS 1 minute measurements AT, RH, WS, WD, P, soil temp, incoming solar radiation ... 
aws_data_1min <- read.csv("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/meteo_sites/20191201_20200131_Novo.txt", header = TRUE, sep = ";")
aws_data_1min$Timestamp_UTC = as.POSIXct(aws_data_1min$Dates_Times, format="%Y-%m-%d %H:%M")

p <- ggplot(aws_data_1min, aes(x = Timestamp_UTC, y = Relative_Humidity)) +
       geom_line(size=0.5, color="blue") +
       labs(x = "Date", y = "RH (%)", title = "") 
ggsave("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/RelativeHumidity_GL.png", plot = p, width = 8, height = 4, dpi = 300)

specific_time <- as.POSIXct("2019-12-25 23:30:00", format = "%Y-%m-%d %H:%M:%S")
p <- ggplot(aws_data_1min, aes(x = Timestamp_UTC, y = Relative_Humidity)) +
  geom_line(size=0.1, color="blue") +                                 # Line for Relative_Humidity from aws_data_1min
  geom_line(data = ECH_LG, aes(x = Timestamp_UTC, y = RH),                # Line for RH from df
            size=0.2, color="red") +                                  # Customize size and color for the second variable
  geom_vline(xintercept = specific_time, linetype = "dashed",        # Vertical line at the specified time
             color = "black", size = 0.4) +  
  labs(x = "Date", y = "RH (%)", title = "") +
  ylim(0, 100) +
  theme_minimal()
ggsave("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/RelativeHumidity_GL.png", plot = p, width = 8, height = 4, dpi = 300)


# to 30m average
aws_data_30m <- aws_data_1min %>%
  group_by(time = floor_date(Timestamp_UTC, "30 minutes")) %>%
  summarize(
    AT = mean(Air_temperature, na.rm=TRUE),
    ST = mean(Soil_Temperature, na.rm=TRUE),
    RH = mean(Relative_Humidity, na.rm=TRUE),
    ISR = mean(Incoming_solar_radiation, na.rm=TRUE) )


# subset for the period of the experiment
aws_data_30m_38d <- subset(aws_data_30m, aws_data_30m$time >= min(Irg_LG$Timestamp_UTC,na.rm=TRUE) & aws_data_30m$time <= max(Irg_LG$Timestamp_UTC,na.rm=TRUE))   ## sub setting for the EC data

# to 1hour average
aws_data_60m <- aws_data_1min %>%
  group_by(time = floor_date(Timestamp_UTC, "60 minutes")) %>%
  summarize(
    AT = mean(Air_temperature, na.rm=TRUE),
    ST = mean(Soil_Temperature, na.rm=TRUE),
    RH = mean(Relative_Humidity, na.rm=TRUE),
    ISR = mean(Incoming_solar_radiation, na.rm=TRUE)
  )

aws_data_60m_38d <- subset(aws_data_60m, aws_data_60m$time >= min(Irg_LG$Timestamp_UTC,na.rm=TRUE) & aws_data_60m$time <= max(Irg_LG$Timestamp_UTC,na.rm=TRUE))   ## sub setting for the EC data

p <- ggplot(aws_data_60m_38d, aes(x = time, y = RH)) +
  geom_line(size=0.5, color="blue") +
  labs(x = "Date", y = "RH (%)", title = "") 
ggsave("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/RelativeHumidity_GL.png", plot = p, width = 8, height = 4, dpi = 300)

p <- ggplot(aws_data_60m_38d, aes(x = time, y = RH)) +
  geom_line(size=0.5, color="blue") +
  labs(x = "Date", y = "RH (%)", title = "") 
ggsave("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/RelativeHumidity_GL.png", plot = p, width = 8, height = 4, dpi = 300)


p <- ggplot(aws_data_60m_38d, aes(x = time, y = ISR)) +
  geom_line(size=0.5, color="orange") +
  labs(x = "Date", y = "Incoming solar radiation (W m-2)", title = "")  
ggsave("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/IncomingSolarRad_GL.png", plot = p, width = 8, height = 4, dpi = 300)


## NOT USED!!!
p <- ggplot(aws_data_30m_38d, aes(x = time)) +
  geom_line(aes(y = AT, color = "AT")) +                  # First variable
  geom_line(aes(y = ST, color = "ST")) +               # Scaled second variable
#  scale_y_continuous(
#   name = "Air Temp, (°C)",                              # Left y-axis label
#    sec.axis = sec_axis(~.*1, name = "Soil Temp, (°C)")  # Right y-axis label with reverse scaling
# ) +
# scale_color_manual(values = c("AT" = "blue", "ST" = "red")) +
  labs(x = "Date", color = "Variables") +
#  theme_minimal()

# Print the plot
print(p)


aws_data_30m_38d$Timestamp_UTC <- aws_data_30m_38d$time
RH_LG <- merge(Irg_LG, aws_data_30m_38d, by = "Timestamp_UTC")
RH_LG <- RH_LG[,c("Timestamp_UTC","RH.x", "RH.y")]

pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/RH_IRGvsNovo_GL.pdf")
par(family = "serif", cex = 1.0)
par(xpd = TRUE)
par(mar = par()$mar + c(0, 0, 0, 7))
plot(RH_LG$Timestamp_UTC, RH_LG$RH, type="l", xlab="Mounth / Day", ylab="RH, %", ylim=c(10,100),col="black",lwd=0.5)
#lines(RH_LG$Timestamp_UTC, RH_LG$RH1, col="blue", lwd=1, lty=2)
lines(RH_LG$Timestamp_UTC, RH_LG$Relative_Humidity, col="red",lwd=0.5, lty=2)
#axis.POSIXct(1, at = seq(min(aws_data_30m$Timestamp, na.rm=TRUE), max(aws_data_30m$Timestamp, na.rm=TRUE), "15 day"), labels = TRUE, format="%m/%d")
grid()
#legend("topright", inset = c(-0.4, 0), legend=c("HOBO", "iBUT", "PanoIMG", "IRG[ON/OFF]"), col=c("blue", "green" ,"red", "black"), lwd=1)
dev.off()


#### figure 5 in Shevnina et al., 2022 ##########################################################################
#### for each method one graph##########################################################

pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/RH_IRGvsNovo_QQ_GL.pdf")
plot(RH_LG$RH.y, RH_LG$RH.x, col="black", pch=20, type="p", xlim=c(10,100), ylim=c(10, 100), xlab="RH_Novo, %", ylab="RH_EC, %")
#points(RH_LG$RH.y, RH_LG$RH1, col="red", pch=0)
abline(0,1, col="grey", lty=1) 

#abline(lm(RH_LG$RH.y ~ RH_LG$RH.x), col="black")
#abline(lm(RH_LG$RH.y ~ RH_LG$RH1), col="red")
grid()

m1<-summary(lm(RH_LG$RH.y ~ RH_LG$RH.x))
#m2<-summary(lm(RH_LG$RH.y ~ RH_LG$RH1))

t<-c("R2 =")
r1<-as.character(round(m1$r.squared, digits=2))
#r2<-as.character(round(m2$r.squared, digits=2))

text(1,2,paste(t,r1), col="black")
#text(1,4,paste(t,r2), col="red")
dev.off()
########################################################################################

#############################################################################################

###  LSWT = AirTmp against Evaporation rate
pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/LSWTvsEvap_LG.pdf")
plot(inBA_LG$TW, inBA_LG$Evap, col="black", pch=20, type="p", ylim=c(-0.05,0.2), xlim=c(-5, 15), ylab="E, mm/30min", xlab=expression(paste("LWST, ",degree,"C")))
#points(inBA_ZB$dTMP, inBA_ZB$Evap, col="red", pch=1, cex = 0.5)
abline(0,1, col="grey", lty=2) 
abline(1,0, col="grey", lty=2) 
grid()
dev.off()

pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/LSWTvsEvap_ZB.pdf")
plot(inBA_ZB$TW, inBA_ZB$Evap, col="red", pch=20, type="p", ylim=c(-0.05,0.2), xlim=c(-5, 15), ylab="E, mm/30min", xlab=expression(paste("LWST, ",degree,"C")))
#points(inBA_ZB$dTMP, inBA_ZB$Evap, col="red", pch=1, cex = 0.5)
abline(0,1, col="grey", lty=2) 
abline(1,0, col="grey", lty=2) 
grid()
dev.off()

# WS against Evap
pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/WindSpeedvsEvap_LG.pdf")
plot(inBA_LG$wind_speed, inBA_LG$Evap, col="black", pch=20, type="p", ylim=c(-0.05,0.2), xlim=c(0, 20), ylab="E, mm/30min", xlab="WS, m/s")
#points(inBA_ZB$wind_speed, inBA_ZB$Evap, col="red", pch=1, cex = 0.5)
abline(0,1, col="grey", lty=2) 
#abline(1,0, col="grey", lty=2) 
grid()
dev.off()

pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/WindSpeedvsEvap_ZB.pdf")
plot(inBA_ZB$wind_speed, inBA_ZB$Evap, col="red", pch=20, type="p", ylim=c(-0.05,0.2), xlim=c(0, 20), ylab="E, mm/30min", xlab="WS, m/s")
abline(0,1, col="grey", lty=2) 
grid()
dev.off()


