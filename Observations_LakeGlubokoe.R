# The supplement to Shevnina, Potes, Vihma and Naakka 2025
# Authorship: eshevnina@gmail.com
# phone: +358449185446
# December, 2024

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
# P is atmospheric pressure, Pa
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
#Figure 5
par(mgp = c(2, 0.5, 0))  # Adjust the margins
plot(dat1$Timestamp_UTC, dat1$Evap, type="p", pch=20, col="red", ylab=expression(Evaporation~(L~m^-2)), xlab="2019 - 2020", ylim=c(0,0.1))
points(dat1$Timestamp_UTC, dat1$Evap_filter, col="black", pch=20)
#abline(h=0, col="black", lwd=0.5)
grid()


###########################################################################
######################################
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
write.csv(evap_daily_LG, "/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/result_tab/evap_daily_LG.csv") # to write the EC output for the indirect methods
rm(dat1, evap_ECbyMean, cor_angle,filt,all)

##################################################################################################################
## LSWT
# Lake Glubokoe
dat4 <- read.table("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/eva_wl_net2020/LAKSHMI_season2020.csv", sep=",", skip=2, header=FALSE)[, 2:5]
colnames(dat4) <- c("Timestamp", "Abs_Press", "TW", "WL")

# Correct the date format
dat4$Timestamp <- paste(substr(dat4$Timestamp, 1, 8),
                        str_replace_all(substr(dat4$Timestamp, 14, 24), c("ip" = "PM", "ap" = "AM")),
                        sep=" ")

dat4$Timestamp_EET <- as.POSIXlt(dat4$Timestamp, tz="Europe/Helsinki", format="%m.%d.%y %I.%M.%S %p")
dat4$Timestamp_UTC <- dat4$Timestamp_EET - 2*60*60  # Convert EET to UTC

Hobo_LG <- dat4

######### to calculate the daily mean of water temperature
## Copilot

dat4 <- Hobo_LG
dat4$one_min <- strptime(dat4$Timestamp_UTC, format="%Y-%m-%d %H:%M:%S")

hobo24hour <- dat4 %>%
  group_by(Timestamp = cut(one_min, "24 hour")) %>%
  summarize(
    TW_mean = mean(TW),
    TW_min = min(TW),
    TW_max = max(TW)
  )

hobo24hour$newTS <- as.Date(hobo24hour$Timestamp, "%Y-%m-%d")
hobo_day <- subset(hobo24hour, newTS >= "2019-12-08" & newTS <= "2020-02-12")

hobo_Glubokoe <- hobo_day                                             # LWST, Lake Gluboke 07.12.2019 - 28.02.2020
Hobo_LG <- dat4                                                       # to merge further with the Irgason observations
rm(hobo24hour, dat4, hobo_day)
#######################################################################################################################################

#######################################################################################################################
#### to calculate the relative humidity
Irg_LG$RH <- to_calulate_RH(Irg_LG$Temp_amb, Irg_LG$H2O_conc)                           # Hoeltgebaum et al. (2020)
#Irg_LG$RH[Irg_LG$RH > 100] <- NA
##summary(Irg_LG$RH)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  31.18   49.60   56.73   57.27   64.32   89.85      13

#### to merge the LSWT and EC
ECH_LG <- merge(Irg_LG, Hobo_LG, by = "Timestamp_UTC")


# the Supplements
write.csv(ECH_LG,"/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/result_tab/ECH_GL.csv")

# to write the EC output for the indirect methods
BA_LG <- ECH_LG[,c("Timestamp_UTC","RH","Temp_amb","wind_speed","Amb_Press","Evap","Evap_filter","TW")]
write.csv(BA_LG,"/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/result_tab/BA_LG.csv")

###### EOF of Observations
##########################################################################################################################


#############################################################################################
## Sub-daily cycle
#### hourly aggregation
ECH_GL_60m <- ECH_LG %>%
  group_by(time = floor_date(Timestamp_UTC, "60 minutes")) %>%
  summarize(
    AT=  mean(Temp_amb, na.rm = TRUE),
    WS = mean(wind_speed, na.rm = TRUE),
    WD = mean(wind_dir_sonic, na.rm = TRUE),
    E60 = sum(Evap, na.rm = TRUE),
    LSWT = mean(TW), na.rm = TRUE)

# non-filtered data
#####################################################
df<-ECH_GL_60m 
# Extract time of day
df <- df %>%
  mutate(TimeOfDay = format(time, "%H:%M:%S"))

# Convert the data to long format for easy plotting of multiple variables
df_long <- df %>%
  pivot_longer(cols = c(AT, WS, LSWT, WD), names_to = "Variable", values_to = "Value")

##################################################### FIGURES, fig. 5 in Shevnina et al., 2025
# Boxplot of intradaily cycle for multiple variables
p_E_hour <- ggplot(df, aes(x = TimeOfDay, y = E60)) +
  geom_boxplot(fill = "white", color = "black") +
  labs(x = "Time of Day", y = "Evaporation rate (mm h-1)", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for better readability
ggsave("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/result_fig/InDayEvap_GL.png", plot = p_E_hour, width = 8, height = 6, dpi = 300)

rm(df)
rm(df_long)
rm(ECH_GL_60m)
#######################################################


##################################### FIGURES, fig. 4 in Shevnina et al., 2025
HL_daily_stat <-tmp_stat_HcLe_flux(ECH_LG)
pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/result_fig/Hcr_LE_TSplot_GL.pdf")
par(family = "serif", cex = 1.0)
par(xpd = TRUE)
par(mar = par()$mar + c(0, 0, 0, 7))

plot(HL_daily_stat$newTS, HL_daily_stat$Hcr_mean, type="l", xlim=c(as.Date("2019-12-01"),as.Date("2020-02-28")), ylab= "Flux30m, Wm-3", xlab="Date", ylim=c(-110,290),col = "orange",lwd=1)
lines(HL_daily_stat$newTS, HL_daily_stat$Hc_max, col="orange", lwd=0.5, lty=2)
lines(HL_daily_stat$newTS, HL_daily_stat$Hc_min, col="orange",lwd=0.5, lty=2)

par(new=TRUE)
lines(HL_daily_stat$newTS, HL_daily_stat$Le_mean, type="l", xlim=c(as.Date("2019-12-01"),as.Date("2020-02-28")), col="green", lwd=1)
lines(HL_daily_stat$newTS, HL_daily_stat$Le_max, col="green", lwd=0.5, lty=2)
lines(HL_daily_stat$newTS, HL_daily_stat$Le_min, col="green",lwd=0.5, lty=2)

abline(v=as.Date("2019-12-07"), col="black",lty=2, lwd=1)                # camera + IRGASON ON
#abline(v=as.Date("2020-01-03"), col="red",lty=2, lwd=0.5)                # camera, panoramic photo
#abline(v=as.Date("2020-02-25"), col="red",lty=2, lwd=0.5)                # camera, panoramic photo
abline(v=as.Date("2020-01-08"), col="black",lty=2, lwd=1)                # Irgason OFF

abline(h=0, col="black", lwd=0.5) 
grid()
#legend("topright", inset = c(-0.4, 0), legend=c("HOBO", "iBUT", "PanoIMG", "IRG[ON/OFF]"), col=c("blue", "green" ,"red", "black"), lwd=1)
dev.off()


#################################################################################################################################################
# Novo AWS 1 minute measurements AT, RH, WS, WD, P, soil temp, incoming solar radiation ... 
aws_data_1min <- read.csv("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/meteo_sites/20191201_20200131_Novo.txt", header = TRUE, sep = ";")
aws_data_1min$Timestamp_UTC <- as.POSIXct(aws_data_1min$Dates_Times, format="%Y-%m-%d %H:%M")

##### to 30m average
aws_data_30m <- aws_data_1min %>%
  group_by(time = floor_date(Timestamp_UTC, "30 minutes")) %>%
  summarize(
    AT = mean(Air_temperature, na.rm=TRUE),
    ST = mean(Soil_Temperature, na.rm=TRUE),
    RH = mean(Relative_Humidity, na.rm=TRUE),
    ISR = mean(Incoming_solar_radiation, na.rm=TRUE) )

# subset for the period of the experiment
aws_data_30m_38d <- subset(aws_data_30m, aws_data_30m$time >= min(Irg_LG$Timestamp_UTC,na.rm=TRUE) & aws_data_30m$time <= max(Irg_LG$Timestamp_UTC,na.rm=TRUE))   ## sub setting for the EC data

#### to 1hour average
aws_data_60m <- aws_data_1min %>%
  group_by(time = floor_date(Timestamp_UTC, "60 minutes")) %>%
  summarize(
    AT = mean(Air_temperature, na.rm=TRUE),
    ST = mean(Soil_Temperature, na.rm=TRUE),
    RH = mean(Relative_Humidity, na.rm=TRUE),
    ISR = mean(Incoming_solar_radiation, na.rm=TRUE)
  )

aws_data_60m_38d <- subset(aws_data_60m, aws_data_60m$time >= min(ECH_LG$Timestamp_UTC,na.rm=TRUE) & aws_data_60m$time <= max(ECH_LG$Timestamp_UTC,na.rm=TRUE))   ## sub setting for the EC data



##### LSWT, AT observations
AWTemp_GL <-tmp_stat_AirWater(ECH_LG)   # statistic of the LSWT and AT for the period of the experiment
Hobo_LG <- hobo_Glubokoe   # statistic of the LSWT for the period 07.12.2019 - 15.02.2020
# min/mean/max of LWST, Lake Glubokoe + IRgason Air Temp [Dec, 2019 - Feb,2020]

#####  AWS Novo, 24 hours average
aws_data24h <- aws_data_1min %>%
  group_by(time = floor_date(Timestamp_UTC, "24 hours")) %>%
  summarize(
    AT = mean(Air_temperature, na.rm=TRUE),
    ST = mean(Soil_Temperature, na.rm=TRUE),
    RH = mean(Relative_Humidity, na.rm=TRUE),
    ISR = mean(Incoming_solar_radiation, na.rm=TRUE),
    WS = mean(Wind_speed, na.rm=TRUE)
  )

aws_data24h$newTS <- as.Date(aws_data24h$time, "%Y-%m-%d") #
##################### EOF NOVO AWS ############################################################################################


################################## FIGURES, Fig. 3 in Shevnina et al., 2025
# to plot the LSWT, AT, AWS data

pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/result_fig/LSWTbyHobo_TSplot_GL.pdf")
par(family = "serif", cex = 1.0)
par(xpd = TRUE)

plot(AWTemp_GL$newTS, AWTemp_GL$TW_mean, type="l", xlim=c(as.Date("2019-12-01"), as.Date("2020-02-28")), ylab="T, C", xlab="Date", ylim=c(-10, 10), col="blue", lwd=2)
lines(AWTemp_GL$newTS, AWTemp_GL$TW_max, col="blue", lwd=0.5, lty=2)
lines(AWTemp_GL$newTS, AWTemp_GL$TW_min, col="blue", lwd=0.5, lty=2)

lines(AWTemp_GL$newTS, AWTemp_GL$AT_mean, col="red", lwd=2)
lines(AWTemp_GL$newTS, AWTemp_GL$AT_max, col="red", lwd=0.5, lty=2)
lines(AWTemp_GL$newTS, AWTemp_GL$AT_min, col="red", lwd=0.5, lty=2)

lines(aws_data24h$newTS, aws_data24h$AT, col="red", lwd=0.5)

lines(Hobo_LG$newTS, Hobo_LG$TW_mean, col="blue", lwd=2)
lines(Hobo_LG$newTS, Hobo_LG$TW_max, col="blue", lwd=0.5, lty=2)
lines(Hobo_LG$newTS, Hobo_LG$TW_min, col="blue", lwd=0.5, lty=2)

abline(v=as.Date("2019-12-07"), col="black", lty=2, lwd=1)               # camera + IRGASON ON

#abline(v=as.Date("2019-12-18"), col="red", lty=2, lwd=0.5)               # image ice cover
abline(v=as.Date("2019-12-25"), col="grey", lty=2, lwd=0.5)               # image fog !!

abline(v=as.Date("2020-01-03"), col="red", lty=2, lwd=0.5)               # image ice cover
abline(v=as.Date("2020-01-14"), col="red", lty=2, lwd=0.5)               # image ice cover
#abline(v=as.Date("2020-01-24"), col="red", lty=2, lwd=0.5)               # image ice cover
abline(v=as.Date("2020-01-08"), col="black", lty=2, lwd=1)               # IrGASON OFF
#abline(v=as.Date("2020-02-15"), col="red", lty=2, lwd=0.5)               # image ice cover
abline(v=as.Date("2020-02-25"), col="red", lty=2, lwd=0.5)               # image ice cover

abline(h=0, col="black", lwd=0.5)
grid()
#legend("topright", inset = c(-0.4, 0), legend=c("HOBO", "iBUT", "PanoIMG", "IRG[ON/OFF]"), col=c("blue", "green" ,"red", "black"), lwd=1)
dev.off()

###
rm(aws_data_30m,aws_data24h,aws_data_1min,aws_data_60m_38d)


##################################################
library(openair)

# Generate wind rose plot
p <- windRose(Irg_LG, ws = "wind_speed", wd = "wind_dir_sonic",
          paddle = FALSE, key.position = "bottom",
          angle=20, grid.line = 15,
          par.settings = list(fontsize=list(text=14)),
          breaks = c(0, 2, 4, 6, 8, 10, 13),
          col = "RdYlBu")


ggsave(p, "/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/result_fig/WindRose_GL.png", width = 8, height = 6, dpi = 300)

