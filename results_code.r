# Supplement to Shevnina et al., 2021: Evaporation over glacial lakes
# Authorship: eshevnina@gmail.com
# phone: +358449185446

library(lubridate)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting

#################################### DATA: Observations  ##################################################
#############################################################################################################################################

#################### Meteoparameters according to Novo site: Contribution by RAE  ###########################################################
################### daily aggregation #######################################################################################################
path<-"/home/shevnina/Manuscripts/2020/Antarctica/evaporation/supplement"
setwd(path)
datOne<-read.table("20171201_20180210_Novo.txt", sep=";", header=T) # meteo DJF 2017/2018
datOne$one_min<-strptime(datOne$Dates_Times, format="%Y-%m-%d %H:%M")
dat24hour<- aggregate(list(AT = datOne$Air_temperature, ST = datOne$Soil_Temperature, RH = datOne$Relative_Humidity, WS = datOne$Wind_speed, NSOL = datOne$Incoming_solar_radiation),list(Timestamp= cut(datOne$one_min, "24 hour" )), mean)
dat24hour$newTS = as.Date(dat24hour$Timestamp, "%Y-%m-%d")
novo1718daily<- dat24hour         # 2017/2018
#  Timestamp         AT       ST       RH       WS     NSOL      newTS
#1 2017-12-01  2.6711111       NA 41.85069 5.307500 376.5026 2017-12-01
#2 2017-12-02  0.7939583       NA 37.02500 3.787847 384.4626 2017-12-02

############## Meteoparameters according to Maitri site: Contribution from P.Dhote and P. Thakur (India) ##############################################
path<-"/home/shevnina/Manuscripts/2021/AntarcticEVA/supplement/evaporation/"
setwd(path)
meteo_m<-read.csv("MeteoMaitri.csv", header = T, sep=",") 
meteo_m$TS<-as.POSIXlt(meteo_m$Date, format = "%m/%d/%Y")
meteo_m$newTS<-as.Date(meteo_m$TS)

########## Meteoparameters and evaporation according to ERA5 at grid node nearest to Novo site: contribution by T. Naakka ##############################
path<-"/home/shevnina/Manuscripts/2021/AntarcticEVA/supplement/evaporation/"
setwd(path)
eva_era<-read.csv("Evaporation_Schirmacher_Oasis_from_ERA5.csv", header = F, sep=",", skip=3)             ##### results from ERA
eva_era$ts<-paste(eva_era$V3,eva_era$V2, eva_era$V1, sep="-")
eva_era$TS_UTS<-as.POSIXlt(eva_era$ts, format = "%Y-%m-%d")
eva_era$newTS<-as.Date(eva_era$TS_UTS)

########### Evaporation after the bulk aerodynamic method: contribution by T. Vihma ####################################################
path<-"/home/shevnina/Manuscripts/2021/AntarcticEVA/supplement/evaporation/"
setwd(path)
eva_bulk<-read.table("Bulk_method_results.txt", header = T, sep=" ")  ### input data is the measurements at site
eva_bulk$TS<-dat1$Timestamp_UTC
evap_BA<-aggregate(list(Evap=eva_bulk$E_mm_day),list(Timestamp= cut(eva_bulk$TS,"1440 min")), mean, na.rm=T)

eva_bulk_m<-read.table("Bulk_method_results_Maitri_input.txt", header = T, sep=" ") ### input data is the measurements at Maitri site
eva_bulk_m$newTS<-dat1$Timestamp_UTC
evap_BA_m<-aggregate(list(Evap_m=eva_bulk_m$E),list(Timestamp= cut(eva_bulk_m$newTS,"1440 min")), mean, na.rm=T)


############## Meteoparameters and evaporation after the EC measurements: contribution by M. Potes #########################################
####### data from the eddy covariance: the fluxes and evaporation after Potes et al., 2017  ############################################################
# Date Time, is the time stamp (TS) with the format “yy/mm/dd HH:MM:SS”
# npoints, is number of measurements
# u_starr (m/s), taur Kg(m s^2) are u_star and momentum flux
# Hcr (W/m^2) is a sensible heat flux humidity 
# LE_wplr (W/m^2), is a latent heat flux //  J m-2 s-1
# Evap (L/m^2) is an evaporation by Potes M. ## equal to mm of evaporation
# Fc_wplr (mg/(m^2 s)), Fc_wpl_umol (umol/(m^2 s)), is carbon flux with WPL correction, mg/(m^2 s) and umol/(m^2 s)
# CO2_conc (mg/m^3), CO2_PPM (PPM), CO2_sig_str (arb), are CO2 concentration (two unit) and the Irgason’s signal strength 
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

###################################### Lake Zub/Priyadarshini, experiment 2017-2018 ########################################################
############################################################################################################################################
path<-"/home/shevnina/Manuscripts/2020/Antarctica/evaporation/supplement/"
setwd(path)
Irgason<-"/home/shevnina/Manuscripts/2020/Antarctica/evaporation/supplement/20180101_20180207_EC_FLUX.txt"

variable_unit <- scan(Irgason, nlines = 2, sep=",", what = character())
dat1 <- read.csv(Irgason, skip = 2, header = FALSE, fill=TRUE)
dat1 <- dat1[,1:24]
colnames(dat1) <- c("TS","npoints","u_star","taur","Hcr","LE_wplr",
                    "Evap","Fc_wplr","Fc_wpl_umol","CO2_conc","CO2_PPM","CO2_sig_str",
                    "H2O_conc","H2O_sig_str","Temp_amb","Amb_Press","wind_speed","wind_dir_sonic",
                    "obukhov","zeta","rou","sigmaw_ustar","Xmax","XR90")

# to create proper UTC Timestamp, values corresponded to the end of time interval
dat1$Timestamp_UTC <- strptime(substr(dat1$TS,1,15),tz="UTC",format="%y/%m/%d %H:%M")



#correct for the wrong wind direction: added by Daniela Franz
dat1$wind_dir_sonic <- (dat1$wind_dir_sonic+43)%%360

######################### filtering the fluxes as suggested by Daniela Franz ###############################################
#Flag conc and fluxes: to flag signal strengths lower than 0.7 + wind direction from the land + individually outliers
### and to count the percentage of the filtered data

# to account to the percentage of data filtered due to the signal strengths (lower than 0.7)
all <- nrow(dat1[!is.na(dat1$H2O_conc),])
filt <- nrow(dat1[(!is.na(dat1$H2O_sig_str) & dat1$H2O_sig_str >= 0.7),])
(all-filt)*100/all # 0.11 %

# to wind direction from the land 
filt <- nrow(dat1[(dat1$wind_dir_sonic >= 105 & dat1$wind_dir_sonic <= 240),])
(all-filt)*100/all # 17.8 %

# to account the percentage of the data after the filtering 
# total after the Latent sensible flux and Evaporation
filt <- nrow(dat1[!is.na(dat1$Evap) & dat1$H2O_sig_str >= 0.7 & (dat1$wind_dir_sonic >= 105 & dat1$wind_dir_sonic <= 240),])
100-((all-filt)*100/all)
# H 
filt <- nrow(dat1[!is.na(dat1$Hcr) & (dat1$wind_dir_sonic >= 105 & dat1$wind_dir_sonic <= 240),])
100-((all-filt)*100/all)

# to calculate the Table XX with the filtering
dat1$H2O_conc_filter <- ifelse(dat1$H2O_sig_str >= 0.7 &
                                 (dat1$wind_dir_sonic >= 105 & dat1$wind_dir_sonic <= 240),
                               dat1$H2O_conc,
                               NA)
dat1$LE_filter <- ifelse(dat1$H2O_sig_str >= 0.7 &
                           (dat1$wind_dir_sonic >= 105 & dat1$wind_dir_sonic <= 240),
                         dat1$LE_wplr,
                         NA)

dat1$Evap_filter <- ifelse(dat1$H2O_sig_str >= 0.7 &
                             (dat1$wind_dir_sonic >= 105 & dat1$wind_dir_sonic <= 240),
                           dat1$Evap,
                           NA)

dat1$H_filter <- ifelse((dat1$wind_dir_sonic >= 105 & dat1$wind_dir_sonic <= 240),
                        dat1$Hcr,
                        NA)
dat1$tau_filter <- ifelse((dat1$wind_dir_sonic >= 105 & dat1$wind_dir_sonic <= 240),
                          dat1$taur,
                          NA)

dat1$u_star_filter <- ifelse((dat1$wind_dir_sonic >= 105 & dat1$wind_dir_sonic <= 240),
                             dat1$u_star,
                             NA)

############################## sub-results for the EC method ############################################################################################3
# the summary statistic for the filtered fluxes and evaporation
summary(dat1$H2O_conc,na.rm=T)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#1.024   2.085   2.479   2.530   2.896   6.068      13

summary(dat1$H2O_conc_filter,na.rm=T)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  1.101   2.076   2.471   2.512   2.871   6.068     332

summary(dat1$LE_filter,na.rm=T)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -17.96   52.82   83.55   90.64  125.94  216.05     337
summary(dat1$H_filter,na.rm=T)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -26.12   28.49   49.30   51.62   70.11  155.80     336

summary(dat1$Evap_filter,na.rm=T)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-0.0120  0.0373  0.0594  0.0644  0.0892  0.1546     337


##################### daily EC evaporation calculated with NA in the values: used further for the errors estimations(R, NS and SSigma)
evap_EC<-aggregate(list(Evap=dat1$Evap_filter),list(Timestamp= cut(dat1$Timestamp,"1440 min")), sum, na.rm=T)   
####################################################################################################################################
double1<-dat1$Evap_filter  # by mean
double2<-dat1$Evap_filter  # by median
## To fill the 30min evaporation missing data (337 value) with average 30min value = 0.0644
double1[is.na(double1)]<-mean(dat1$Evap_filter, na.rm=T)
#new summary after the fitting of evaporation missing
# summary(double1,na.rm=T)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.01205  0.04177  0.06438  0.06438  0.08113  0.15459
#sum = 115 mm per whole period
## To fill the 30min evaporation missing data (337 value) with median
double2[is.na(double2)]<-median(dat1$Evap_filter, na.rm=T)
dat1$Evap_filter<-double2
evap_ECbyMedian<-aggregate(list(Evap=dat1$Evap_filter),list(Timestamp= cut(dat1$Timestamp,"1440 min")), sum, na.rm=T)
# sum(evap_ECbyMedian$Evap) [1] 114.133, the daily evaporation measured

# irg_day consists of the daily averages for the wind_speed, wind_direction, air_temperature, H2O concentration, fluxes and the evaporation 
##Timestamp        LE[Wm-2]          H[Wm-2]        WS[ms-1]         AT[C]      H2O[g sm-3]        P[PA]
###eva[mm day-1] SpecHumid[Pa]   SatPre[Pa] SatSpecHumid[Pa]       RH[%]   VapPre[Pa]

irg_day <- strptime(dat1$Timestamp_UTC, format="%Y-%m-%d")
irg_day <-aggregate(list(LE = dat1$LE_filter, H = dat1$H_filter, WS = dat1$wind_speed, AT = dat1$Temp_amb, H2O = dat1$H2O_conc_filter, P= dat1$Amb_Press*10),list(Timestamp= cut(dat1$Timestamp,"1440 min")), mean, na.rm=T)         ### Amp_Pres is given in Pa due to dividing by 10  ###

# to calculate the daily evaporation from the open water surface of the lakes
irg_day$eva<-"NA"
eva_day <-aggregate(list(Evap=dat1$Evap_filter),list(Timestamp= cut(dat1$Timestamp,"1440 min")), sum, na.rm=T)
irg_day$eva<-eva_day$Evap

# to calculate the relative humidity from the H2O concentration
### 1. the specific humidity ?g kg-1? #############
R<-287.05                                                                  #############the universal gas constant for a dry air, J kg-1 K-1
irg_day$SpecHumid <- R*irg_day$H2O*(irg_day$AT+273.15)/irg_day$P           ######### specific humidity

# to calculate the saturated vapor pressure (Pa) after Teten's formula
b1<- 611.3                                                                  # the coefficients according to Stull, 2017
b2<- 17.2694                                                                # #---#
b3<- 35.86                                                                  # #---#
irg_day$SatPre<-"NA"
irg_day$SatPre<-b1*exp((b2*irg_day$AT)/(irg_day$AT + 273.16-b3))

# to calculate the saturation specific humidity, g kg-1 ############### 
irg_day$SatSpecHumid<-0.622 *irg_day$SatPre/(irg_day$P-(1-0.622)*irg_day$SatPre)*1000   #### to consider in units with the specific humidity
irg_day$RH <-"NA"
irg_day$RH<-100*irg_day$SpecHumid/irg_day$SatSpecHumid
irg_day$VapPre<-irg_day$RH*irg_day$SatPre/100
irg_day$newTS = as.Date(irg_day$Timestamp, "%Y-%m-%d")


###############################################################
#### new RH for Irgason after Hoeltgebaum et al. (2020) 
SpecHumid <- R*irg_day$H2O*0.001*(irg_day$AT+273.15)/622
SatPre<-0.611*exp((12.27*irg_day$AT)/(irg_day$AT + 273.16))
RH<-100*SpecHumid/SatPre


##############################################################################################################################################
#### Evaporation after empirical formulas: contribution by E. Shevnina



############################################################################################################
############################################## RESULTS #####################################################

######################## Figure 3 ####################################################################
### (a)
plot(novo1718daily$newTS, novo1718daily$AT, pch=10, axes=FALSE, ylim=c(-15,5), xlab="Date", ylab="", type="l",col="blue")    ##### Novo air temperature
lines(meteo_m$newTS,meteo_m$AT, col="red")
axis(2, ylim=c(0,20),col="black",las=1)
mtext("C",side=2,line=2.5)
axis.Date(1, at=seq(as.Date("2017/12/01"), as.Date("2018/02/28"),by="11 day"), format="%Y/%m/%d", cex.lab=1.5)
box()
grid()

### (b)
plot(novo1718daily$newTS, novo1718daily$RH, pch=10,  xlab="", ylab="", ylim=c(0,100), axes=FALSE, type="l", col="blue")
lines(meteo_m$newTS,meteo_m$RH, col="red")
lines(irg_day$newTS, irg_day$RH, col="green") 
mtext("%",side=2,line=2.5) 
axis(2,  ylim=c(20,80), col="black",col.axis="black",las=1)
axis.Date(1, at=seq(as.Date("2017/12/01"), as.Date("2018/02/28"),by="11 day"), format="%Y/%m/%d", cex.lab=1.5)
grid()
box()


### (c)
plot(novo1718daily$newTS, novo1718daily$WS, pch=10,  xlab="", ylab="", ylim=c(0,max(novo1718daily$WS)), axes=FALSE, type="l", col="blue")
lines(meteo_m$newTS,meteo_m$WS, col="red")
mtext("ms-1",side=2,line=2.5) 
axis(2,  ylim=c(20,80), col="black",col.axis="black",las=1)
axis.Date(1, at=seq(as.Date("2017/12/01"), as.Date("2018/02/28"),by="11 day"), format="%Y/%m/%d", cex.lab=1.5)
grid()
box()

### (d)

plot(novo1718daily$newTS, novo1718daily$NSOL, pch=10,  xlab="", ylab="", ylim=c(50,max(novo1718daily$NSOL)), axes=FALSE, type="l", col="blue")
mtext("W m-2",side=2,line=2.5) 
axis(2,  ylim=c(20,80), col="black",col.axis="black",las=1)
axis.Date(1, at=seq(as.Date("2017/12/01"), as.Date("2018/02/28"),by="11 day"), format="%Y/%m/%d", cex.lab=1.5)
grid()
box()
###########################################################################################################################



##### Figure 5 
### (a)
plot(irg_day$newTS, irg_day$AT, xlab="", ylab="", ylim=c(min(irg_day$A),max(irg_day$AT)), axes=FALSE, type="l", col="black")
lines(meteo_m$newTS, meteo_m$AT, xlab="", ylab="", col="red")
#mtext("C",side=2,line=2.5) 
axis(2,  ylim=c(20,80), col="black",col.axis="black",las=1)
axis.Date(1, at=seq(as.Date("2018/01/01"), as.Date("2018/02/08"), by="5 day"), format="%y/%m/%d", cex.lab=1.5)
grid()
box()

### (b)
plot(irg_day$newTS, irg_day$WS, xlab="", ylab="", ylim=c(0,max(irg_day$WS)), axes=FALSE, type="l", col="black")
lines(meteo_m$newTS, meteo_m$WS, xlab="", ylab="", col="red")
mtext("ms-1",side=2,line=2.5) 
axis(2,  ylim=c(20,80), col="black",col.axis="black",las=1)
axis.Date(1, at=seq(as.Date("2018/01/01"), as.Date("2018/02/08"),by="5 day"), format="%Y/%m/%d", cex.lab=1.5)
grid()
box()

#### (c)
plot(irg_day$newTS, irg_day$RH, xlab="", ylab="", ylim=c(0,100), axes=FALSE, type="l", col="black")
lines(meteo_m$newTS, meteo_m$RH, xlab="", ylab="", col="red")
mtext("%",side=2,line=2.5) 
axis(2,  ylim=c(20,80), col="black",col.axis="black",las=1)
axis.Date(1, at=seq(as.Date("2018/01/01"), as.Date("2018/02/08"),by="5 day"), format="%Y/%m/%d", cex.lab=1.5)
grid()
box()

####### Figure 8
### (b)
plot(evap_EC$Evap, sub_maitri$OD, xlab="mm day-1", ylab="mm day-1", pch=12, col="black", xlim=c(0,6), ylim=c(0,6),cex.axis=1.5, cex.main=1.5)  ## ERA
points(evap_EC$Evap, eva_emp$IR_evaOD,  xlab="", ylab="", pch=12, col="red")
grid()

### (c)
plot(evap_EC$Evap, sub_maitri$PN, xlab="mm day-1", ylab="mm day-1", pch=12, col="black", xlim=c(0,6), ylim=c(0,6),cex.axis=1.5, cex.main=1.5)  ## ERA
points(evap_EC$Evap, eva_emp$IR_evaPN,  xlab="", ylab="", pch=12, col="red")
grid()

### (d)
plot(evap_EC$Evap, sub_maitri$DP, xlab="mm day-1", ylab="mm day-1", pch=12, col="black", xlim=c(0,6), ylim=c(0,6),cex.axis=1.5, cex.main=1.5)  ## ERA
points(evap_EC$Evap, eva_emp$IR_evaDP,  xlab="", ylab="", pch=12, col="red")
grid()

### (a)
plot(evap_EC$Evap, evap_BA_i$Evap, xlab="mm day-1", ylab="mm day-1", pch=12, col="black", xlim=c(0,6), ylim=c(0,6),cex.axis=1.5, cex.main=1.5)  ## ERA
points(evap_EC$Evap, evap_BA_m$Evap_m,  xlab="", ylab="", pch=12, col="red")
grid()
########################################################################################################################








##################### Performance of the indirect methods #####################
path<-"/home/shevnina/Manuscripts/2021/AntarcticEVA/supplement/evaporation/"
setwd(path)
eva_emp<-read.csv("EvaporationZub_Irgason.csv", header = T, sep=",")                                      ##### results after the semi-empirical equations, EC and meteo ay Irgason site
eva_mai<-read.csv("Evaporation_Priyadarshini_Lake_Antarctica_2018_Pankaj.csv", header = T, sep=",")       ##### results after Odrova and meteo Maitri site

eva_era<-read.csv("Evaporation_Schirmacher_Oasis_from_ERA5.csv", header = F, sep=",", skip=3)             ##### results from ERA
eva_era$ts<-paste(eva_era$V3,eva_era$V2, eva_era$V1, sep="-")
eva_era$TS_UTS<-as.POSIXlt(eva_era$ts, format = "%Y-%m-%d")
era_subset<-subset(eva_era$V5, eva_era$TS_UTS>="2018-01-01" & eva_era$TS_UTS <="2018-02-07")              #### subset for the period of observation


############# s/sigma criterio, (Popov,  1979, p. 58), S_Sigma criterio #######################
m1<-mean(eva_emp$IR_eva) 
l1<-length(eva_emp$IR_eva)                        #
s2<-sum((eva_emp$IR_eva-evap_BA$Evap)**2)/(l1-2)            ############# Bulk method ### 1.07
s2<-sum((eva_emp$IR_eva-eva_mai$E..mm.day.)**2)/(l1-2)      ############# Odrova + Maitri ### 1.28
s2<-sum((eva_emp$IR_eva-eva_emp$IR_evaOD)**2)/(l1-2)        ##### Odrova +IR ### 1.64
s2<-sum((eva_emp$IR_eva-eva_emp$IR_evaPN)**2)/(l1-2)        #### Penmann + Irgason ### 1.15
s2<-sum((eva_emp$IR_eva-eva_emp$IR_evaDP)**2)/(l1-2)        ##### Doorenboss + Irgason ### 1.13

############ Nash Sutcliffe index #########################
m1<-mean(eva_emp$IR_eva) 
v1<-sum((evap_BA$Evap-eva_emp$IR_eva)**2)          ########## Bulk method ##### -0.10
1-v1/v2
v1<-sum((eva_mai$E..mm.day.-eva_emp$IR_eva)**2)    #### Odrova + Maitri ##### -0.56
v1<-sum((eva_emp$IR_evaDP-eva_emp$IR_eva)**2)      ##### Doorenboss #### -0.23
v1<-sum((eva_emp$IR_evaPN-eva_emp$IR_eva)**2)      ###### Penmann #### -0.26
v1<-sum((eva_emp$IR_evaOD-eva_emp$IR_eva)**2)      ######## Odrova +Irgason ### -1.63



