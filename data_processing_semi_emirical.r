# Supplement to Shevnina et al., 2021 (Sections 3.1.2 The semi-empirical equations
# Authorship: eshevnina@gmail.com
# phone: +358449185446

library(lubridate)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting

################################################## Observations aggregated to daily values ##################################################
#############################################################################################################################################
path<-"/home/shevnina/Manuscripts/2020/Antarctica/evaporation/supplement"
setwd(path)
datOne<-read.table("20171201_20180210_Novo.txt", sep=";", header=T) # meteo DJF 2017/2018
datTwo<-read.table("20191201_20200131_Novo.txt", sep=";", header=T) # meteo DJF 2019/2020

################################################ Novo station SYNOP meteo ##################### #############################################################################################
################### daily aggregation #######################################################
datOne$one_min<-strptime(datOne$Dates_Times, format="%Y-%m-%d %H:%M")
dat24hour<- aggregate(list(AT = datOne$Air_temperature, ST = datOne$Soil_Temperature, RH = datOne$Relative_Humidity, WS = datOne$Wind_speed, NSOL = datOne$Incoming_solar_radiation),list(Timestamp= cut(datOne$one_min, "24 hour" )), mean)
dat24hour$newTS = as.Date(dat24hour$Timestamp, "%Y-%m-%d")
novo1718daily<- dat24hour         # 2017/2018
#  Timestamp         AT       ST       RH       WS     NSOL      newTS
#1 2017-12-01  2.6711111       NA 41.85069 5.307500 376.5026 2017-12-01
#2 2017-12-02  0.7939583       NA 37.02500 3.787847 384.4626 2017-12-02

# season 2019/2020 aggregation to daily values
datTwo$one_min<-strptime(datTwo$Dates_Times, format="%Y-%m-%d %H:%M")
dat24hour<- aggregate(list(AT = datTwo$Air_temperature, ST = datTwo$Soil_Temperature, RH = datTwo$Relative_Humidity, WS = datTwo$Wind_speed, NSOL = datTwo$Incoming_solar_radiation),list(Timestamp= cut(datTwo$one_min,"24 hour")), mean)
dat24hour$newTS = as.Date(dat24hour$Timestamp, "%Y-%m-%d")
novo1920daily <- dat24hour                  ## daily averages meteo from NOVO

#head(novo1920daily)
#   Timestamp         AT        ST       RH        WS     NSOL      newTS
#1 2019-12-01 -3.4656944  4.787917 57.32708  3.601528 283.4272 2019-12-01
#2 2019-12-02 -4.1028472  5.186389 64.80694  2.705208 277.3124 2019-12-02

############################################## Results #########################
# 4.1: Weather conditions during the summer seasons 2017 – 2018 and 2019 – 2020 
######################## Figure XX a, b ####################################################################
# daily air temp and relative humidity
pdf("/home/shevnina/Manuscripts/2020/Antarctica/evaporation/supplement/AirTemp_RH1920.pdf")
par(mfrow=c(2,1)) 
par(mar=c(5, 4, 4, 6) + 0.1)
plot(novo1920daily$newTS, novo1920daily$AT, pch=10, axes=FALSE, ylim=c(-8,5), xlab="Date", ylab="", type="b",col="blue")
axis(2, ylim=c(0,20),col="blue",las=1)
mtext("Air Temp, C",side=2,line=2.5)
axis.Date(1, at=seq(as.Date("2019/12/01"), as.Date("2020/02/28"),by="11 day"), format="%y/%m/%d", cex.lab=1.5)
box()
grid()
par(mar=c(5, 4, 4, 6) + 0.1)
#par(new=TRUE)
plot(novo1920daily$newTS, novo1920daily$RH, pch=10,  xlab="", ylab="", ylim=c(0,100), axes=FALSE, type="b", col="red")
mtext("RH, %",side=2,line=2.5) 
axis(2,  ylim=c(20,80), col="red",col.axis="red",las=1)
axis.Date(1, at=seq(as.Date("2019/12/01"), as.Date("2020/02/28"),by="11 day"), format="%y/%m/%d", cex.lab=1.5)
grid()
box()
dev.off()

pdf("/home/shevnina/Manuscripts/2020/Antarctica/evaporation/supplement/AirTemp_RH1718.pdf")
par(mfrow=c(2,1)) 
par(mar=c(5, 4, 4, 6) + 0.1)
plot(novo1718daily$newTS, novo1718daily$AT, pch=10, axes=FALSE, ylim=c(-8,5), xlab="Date", ylab="", type="b",col="blue")
axis(2, ylim=c(0,20),col="blue",las=1)
mtext("Air Temp, C",side=2,line=2.5)
axis.Date(1, at=seq(as.Date("2017/12/01"), as.Date("2018/02/28"),by="11 day"), format="%y/%m/%d", cex.lab=1.5)
box()
grid()
par(mar=c(5, 4, 4, 6) + 0.1)
#par(new=TRUE)
plot(novo1718daily$newTS, novo1718daily$RH, pch=10,  xlab="", ylab="", ylim=c(0,100), axes=FALSE, type="b", col="red")
mtext("RH, %",side=2,line=2.5) 
axis(2,  ylim=c(20,80), col="red",col.axis="red",las=1)
axis.Date(1, at=seq(as.Date("2019/12/01"), as.Date("2020/02/28"),by="11 day"), format="%y/%m/%d", cex.lab=1.5)
grid()
box()
dev.off()

# net solar radiation
pdf("/home/shevnina/Manuscripts/2020/Antarctica/evaporation/supplement/Net_radiation.pdf")
par(mfrow=c(2,1)) 
par(mar=c(5, 4, 4, 6) + 0.1)
plot(novo1718daily$newTS, novo1718daily$NSOL, pch=10, axes=FALSE, ylim=c(0,500), xlab="Date", ylab="", type="b",col="orange")
axis(2, ylim=c(0,20),col="black",las=1)
mtext("NetSol, W m-2",side=2,line=2.5)
axis.Date(1, at=seq(as.Date("2017/12/01"), as.Date("2018/02/28"),by="11 day"), format="%y/%m/%d", cex.lab=1.5)
box()
grid()
par(mar=c(5, 4, 4, 6) + 0.1)
#par(new=TRUE)
plot(novo1920daily$newTS, novo1920daily$NSOL, pch=10,  xlab="", ylab="", ylim=c(0,500), axes=FALSE, type="b", col="orange")
mtext("NetSOL, W m-2",side=2,line=2.5) 
axis(2,  ylim=c(20,80), col="black",col.axis="black",las=1)
axis.Date(1, at=seq(as.Date("2019/12/01"), as.Date("2020/02/28"),by="11 day"), format="%y/%m/%d", cex.lab=1.5)
grid()
box()
dev.off()

### results for the section of 4.1: weather conditions in seasons 2017/2018 and 2019/2020
summary(novo1718daily)	
summary(novo1920daily)

## daily wind speed and direction : WIND ROSE 


##############################################################################################################################################
##################################################### Data and Methods ###############################

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
###################################### Lake Glubokoe, experiment 2019/2020 not yet done !!! #################################################

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



######################### filtering the fluxes by Daniela Franz ###############################################
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

# to calculate the Table XX with the filetring
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

############################## results ############################################################################################3
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


##################### daily EC evaporation calculated with NA in the values!!!! used for the errors estimations(NS and SSigma)
evap_EC<-aggregate(list(Evap=dat1$Evap_filter),list(Timestamp= cut(dat1$Timestamp,"1440 min")), sum, na.rm=T)   #!!!! TO VERYFY!!!!!!!
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
R<-287.05      #############the universal gas constant for a dry air, J kg-1 K-1
irg_day$SpecHumid <- R*irg_day$H2O*(irg_day$AT+273.15)/irg_day$P
# to calculate the saturated vapoir pressure (Pa) after Teten's formula
b1<- 611.3     # the coefficients according to Stull, 2017
b2<- 17.2694   # #---#
b3<- 35.86     # #---#
irg_day$SatPre<-"NA"
irg_day$SatPre<-b1*exp((b2*irg_day$AT)/(irg_day$AT + 273.15-b3))
# to calculate the saturation specific humidity, g kg-1 ############### 
irg_day$SatSpecHumid<-0.622 *irg_day$SatPre/(irg_day$P-(1-0.622)*irg_day$SatPre)*1000  #### to conside in units with the specific humidity
irg_day$RH <-"NA"
irg_day$RH<-100*irg_day$SpecHumid/irg_day$SatSpecHumid
irg_day$VapPre<-irg_day$RH*irg_day$SatPre/100


################################# figures ####################################
############################################ figures ###############################################################################################
# Fig. XX. 30 minute series of the sensible heat flux humidity (H, W m-2), the latent heat flux (LE, W m-2), and H2O concentration (g m-3) measured by the flux tower (red dots are the values filtered out after the processing). 
# Fig. 4 in Shevnina et al. 2020.
pdf("/home/shevnina/Manuscripts/2020/Antarctica/evaporation/supplement/FigureXX.pdf")
par(mfrow=c(3,1)) 
# the sensible heat flux humidity (H, W m-2)
plot(as.POSIXct(dat1$Timestamp_UTC,tz="UTC"),dat1$Hcr, col="red", pch =16, xlab="2018", ylab="", xaxt="n") # the points to be excluded
points(as.POSIXct(dat1$Timestamp_UTC,tz="UTC"),dat1$H_filter,col="black", pch=16, xlab="", ylab="", xaxt="n") # the points to be used in the further calculatons
axis.POSIXct(1, at = seq(min(dat1$Timestamp_UTC), max(dat1$Timestamp_UTC), "day"), labels = TRUE, format="%m/%d")
axis(2, ylim=c(0,250),col="black")
mtext("Hc, W m-2",side=2,line=2.5)
#box()
grid()

# the latent heat flux (LE, W m-2)
plot(as.POSIXct(dat1$Timestamp_UTC,tz="UTC"),dat1$LE_wplr, col="red", pch =16, xlab="2018", ylab="", xaxt="n")
points(as.POSIXct(dat1$Timestamp_UTC,tz="UTC"),dat1$LE_filter,col="orange", pch=16, xlab="", ylab="", xaxt="n") # the points to be used in the further calculatons
axis.POSIXct(1, at = seq(min(dat1$Timestamp_UTC), max(dat1$Timestamp_UTC), "day"), labels = TRUE, format="%m/%d")
#axis(2, ylim=c(0,250),col="black",las=1)
mtext("LE, W m-2",side=2,line=2.5)
grid()

# H2Oconcentration
plot(as.POSIXct(dat1$Timestamp,tz="UTC"),dat1$H2O_conc, col="red", pch =16, xlab="Date", ylab="",xaxt="n") # the points to be excluded
points(as.POSIXct(dat1$Timestamp,tz="UTC"),dat1$H2O_conc_filter,col="blue", pch=16, xlab="", ylab="", xaxt="n") # the points to be used in the further calculatons
axis.POSIXct(1, at = seq(min(dat1$Timestamp_UTC), max(dat1$Timestamp_UTC), "day"), labels = TRUE, format="%m/%d")
mtext("a, g m-3",side=2,line=2.5)
grid()
dev.off()



# Fig. XX: The daily values of the evaporation, Lake Zub/Priyadarshini
pdf("/home/shevnina/Manuscripts/2020/Antarctica/evaporation/supplement/FigureXX.pdf")
par(mar=c(5, 4, 4, 6) + 0.1)
plot(as.POSIXct(irg_day$Timestamp,tz="UTC"),irg_day$eva, col="black", pch =16, xlab="2018", ylab="",xaxt="n", yaxt="n", type="b")
axis(2, ylim=c(0,8), col="black",col.axis="black",las=1)
axis.POSIXct(1, at = seq(min(as.Date(irg_day$Timestamp)), max(as.Date(irg_day$Timestamp)), "day"), labels = TRUE, format="%m/%d")
mtext("Evaporation, mm day-1",side=2,line=2.5)
par(new=TRUE)
plot(as.POSIXct(irg_day$Timestamp,tz="UTC"), irg_day$WS, col="blue", pch =16, xlab="", ylab="",xaxt="n", yaxt="n",type="b")
axis(4, ylim=c(0,20), col="blue",col.axis="blue",las=1)
mtext("WS, m s-1",side=4,line=2.5)
grid()
dev.off()

##############################################################################################################################################
####################################### Evaporation after semi-empirical formulas #############################################################
##################### CONSTANTS #########################################
b1<- 611.3     # the coefficients according to Stull, 2017
b2<- 17.2694   # #---#
b3<- 35.86     # #---#

########################## the forcing calculated from the Novo SYNOP (the case SSR, with data "novo1718daily") #############################
## Penmann is according to Penmann, 1948 after Tannay et al., 2007 
## Dooren is according to Doorenbos and Pruitt, 1975 
## to calculate the water vapor saturation pressure, [Pa] ###########################
## SatPres is the saturation pressure, Pa , Stull, 2017 ##########################
## VapPress is the actual vapor pressure, Pa ###################################################################
novo1920daily$SatPre<-"NA"
novo1920daily$SatPre<-b1*exp((b2*novo1920daily$AT)/(novo1920daily$AT + 273.15-b3))
novo1920daily$VapPre<-novo1920daily$RH*novo1920daily$SatPre/100

gl_eva<-list("date"=novo1920daily$Timestamp)   #### evaporation over the lake Glubokoe, season 2019/2020
gl_eva$Penmann<- "NA"
gl_eva$Penmann<- 0.26*(1+0.54*novo1920daily$WS)*(novo1920daily$SatPre/100-novo1920daily$VapPre/100) ######## vapour pressure is given in millibars
gl_eva$Dooren<- 0.26*(1+0.86*novo1920daily$WS)*(novo1920daily$SatPre/100-novo1920daily$VapPre/100)

pr_eva<-list("date"=novo1718daily$Timestamp)   #### evaporation over the lake Priyadarshini, season 2017/2018
novo1718daily$SatPre<-"NA"
novo1718daily$SatPre<-b1*exp((b2*novo1718daily$AT)/(novo1718daily$AT + 273.15-b3))
novo1718daily$VapPre<-novo1718daily$RH*novo1718daily$SatPre/100

pr_eva$PN<- 0.26*(1+0.54*novo1718daily$WS)*(novo1718daily$SatPre/100-novo1718daily$VapPre/100) ######## vapour pressure is given in millibars
pr_eva$DN<- 0.26*(1+0.86*novo1718daily$WS)*(novo1718daily$SatPre/100-novo1718daily$VapPre/100)
pr_eva$day<-as.POSIXct(pr_eva$date)

# selection of the data covering the only observational period
date1 <- as.POSIXct(irg_day$Timestamp[1], format="%Y-%m-%d")
date2 <- as.POSIXct(irg_day$Timestamp[length(irg_day$Timestamp)])

subset1<-list(subset(pr_eva$day, pr_eva$day>=date1 & pr_eva$day<=date2))  ##### the Novo SYNOP for the period of 01.01.2018-07.02.2018
subset1$PN<-"NA"
subset1$DN<-"NA"
subset1$PN<-subset(pr_eva$PN, pr_eva$day>=date1 & pr_eva$day<=date2)
subset1$DN<-subset(pr_eva$DN, pr_eva$day>=date1 & pr_eva$day<=date2)

#int <- interval(date1, date2)
#selXX<-pr_eva[pr_eva$day %within% int,]

####### from the Irgason data ("irg_day") season 2017/2018
#irg_day$date<-as.POSIXct(irg_day$Timestamp, format="%Y-%m-%d")
prlake_EC<-list("date" =irg_day$Timestamp)
prlake_EC$PI<-"NA"
prlake_EC$DI<-"NA"
prlake_EC$PI<-0.26*(1+0.54*irg_day$WS)*(irg_day$SatPre/100-irg_day$VapPre/100) ### vapour pressure is given in millibars
prlake_EC$DI<- 0.26*(1+0.86*irg_day$WS)*(irg_day$SatPre/100-irg_day$VapPre/100)

########################### results #########################################################
## comparison the daily evaporation EC (measured) and semi-empirical formulas (modeled): period of 01.01.2018 to 07.02.2018
## prlake_EC inludes the modeled daily evaporation mm day-1 (prlake_EC$eva and modeled prlake_EC$PI and prlake_EC$DI
## subset1 includes the modeled daily evaporation mm day-1 after semi-empirical methods (pr_eva$PI and pr_eva$DI)
# cor is the correlation to the EC evaporationmean

#summary(prlake_EC$PI)                               3.9 +/- 1.9  Penman , 1948 + Irgason observations (cor=0.62, k=mod/obs =1.3) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.013   2.368   3.935   3.908   5.375   8.558 
#> sd(prlake_EC$PI) [1] 1.887193
# sum(prlake_EC$PI) [1] 148.4944

#> sd(prlake_EC$DI) [1] 2.91594                      5.7 +/- 2.9 Dooren and T, 1975 + Irgason observations (cor=0.62, k=mod/obs = 1.9)
# summary(prlake_EC$DI)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.300   3.312   5.707   5.699   7.867  13.056 
# sum(prlake_EC$DI) [1] 216.5591

#summary(subset1$DN)                                    4.7 +/- 2.5  Dooren and .. 1975, + Novo SYNOP observations (cor = 0.67, k=mod/obs =0.9)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4367  2.3974  4.4851  4.7047  6.7728  9.1370 
# sum(subset1$DN) [1] 178.7773

#summary(subset1$PN)                              3.2 +/- 1.7    Penman, 1948 + Novo SYNOP observations (cor = 0.0.68, k=mod/obs =1.2) BEST???
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3328  1.6677  3.0823  3.1943  4.5193  5.9358 
# sd(subset1$PN) +/- 1.652359
# sum(subset1$PN) [1] 121.3838

# summary(irg_day$eva)					3.0 +/- 1.0 Eddy Covariance reference 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.502   2.180   2.731   3.003   4.007   4.999
# sd(irg_day$eva) [1] = 1.054796
# sum(irg_day$eva) [1] = 114.133 mm per period


################ for the DJF months of 2017/2018
#summary(pr_eva$DN)                                      5.2 +/- 2.5  
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4367  3.0832  5.4830  5.2303  6.6556 12.3326
# sum(pr_eva$DN)  = 470.7291 mm per season (DJF)
 
#summary(pr_eva$PN)                                     3.5 +/- 1.7
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3328  2.1633  3.7174  3.5424  4.4588  8.1062
# sum(pr_eva$PN)   = 318.8186

####### the Nash-Sutcliffe efficiency coefficient, Nash and Sutcliffe, 1970 in Tanny et al., 2008 q(
q1<-(irg_day$eva-subset1$PN)**2
q2<-(irg_day$eva-mean(irg_day$eva))**2
-1-sum(q1)/sum(q2)                        ########### -2.30

q1<-(irg_day$eva-subset1$DN)**2
1-sum(q1)/sum(q2)                         ########### -5.14

q1<-(irg_day$eva-prlake_EC$PI)**2         ###########  -1.65
1-sum(q1)/sum(q2)

q1<-(irg_day$eva-prlake_EC$DI)**2         ########### -10.8
1-sum(q1)/sum(q2)


############# s/sigma criterio, (Popov,  1979, p. 58), S_Sigma criterio #######################
m1<-mean(irg_day$eva) 
l1<-length(irg_day$eva)      #

s2<-sum((irg_day$eva-subset1$PN)**2)/l1     ########### 1.14
sigma2<-sum((irg_day$eva-m1)**2)/l1
sqrt(s2)/sqrt(sigma2)

s2<-sum((irg_day$eva-subset1$DN)**2)/l1     ########### 2.48
sqrt(s2)/sqrt(sigma2)

s2<-sum((irg_day$eva-prlake_EC$PI)**2)/l1   ########### 1.62
sqrt(s2)/sqrt(sigma2)

s2<-sum((irg_day$eva-prlake_EC$DI)**2)/l1   ########### 3.44
sqrt(s2)/sqrt(sigma2)

########################################################################################################################
# sum(pr_eva$PN) [1] 318.8186 Total evaporation over Dec-Feb
## December 2017
date1<-as.POSIXct("2017-12-01", format="%Y-%m-%d")
date2<-as.POSIXct("2017-12-31", format="%Y-%m-%d")
dec_eva<-list(subset(pr_eva$day, pr_eva$day>=date1 & pr_eva$day<=date2))
dec_eva$PN<-"NA"
dec_eva$PN<-subset(pr_eva$PN, pr_eva$day>=date1 & pr_eva$day<=date2)
# sum(dec_eva$PN) [1] 123.9243
date1<-as.POSIXct("2018-01-01", format="%Y-%m-%d")
date2<-as.POSIXct("2018-01-31", format="%Y-%m-%d")
jan_eva<-list(subset(pr_eva$day, pr_eva$day>=date1 & pr_eva$day<=date2))
jan_eva$PN<-"NA"
jan_eva$PN<-subset(pr_eva$PN, pr_eva$day>=date1 & pr_eva$day<=date2)
# sum(jan_eva$PN) [1] 89.8

date1<-as.POSIXct("2018-02-01", format="%Y-%m-%d")
date2<-as.POSIXct("2018-02-28", format="%Y-%m-%d")
feb_eva<-list(subset(pr_eva$day, pr_eva$day>=date1 & pr_eva$day<=date2))
feb_eva$PN<-"NA"
feb_eva$PN<-subset(pr_eva$PN, pr_eva$day>=date1 & pr_eva$day<=date2)
#sum(feb_eva$PN) [1] 105.0635
#jan+feb = 195 mm or about 28 % more than in Dhote et al., 2020 (167 mm for Jan-Feb)




