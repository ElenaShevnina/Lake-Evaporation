# Supplement to Shevnina, Potes, Vihma and Naakka 2025
# Authorship: eshevnina@gmail.com
# phone: +358449185446
# October 2025

library(stringr)
library(dplyr)
library(tidyr) 

######################### DATA #####################################################################################
# The data for the bulk-aerodynamic method: 30 min data in the format
#"Timestamp_UTC": time stamp
#"Evap": evaporation mm 30min-1 
#"Temp_amb": the air temperature, C, irgason
#"wind_speed": m s-1 Irgason
#"Amb_Press": ambient pressure, kPa,Irgason
#"RH": relative humidity, % Irgason  after Hoeltgebaum et al. (2020)
#"TW": LSWT, C, Hobo

######################################## METHODS ##################################################################################################
to_calulate_SVD <- function (T, LWST) {
  SVD <- 0.6113 * exp(17.27*LWST/(LWST + 237.3)) - 0.6113 * exp(17.27*T/(T + 237.3))                                      # saturation vapor pressure deficit, kPa. Tetens formula in Stull 2007
  return(SVD)
}


### converting 39 min data to daily data
to_daily_evap_BA <- function(dat)
{
  dat$hh<-strptime(dat$Timestamp_UTC, format="%Y-%m-%d %H:%M")
  d24h <- aggregate(list(EEC = dat$Evap, Ehk = dat$EBA_hk, Ewd = dat$EBA_wd, Ean = dat$EBA_an, Eaf = dat$EBA_af), list(Timestamp= cut(dat$hh, "24 hour" )), sum, na.rm=TRUE)
  d24h$newTS <- as.Date(d24h$Timestamp, "%Y-%m-%d")
  return(d24h)
}

############ Scores for in direct methods


### Table 4 Scores by Copilot
to_scores_BA_CP <- function(eva_day) {
  eva_day <- setNames(eva_day, sub("^Ehk$", "eva_A", names(eva_day)))
  eva_day <- setNames(eva_day, sub("^Ewd$", "eva_B", names(eva_day)))
  eva_day <- setNames(eva_day, sub("^Ean$", "eva_C", names(eva_day)))
  eva_day <- setNames(eva_day, sub("^Eaf$", "eva_D", names(eva_day)))
  eva_day <- setNames(eva_day, sub("^Evap$", "eva_ec", names(eva_day)))

  vars <- c("eva_A", "eva_B", "eva_C", "eva_D")
  rmse <- sapply(vars, function(v) sqrt(mean((eva_day$eva_ec - eva_day[[v]])^2, na.rm=TRUE)))
  m1 <- mean(eva_day$eva_ec, na.rm=TRUE)
  l1 <- length(eva_day$eva_ec)
  sigma <- sqrt(sum((eva_day$eva_ec - m1)^2, na.rm=TRUE) / l1)
  ssg <- sapply(vars, function(v) sqrt(sum((eva_day[[v]] - eva_day$eva_ec)^2, na.rm=TRUE) / (l1 - 2)) / sigma)          #ssigma criterio, Popov 1979
  data.frame(names = c("HK", "WD", "AN", "AF"), RMSE = rmse, SSg = ssg, stringsAsFactors = FALSE)
}

## summary statistics, Table 3 by Copilot
get_stats_BA <- function(eva_day) {
  methods <- c("Ehk", "Ewd", "Ean", "Eaf")
  stats <- sapply(methods, function(m) {
    x <- eva_day[[m]]
    s <- sum(x, na.rm = TRUE)
    m <- mean(x, na.rm = TRUE)
    mx <- max(x, na.rm = TRUE)
    err <- sqrt(sum((x - m)^2, na.rm = TRUE) / (length(x) - 1)) / sqrt(length(x))
    c(SUM = s, MEAN = m, MAX = mx, ERR_MEAN = err)
  })
  stats_df <- as.data.frame(t(stats))
  stats_df$Method <- methods
  stats_df <- stats_df[, c("Method", "SUM", "MEAN", "MAX", "ERR_MEAN")]
  rownames(stats_df) <- NULL
  return(stats_df)
}

######################################## Bulk-aerodynamic method ##########################################################
# E is the evaporation (in kg m−2 s−1) (which we in the following convert to mm d−1 ),
# ρ is the air density (in kg m3 ) == 1.2 kg m3
# CEz is the turbulent transfer coefficient for moisture (unit less)
# qs is the saturation specific humidity at the water surface of the lake (kg kg−1 ) 
# qaz is the air saturation specific humidity (kg kg−1 )
# WS is the wind speed (m s−1 ).

#### Lake Glubokoe
inBA_LG <- read.csv("/home/shevnina/Manuscripts/2025/Ant/september/result_tab/BA_LG.csv")
inBA_LG$Timestamp_UTC <- as.POSIXct(inBA_LG$Timestamp_UTC, format="%Y-%m-%d %H:%M:%S", tz="UTC")

# to read the reference daily evaporation (EC) data
dEVap_LG <- read.csv("/home/shevnina/Manuscripts/2025/Ant/september/result_tab/evap_daily_LG.csv")

##  calculated from the EC measurements
inBA_LG$es <- 0.6113* exp((17.27*inBA_LG$TW)/(inBA_LG$TW + 237.3))                                   # es is in kPa, Tetens formula in Stull 2007
inBA_LG$ea <- 0.01*inBA_LG$RH*0.6113* exp((17.27*inBA_LG$Temp_amb)/(inBA_LG$Temp_amb + 237.3))                       # ea in kPa
inBA_LG$qs <- (0.622*inBA_LG$es) / (inBA_LG$Amb_Press - (0.378*inBA_LG$es))                          # P is in kPa, qs in kg/kg; 0.378 = 1 - 0.622 (Stull
inBA_LG$qa <-  (0.622*inBA_LG$ea) / (inBA_LG$Amb_Press - (0.378*inBA_LG$ea))                         # qa in kg/kg


### Evaporation by BA
rho <- 1.2                                                                     # the air density (in kg m3 ) 
# Heikinheimo et al. / Agricultural and Forest Meteorology, (1999) : Forest lakes
C_3 <- 0.00107                                                                  # Dalton number at 3 m . z_0 <- 0.001                                                                    # z0 roughness length in meters (for water surface)
z_0 <- 0.002
# recalculate the coefficient at 2 m
C_2 <- C_3 * (log(1.8 / z_0) / log(3 / z_0))
inBA_LG$EBA_hk <- rho*C_2*inBA_LG$wind_speed * (inBA_LG$qs - inBA_LG$qa)*30*60   #30*60*E to convert from kg m-2s-1 to mm per 30 min
rm(C_2)

# wind depended coefficient (values from the experiment on Lake Zub/Priyadarshini): Annex A, Timo Vihma 25.8.2024: Material on lake evaporation for EMS 2024
inBA_LG$CEz_wd <- ifelse(inBA_LG$wind_speed < 13, 0.0000119*inBA_LG$wind_speed+0.0014, 0.0018)
inBA_LG$EBA_wd <- rho*inBA_LG$CEz_wd*inBA_LG$wind_speed * (inBA_LG$qs - inBA_LG$qa)*30*60   #30*60*E to convert from kg m-2s-1 to mm per 30 min

# after Arya and Fedorovich: Miguel Potes ...Annex A z=2m
CEz <- 0.001166                                                                 #the turbulent transfer coefficient for moisture (unit less). 
inBA_LG$EBA_af <- rho*CEz*inBA_LG$wind_speed * (inBA_LG$qs - inBA_LG$qa)*30*60  #30*60*E to convert from kg m-2s-1 to mm per 30 min
rm(CEz)

# after  Andreas: Miguel Potes ...Annex A z=2
CEz <- 0.001676                                                                 # 
inBA_LG$EBA_an <- rho*CEz*inBA_LG$wind_speed * (inBA_LG$qs - inBA_LG$qa)*30*60  
rm(CEz, z_0)
##################################################################
### converting to daily data'
eva_by_day <-  to_daily_evap_BA(inBA_LG)
# Annex C: the bulk-aerodynamic, Lake Glubokoe
# output_file_name <- "/home/shevnina/Manuscripts/2025/Ant/september/result_tab/BA_daily_GL.csv"
#write.csv(eva_by_day,output_file_name)
BA_daily_LG <- eva_by_day
rm(eva_by_day)
rm(C_3)

# to merge with the reference evaporation (EC) data
nd <- merge(dEVap_LG, BA_daily_LG, by="newTS", all.x=TRUE)
BA_daily_LG <- nd
rm(nd)


#################################################################################
### Lake Glubokoe
table3_LG <- get_stats_BA(BA_daily_LG)                                         # summary statistics of the BA
table4_LG<-to_scores_BA_CP(BA_daily_LG)                                       # the scores of the BA


######################################################################################################################
################### Lake ZUB: BA
inBA_ZB <-read.csv("/home/shevnina/Manuscripts/2025/Ant/september/result_tab/BA_ZB.csv")
inBA_ZB$Timestamp_UTC <- as.POSIXct(inBA_ZB$Timestamp_UTC, format="%Y-%m-%d %H:%M:%S", tz="UTC")

ECH_ZB<-read.csv("/home/shevnina/Manuscripts/2025/Ant/september/result_tab/ECH_ZB.csv")
# to read the reference daily evaporation (EC) data
dEVap_ZB <- read.csv("/home/shevnina/Manuscripts/2025/Ant/september/result_tab/evap_daily_ZB.csv")

# to calculate qs-qa
inBA_ZB$es <- 0.6113* exp((17.27*inBA_ZB$TW)/(inBA_ZB$TW + 237.3))                                   # es is in kPa, Tetens formula in Stull 2007
inBA_ZB$ea <- 0.01*inBA_ZB$RH*0.6113* exp((17.27*inBA_ZB$Temp_amb)/(inBA_ZB$Temp_amb + 237.3))       # ea in kPa
inBA_ZB$qs <- (0.622*inBA_ZB$es) / (inBA_ZB$Amb_Press - (0.378*inBA_ZB$es))                          # P is in kPa, qs in kg/kg
inBA_ZB$qa <- (0.622*inBA_ZB$ea) / (inBA_ZB$Amb_Press - (0.378*inBA_ZB$ea))                          # qa in kg/kg

### parametrization of CEz
# Heikinheimo et al. / Agricultural and Forest Meteorology, (1999) : Forest lakes
########### contribution by T. Vihma in Shevnina et al., 2022 ##########################################################################################
path<-"/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/"
setwd(path)
eva_bulk<-read.table("Bulk_method_results_Irgason_input.txt", header = T, sep=" ")  ### Irgason input data
#colnames(eva_bulk)
#"T2m"
# "Ts"
# "V2m"
#"H2O_conc"
# "LE"
# "E_mm_day"
# "H"
# "TS"

# R
eva_bulk$Timestamp_UTC <- as.POSIXct(
  paste(
    2018,
    eva_bulk$Month,
    eva_bulk$Day,
    sprintf("%02d:%02d:00", eva_bulk$Hour, eva_bulk$Minute)
  ),
  format = "%Y %m %d %H:%M:%S",
  tz = "UTC"
)

#
merged_ZB <- merge(inBA_ZB, eva_bulk, by = "Timestamp_UTC", all.x = TRUE)
inBA_ZB <- merged_ZB
inBA_ZB$EBA_hk <- inBA_ZB$E_mm_day/48     # 2 * 24 ... to convert in mm 30 min-1

# wind depended coefficient (values from the experiment on Lake Zub/Priyadarshini): Annex A, Timo Vihma 25.8.2024: Material on lake evaporation for EMS 2024
inBA_ZB$CEz_wd <- ifelse(inBA_ZB$wind_speed < 13, 0.0000119*inBA_ZB$wind_speed+0.0014, 0.0018)
inBA_ZB$EBA_wd <- rho*inBA_ZB$CEz_wd*inBA_ZB$wind_speed * (inBA_ZB$qs - inBA_ZB$qa)*30*60   #30*60*E to convert from kg m-2s-1 to mm per 30 min

# CEz is the turbulent transfer coefficient for moisture (unit less).
# after Arya and Fedorovich: Miguel Potes ...Annex A z= 1.8m
# average from the experiment 2019-2020
#CEz <- 0.001166
C_2 <- 0.001166                                                                   # Dalton number at 1.8 m . z_0 <- 0.002                                                                    # z0 roughness length in meters (for water surface)
z_0 <- 0.002
# recalculate the coefficient at 2 m
CEz <- C_2 * (log(2/ z_0) / log(1.8 / z_0))

inBA_ZB$EBA_af <- rho*CEz*inBA_ZB$wind_speed * (inBA_ZB$qs - inBA_ZB$qa)*30*60  #30*60*E to convert from kg m-2s-1 to mm per 30 min
rm(CEz, C_2,z_0)

# after  Andreas: Miguel Potes ...Annex A z=1.8 m
#CEz <- 0.001676
C_2 <- 0.001676                                                                   # Dalton number at 1.8 m . z_0 <- 0.002                                                                    # z0 roughness length in meters (for water surface)
z_0 <- 0.002
# recalculate the coefficient at 2 m
CEz <- C_2 * (log(2/ z_0) / log(1.8 / z_0))
# average from the experiment 2019-2020
inBA_ZB$EBA_an <- rho*CEz*inBA_ZB$wind_speed * (inBA_ZB$qs - inBA_ZB$qa)*30*60
rm(CEz)
rm(C_2, z_0)


### converting to daily data
eva_by_day <-  to_daily_evap_BA(inBA_ZB)
# Annex C: the bulk-aerodynamic, Lake ZUB
write.csv(eva_by_day,"/home/shevnina/Manuscripts/2025/Ant/september/result_tab/BA_daily_ZB.csv")
BA_daily_ZB <- eva_by_day
rm(eva_by_day)

# to merge with the reference evaporation (EC) data
nd <- merge(dEVap_ZB, BA_daily_ZB, by="newTS", all.x=TRUE)
BA_daily_ZB <- nd
rm(nd)
#################################################################################


#################################################################################
### Lake Zub
table3_ZB <- get_stats_BA(BA_daily_ZB)                                         # summary statistics of the BA
table4_ZB<-to_scores_BA_CP(BA_daily_ZB)                                       # the scores of the BA
#####################################################################################################################
######################################### EOF Bulk-aerodynamic method #################################################



################################ METHODS
################################ Combination Equations ################################################################
################################# Functions ###########################################################################
################ Daily aggregation ##########################################################################################
to_daily_meteo_EC <- function(dat)
 {
  d24h <- aggregate(list(TW = dat$TW, AT = dat$Temp_amb, RH = dat$RH, WS = dat$wind_speed, e2 = dat$ea, es = dat$es), list(Timestamp= cut(dat$Timestamp_UTC, "1440 min" )), mean, na.rm=TRUE)
  EEC <-aggregate(list(EEC_mean = dat$Evap_filter),list(Timestamp= cut(dat$Timestamp_UTC,"1440 min" )), sum, na.rm=TRUE)
  d24h$newTS <- as.Date(d24h$Timestamp, "%Y-%m-%d")
  d24h$EEC <- EEC$EEC_mean
  d24h$DSP <- (d24h$es-d24h$e2)                                                # DSP is (es-e2) saturation deficit => kPa
  return(d24h)
}
# format of dat: output for BA 30 min data
# format of d24h: newTS, AT,TW,WS,RH EEC (direct estimate of the evaporation)

#A is the surface area in square metres (m2 )
evap_CC <- function (d24h, A)
{
  d24h$eva_pe <- 0.26*(1+0.54*d24h$WS)*d24h$DSP                             ### penmann, 1948
  d24h$eva_od <- 0.14*(1+0.72*d24h$WS)*d24h$DSP                             ### doorenboss and Pruit, 1975
  d24h$eva_dp <- 0.26*(1+0.86*d24h$WS)*d24h$DSP                             ### odrova, 1979
  d24h$eva_sw <- 2.909*A^(-0.05)*d24h$DSP                                   ### Shuttleworth (1993)
  d24h$eva_sp <- -0.33*(1-1.82*d24h$WS)*d24h$DSP                            ### shevnina, potes etal., 2022
  d24h$eva_ec <- d24h$EEC                                                   ### EC data, reference
  return(d24h)
}
##

# combination equations CE
to_scores_CE <- function (eva_day)
{
### Rename columns
  names(eva_day)[names(eva_day) == "eva_pe"] <- "eva_A"
  names(eva_day)[names(eva_day) == "eva_dp"] <- "eva_B"
  names(eva_day)[names(eva_day) == "eva_od"] <- "eva_C"
  names(eva_day)[names(eva_day) == "eva_sp"] <- "eva_D"
  names(eva_day)[names(eva_day) == "eva_sw"] <- "eva_E"
  names(eva_day)[names(eva_day) == "EEC"] <- "eva_ec"

#### RMSE == Root Mean Square Error (RMSE)  ###########################################################
  m_A <- sqrt(mean((eva_day$eva_ec - eva_day$eva_A)^2, na.rm=TRUE))                  ## ODROVA
  m_B <- sqrt(mean((eva_day$eva_ec - eva_day$eva_B)^2, na.rm=TRUE))                  ## PENMANN
  m_C <- sqrt(mean((eva_day$eva_ec - eva_day$eva_C)^2, na.rm=TRUE))                  ## DOORENBOSS
  m_D <- sqrt(mean((eva_day$eva_ec - eva_day$eva_D)^2, na.rm=TRUE))                  ## SP
  m_E <- sqrt(mean((eva_day$eva_ec - eva_day$eva_E)^2, na.rm=TRUE))                  ## Shuttleworth (1993)

#### s/sigma criterio, (Popov,  1979, p. 58) #############################################
  m1 <- mean(eva_day$eva_ec, na.rm=T)                                                          # mean evaporation, EC method
  l1 <- length(eva_day$eva_ec)                                                        # length of time series
  sigma <- sqrt(sum((eva_day$eva_ec-m1)^2, na.rm=T)/l1)

  s<-sqrt(sum((eva_day$eva_D-eva_day$eva_ec)^2, na.rm=T)/(l1-2))                             ## SP
  ssg_D <-s/sigma

  s<-sqrt(sum((eva_day$eva_A-eva_day$eva_ec)^2, na.rm=T)/(l1-2) )                           ## Odrova ### 2.37
  ssg_A <- s/sigma

  s<-sqrt(sum((eva_day$eva_B-eva_day$eva_ec)^2, na.rm=T)/(l1-2))                            ## Penmann  ### 2.20
  ssg_B <- s/sigma

  s<-sqrt(sum((eva_day$eva_C-eva_day$eva_ec)^2, na.rm=T)/(l1-2))                            ## Doorenboss  ### 2.09
  ssg_C <- s/sigma

  s<-sqrt(sum((eva_day$eva_E-eva_day$eva_ec)^2, na.rm=T)/(l1-2))                            ## Shuttleworth (1993)  ###
  ssg_E <- s/sigma

#### Pearson correlation coefficient ###########################################################
  c_A <- cor(eva_day$eva_ec,eva_day$eva_A, method="pearson",use="complete.obs")                         #
  c_C <- cor(eva_day$eva_ec,eva_day$eva_C, method="pearson",use="complete.obs")                         #
  c_B <- cor(eva_day$eva_ec,eva_day$eva_B, method="pearson",use="complete.obs")                         #
  c_D <- cor(eva_day$eva_ec, eva_day$eva_D, method="pearson",use="complete.obs")                        #
  c_E <- cor(eva_day$eva_ec, eva_day$eva_E, method="pearson",use="complete.obs")                        #

  tableCE <- data.frame(names = c("OD","PE","DB","SP","SW"), RMSE = c(m_A,m_B,m_C,m_D,m_E), SSg = c(ssg_A,ssg_B,ssg_C,ssg_D, ssg_E), Pe = c(c_A,c_B,c_C,c_D,c_E), stringsAsFactors = FALSE)
  return(tableCE)
}

#########################################################################################
get_summaryCE <- function(CE_evap_GL)
  {
s_ec <- sum(CE_evap_GL$EEC,na.rm=T)
s_pe <- sum(CE_evap_GL$eva_pe,na.rm=T)
 r_pe <- s_ec/sum(CE_evap_GL$eva_pe,na.rm=T)
s_dp <- sum(CE_evap_GL$eva_dp,na.rm=T)
 r_dp <- s_ec/sum(CE_evap_GL$eva_dp,na.rm=T)
s_od <- sum(CE_evap_GL$eva_od,na.rm=T)
 r_od <- s_ec/sum(CE_evap_GL$eva_od,na.rm=T)
s_sw <- sum(CE_evap_GL$eva_sw,na.rm=T)
 r_sw <- s_ec/sum(CE_evap_GL$eva_sw,na.rm=T)
s_sp <- sum(CE_evap_GL$eva_sp,na.rm=T)
 r_sp <- s_ec/sum(CE_evap_GL$eva_sp,na.rm=T)

ax_pe <- max(CE_evap_GL$eva_pe,na.rm=T)
ax_dp <- max(CE_evap_GL$eva_dp,na.rm=T)
ax_od <- max(CE_evap_GL$eva_od,na.rm=T)
ax_sw <- max(CE_evap_GL$eva_sw,na.rm=T)
ax_sp <- max(CE_evap_GL$eva_sp,na.rm=T)

an_pe <- mean(CE_evap_GL$eva_pe,na.rm=T)
an_dp <- mean(CE_evap_GL$eva_dp,na.rm=T)
an_od <- mean(CE_evap_GL$eva_od,na.rm=T)
an_sw <- mean(CE_evap_GL$eva_sw,na.rm=T)
an_sp <- mean(CE_evap_GL$eva_sp,na.rm=T)

tableCE <- data.frame(names = c("PE","DB","OD","SW","SP"), SUM_CE = c(s_pe,s_dp,s_od,s_sw,s_sp), R_CE = c(r_pe,r_dp,r_od,r_sw, r_sp), MEAN_CE = c(an_pe,an_dp,an_od,an_sw, an_sp), MAX_CE = c(ax_pe,ax_dp,ax_od,ax_sw, ax_sp), stringsAsFactors = FALSE)

return(tableCE)
}
########################### EOF Functions ###########################################################

#####################################################################################################
# Lake Glubokoe
#### Lake Glubokoe
inBA_LG <- read.csv("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/result_tab/BA_LG.csv")
inBA_LG$Timestamp_UTC <- as.POSIXct(inBA_LG$Timestamp_UTC, format="%Y-%m-%d %H:%M:%S", tz="UTC")

##  calculated from the EC measurements
inBA_LG$es <- 0.611* exp((12.27*inBA_LG$TW)/(inBA_LG$TW + 273.15))*10                                # es is in kPa, Tetens formula in Stull 2007
inBA_LG$ea <- 0.611* exp((12.27*inBA_LG$Temp_amb)/(inBA_LG$Temp_amb + 273.15))*10                    # ea in kPa

#inBA_LG$es <- 0.6112* exp((17.27*inBA_LG$TW)/(inBA_LG$TW + 243.5))                                   # es is in kPa, WMO guide, 2008, p. I.4–29
#inBA_LG$ea <- 0.01*inBA_LG$RH*0.6112* exp((17.27*inBA_LG$Temp_amb)/(inBA_LG$Temp_amb + 243.5))       # ea in kPa

# to aggregate to
CE_daily_LG <- to_daily_meteo_EC(inBA_LG)

# to calculate the evaporation by the combination equations
A <- 147000
CE_evap_GL <- evap_CC(CE_daily_LG, A)

# to get summery of the CE
CE_sum_LG <- get_summaryCE(CE_evap_GL)
write.csv(CE_sum_LG, "/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/result_tab/CE_summary_GL.csv")

eva_day <- CE_evap_GL
CE_scope_LG <- to_scores_CE(eva_day)
write.csv(CE_scope_LG, "/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/supplement/result_tab/CE_scope_GL.csv")
##################### EOF of the combination equations, Lake Glubokoe ##############################################


############### FIGURES #################################################################
### Figure 11 b, LAKE GLUBOKOE
eva_day <- inBA_LG
# Rename columns
names(eva_day)[names(eva_day) == "EBA_hk"] <- "eva_A"
names(eva_day)[names(eva_day) == "EBA_wd"] <- "eva_B"
names(eva_day)[names(eva_day) == "EBA_af"] <- "eva_C"
names(eva_day)[names(eva_day) == "EBA_an"] <- "eva_D"
names(eva_day)[names(eva_day) == "Evap"] <- "eva_ec"
#str(eva_day)

# QQ EC vs BA
pdf("/home/shevnina/Manuscripts/2025/Ant/september/result_fig/ECvsBA_QQ_GL.pdf")                                     # Figure 11 in the manuscript
par(family = "serif", cex = 1.0)
plot(eva_day$eva_ec,eva_day$eva_A, xlim=c(-0.02,0.2), ylim=c(-0.02,0.2), pch=3, xlab="E_EC, mm 30m-1", ylab="E_BA, mm 30m-1", col="blue")
#points(eva_day$eva_ec,eva_day$eva_D,col="grey", pch=3)
#points(eva_day$eva_ec,eva_day$eva_C, pch=3, col = "orange")
points(eva_day$eva_ec,eva_day$eva_B,col="red", pch=3)
grid()
abline(a = 0, b = 1, col="black", lwd = 0.5)
abline(lm(eva_day$eva_B ~ eva_day$eva_ec), col="red", lwd=0.5)
abline(lm(eva_day$eva_A ~ eva_day$eva_ec), col="blue", lwd=0.5)
#abline(lm(eva_day$eva_C ~ eva_day$eva_ec), col="orange", lwd=0.5)
#abline(lm(eva_day$eva_D ~ eva_day$eva_ec), col="grey", lwd=0.5)
#legend("topleft", legend=c("HK", "WD","AF","AN"), col=c("blue", "red","orange","grey"), pch=3)
legend("topleft", legend=c("HK", "WD"), col=c("blue", "red"), pch=3)
dev.off()

### Figure 11a, LAKE ZUB
eva_day <- inBA_ZB
# Rename columns
names(eva_day)[names(eva_day) == "EBA_hk"] <- "eva_A"
names(eva_day)[names(eva_day) == "EBA_wd"] <- "eva_B"
names(eva_day)[names(eva_day) == "EBA_af"] <- "eva_C"
names(eva_day)[names(eva_day) == "EBA_an"] <- "eva_D"
names(eva_day)[names(eva_day) == "Evap"] <- "eva_ec"
#str(eva_day)

# QQ EC vs BA
pdf("/home/shevnina/Manuscripts/2025/Ant/september/result_fig/ECvsBA_QQ_ZB.pdf")
par(family = "serif", cex = 1.0)
plot(eva_day$eva_ec,eva_day$eva_A, xlim=c(-0.02,0.2), ylim=c(-0.02,0.2), pch=3, xlab="E_EC, mm 30m-1", ylab="E_BA, mm 30m-1", col="blue")
points(eva_day$eva_ec,eva_day$eva_D,col="grey", pch=3)
points(eva_day$eva_ec,eva_day$eva_C, pch=3, col = "orange")
#points(eva_day$eva_ec,eva_day$eva_B,col="red", pch=3)
grid()
abline(a = 0, b = 1, col="black", lwd = 0.5)
#abline(lm(eva_day$eva_B ~ eva_day$eva_ec), col="red", lwd=0.5)
abline(lm(eva_day$eva_A ~ eva_day$eva_ec), col="blue", lwd=0.5)
abline(lm(eva_day$eva_C ~ eva_day$eva_ec), col="orange", lwd=0.5)
abline(lm(eva_day$eva_D ~ eva_day$eva_ec), col="grey", lwd=0.5)
# Add a legend
#legend("topleft", legend=c("HK", "WD","AF","AN"), col=c("blue", "red","orange","grey"), pch=3)
legend("topleft", legend=c("HK", "AF","AN"), col=c("blue", "orange", "grey"), pch=3)
dev.off()
###############################################################################################################################


########################################################################################################################
# Figure XX Daily evaporation 
# Lake Glubokoe  
eva_day <- eva_GL_daily
# Rename columns
names(eva_day)[names(eva_day) == "eva_pe"] <- "eva_A"
names(eva_day)[names(eva_day) == "eva_dp"] <- "eva_B"
names(eva_day)[names(eva_day) == "eva_od"] <- "eva_C"
names(eva_day)[names(eva_day) == "eva_sp"] <- "eva_D"

pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/EVA_CEq_TSplot_GL.pdf", width = 6, height = 5)
par(family = "serif", cex = 1.0)
par(xpd = TRUE)
par(mar = par()$mar + c(0, 0, 0, 5.5))
plot(as.Date(eva_day$newTS), eva_day$eva_ec, type="p",pch = 20, cex=0.7, col="black", xlim=c(as.Date("2019-12-01"), as.Date("2020-01-28")), ylim=c(-1,6), ylab="E, mm/d", xlab="Date")
points(as.Date(eva_day$newTS), eva_day$eva_C, col="red", pch = "+", cex=0.5)                 # Odrova
points(as.Date(eva_day$newTS), eva_day$eva_B, col="blue", pch = "+", cex=0.5)                         #DP
points(as.Date(eva_day$newTS),eva_day$eva_D, col="orange", pch = "+", cex=0.5)                        #SP.22
points(as.Date(eva_day$newTS), eva_day$eva_A, col="grey", pch = "+", cex=0.5)                         #penmann

lines(as.Date(eva_day$newTS), eva_day$eva_ec, col="black",lty=2, lwd=0.5) 
lines(as.Date(eva_day$newTS), eva_day$eva_C, col="red",lty=2, lwd=0.5)                           # Odrova
lines(as.Date(eva_day$newTS), eva_day$eva_B, col="blue", lty=2, lwd=0.5)                         #DP
lines(as.Date(eva_day$newTS),eva_day$eva_D, col="orange", lty=2, lwd=0.5)                        #SP.22
lines(as.Date(eva_day$newTS), eva_day$eva_A, col="grey", lty=2, lwd=0.5)                         #penmann
abline(h=0)
abline(v=as.Date("2019-12-11"), col="black",lty=2, lwd=0.3)                                     # Air temp is positive until this day
grid()
#legend("topright", inset = c(-0.2, 0), legend=c("EC", "OD","DP","SP","PE"), col=c("black","red", "blue","orange","grey"),lwd=1)
dev.off()

### LAKE GLUBOKOE
eva_day <- eva_GL_daily
# Rename columns
names(eva_day)[names(eva_day) == "eva_pe"] <- "eva_A"
names(eva_day)[names(eva_day) == "eva_dp"] <- "eva_B"
names(eva_day)[names(eva_day) == "eva_od"] <- "eva_C"
names(eva_day)[names(eva_day) == "eva_sp"] <- "eva_D"
str(eva_day)

#!Figure QQplot
# EC + empiric
pdf("/home/shevnina/Manuscripts/2024/ant/Ant_Evaporation/figures/ECvsComEq_QQplot_GL.pdf")
par(family = "serif", cex = 1.0)
plot(eva_day$eva_ec,eva_day$eva_A, xlim=c(0,6), ylim=c(0,6), pch=3, xlab="E_EC, mm/d", ylab="E_CE, mm/d", col="blue") # Water
points(eva_day$eva_ec,eva_day$eva_D,col="grey", pch=3)                                             
points(eva_day$eva_ec,eva_day$eva_C, pch=3, col = "orange")    
points(eva_day$eva_ec,eva_day$eva_B,col="red", pch=3)
grid()

abline(a = 0, b = 1, col="black", lwd=0.1)
#abline(lm(eva_day$eva_B~eva_day$eva_ec), col="red")
#abline(lm(eva_day$eva_A~eva_day$eva_ec), col="blue")
#abline(lm(eva_day$eva_C~eva_day$eva_ec), col="orange")
#abline(lm(eva_day$eva_D~eva_day$eva_ec), col="grey")

# Add a legend
legend("topleft", legend=c("A", "B","C","D"), col=c("blue", "red","orange","grey"), pch=3)
dev.off()
###############################################################################################################################
