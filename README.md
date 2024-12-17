# Lake Mass/Water balance: two glacial lakes in the Schimacher oasis (SA), East Antarctica

# the evaporation (day-by-day, seasonal)

Lake Priyadarshini / Zub (2017-2018) and Lake Glubokoe (2019-2020)

## List of the supplemens 
### The eddy covariance: experiments in 2018 and in 2019 (M. Potes)
20180101_20180207_EC_FLUX.txt and 20191207_20200108_EC_FLUX.txt are the raw data on the fluxes collected during two austral summers (Dec-Feb) 2017-2018 and 2019-2020 with the EC system instaled on the shore on the lakes in the Schirmacher oasis, East Antarctica. This data set includes measurements collected by the integrated CO2 and H2O open-path gas analyzer and 3-D sonic anemometer (Irgason by Campbell Scientific. The Irgason is measured an absolute carbon dioxide and water vapor densities, three-dimensional wind speed, sonic air temperature, air temperature and barometric pressure (auxiliary sensor). The auxiliary sensors are Vaisala PTB110 as an optional enhanced barometer, and temperature sensor by Beta Therm 100K6A1A Thermistor. The post-processing procedure is described in Potes et al., (2017).
List of the variables in the dataset:
Date Time, is the time stamp (TS) with the format “yy/mm/dd HH:MM:SS”
npoints, is number of measurements
u_starr (m/s), taur Kg(m s^2) are u_star and momentum flux
Hcr (W/m^2) is a sensible heat flux humidity 
LE_wplr (W/m^2), is a latent heat flux 
Evap (L/m^2) is an evaporation by Potes M.
Fc_wplr (mg/(m^2 s)), Fc_wpl_umol (umol/(m^2 s)), is carbon flux with WPL correction, mg/(m^2 s) and umol/(m^2 s);
CO2_conc (mg/m^3), CO2_PPM (PPM), CO2_sig_str (arb), are CO2 concentration (two unit) and the Irgason’s signal strength; 
H2O_conc (g/m^3), H2O_sig_str (arb), are H2O concentration and the Irgason’s signal strength; 
Temp_amb (C), Amb_Press (kPa), are ambient air temperature and atmospheric pressure;
wind_speed (m/s), wind dir_sonic (degree), are wind speed and direction;
obukhov (m), zeta (arb), rou (m), sigmaw_ustar (arb), are Obuhkov length, stability parameter, roughness;
Xmax (m), XR90 (m), are 90% of footprint (m) and maximum of footprint (m).


### The bulk mass transfer method (T. Vihma)
MTC.CSV : The data include the daily time series of the air temperature (at, C), wind speed (ws, ms-1), atmosheric pressure (ap, Pa) and H2O concentration (h2o, g sm-3). These variables are measured at height of 2 m, and I have used them to calulate the specific humidity (sp, Pa):
SP = R*H2O*(AT+273.15)/P
where R=287.05 is the universal gas constant for a dry air (J kg-1 K-1); AT is the air temperature (C); H2O is water vapor concentration (g sm-3); P is the barometric pressure, (Pa). Also, the daily series of the lake surface temperature (lwt, C) are included. 
MTC_30min.csv: 30 min series of air temperature, wind speed, specific humidity, h2o concentration.... 

### The empirical equations (E.Shevnina)
20171201_20180210_Novo.txt; 20191201_20200131_Novo.txt are 1 minute meteorological records collected by he Milos station at Novo site
WMO ID 89512. The long term climatology at the site is given here: http://www.nerc-bas.ac.uk/icd/gjma/novol.temps.html
Data format by the columns: number;Dates_Times;Air_temperature;Soil_Temperature;Relative_Humidity;Atmospheric_pressure;Wind_speed;Incoming_solar_radiation
The records are provided by the Arctic Antarctic Research Institute for the period of December–February 2017–2018 and 2019–2020.
Termo7_iButton2141.csv is the measured water temperature on Lake Zub/Priyadarshini, experiment 2017-2018.

### The Maitri site meteorological observations (P. Dhote and P. Thakur)



Bellagamba, A. W., Berkelhammer,  M., Winslow, L., Doran, P. T., Myers,  K. F., Devlin, S. , Hawes, I.: The magnitude and climate sensitivity of isotopic fractionation from ablation of Antarctic Dry Valley lakes, Arctic, Antarctic, and Alpine Research, 53:1, 352-371, DOI: 10.1080/15230430.2021.2001899, 2021. 

Box, J. E., Steffen, K.: Sublimation on the Greenland Ice sheet from automated weather station observations. Journal of Geophysical Research: Atmospheres 106 (D24):33965–81, 2001

Gat, J.R., Gonfiantini, R. (Eds): Stable isotope hydrology, Deuterium and oxygen-18 in the water cycle, International Atomic Energy Agency, Vienna, 356 p.,  1981.

