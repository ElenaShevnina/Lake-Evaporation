# Evaporation over lakes (direct and indirect methods)
Summertime evaporation over two glacial lakes in the Schimacher oasis (SA), East Antarctica (Lake Priyadarshini/Zub and Lake Glubokoe) was measured by the eddy-covariance technique in December – February 2017-2018 and 2019-2020. The evaporation was also calculated applying the indirect (mass transfer) methods.

## Data 
20180101_20180207_EC_FLUX.txt and 20191207_20200108_EC_FLUX.txt are the raw data on the fluxes collected by the EC system (Irgason by Campbell Scientific) installed on the shore on the lakes. The EC system was also equipped with air temperature and barometric pressure auxiliary sensors (Betatherm 100K6A1A Thermistor and Vaisala PTB110). Lake Surface Water Temperature (LSWT) was measured by the HOBO sensor, and HOBO20172018.CSV and HOBO20192020.CSV contain the raw data collected in two experiments.

## Methods
E_direct_method.R contains the code processing the EC fluxes following (Potes et al., 2017). E_indirect_methods.R codes the processing of the daily evaporation with the bulk-aerodynamic methods and combination formulas (Shevnina et 2022).  
