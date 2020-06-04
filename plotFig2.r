#!/usr/bin/Rscript

# The supplement to Shevnina et al., 2020
# Fig. 2. Wind direction and frequency of wind speed anomalies (according to 6 hour synoptic records of Novo meteo site).
args<-commandArgs(TRUE)

# To read table with data from file: 512n_wind_dec.txt or 89512n_wind_jan.txt (BAS, Reader old).
wdspt <- read.table(args[1], header=T, sep=" ")

# To calculate the mean wind speed (ms-1) 
m_wsp<-mean(wdspt$wsd)*0.5144

# To caculate anomaly of wind speed
wdspt$awsp<-c(wdspt$wsd*0.5144-m_wsp)

# to reclass of wind speed anomalies
#Levels: (-100,-10] (-10,-5] (-5,0] (0,5] (5,10] (10,15] (15,23] (23,100]
#Level_labels: 1 2 3 4 5 6 7 8
wdspt$lawsp<-cut(wdspt$awsp, breaks=c(-100,-10,-5,0,5,10,15,23,100))

# to reclass of wind direction
#Levels: [0,10] (10,20] (20,30] (30,40] (40,50] (50,60] (60,70] ... (350,360]
wdbreaks<-seq(0,360, by=10)
#Level_labels:10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 ... 360
wdlabels<-seq(10,360, by=10)
wdspt$lwd<-cut(wdspt$wd,breaks=wdbreaks,labels=c(wdlabels),include.lowest=T)

# to calculate the table for the figure: x-value wd, y-value class of aspt
gdata<-table(wdspt$lawsp,wdspt$lwd)

# to save plot in external file: name for the output file
pdf(args[2])
barplot(gdata, xlab='Wind Direction',ylab='Number of cases',col=c('grey','blue','lightblue','orange','yellow','red','darkred','black'),legend=rownames(gdata))
dev.off()


