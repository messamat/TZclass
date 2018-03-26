#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/22/2018
#Date last updated: 03/24/2018

#Purpose: visualize and assess quality and quantity of hydrological data in the Rufiji river basin, provided by CDMSmith Zach T. Eichenwald

library(ggplot2)
library(data.table)
library(FlowScreen)
library(waterData)
library(prospectr)
setwd("F:/Tanzania/Tanzania/results") #UPDATE
datadir = file.path(getwd(),paste('rufiji_hydrodataraw','20180324',sep='_')) #UPDATE
origdatadir = "F:/Tanzania/Tanzania/data"

setClass('myDate')
setAs("character","myDate", function(from)  as.POSIXlt(from, format= "%Y-%m-%d %H:%M:%S"))
rufidat <- read.csv(file.path(datadir,'ZTE_rufidat.csv'),colClasses = c('character','factor','myDate','numeric','numeric','factor','factor','factor'))
str(rufidat)

#General plotting of time series
rawplot <- ggplot(rufidat, aes(x=Date.Time, y=Calculated.Flow..cms.)) + 
  geom_point(color='#045a8d', size=1) + 
  geom_point(data=rufidat[rufidat$Calculated.Flow..cms.==0,],aes(x=Date.Time, y=Calculated.Flow..cms.), color='#e31a1c', size=1) +
  facet_wrap(~Gage.ID+Station.Name, scale='free') +
  theme_bw() + 
  labs(y='Discharge (m3/s)')
#rawplot
png('rufidat_rawts.png',width = 32, height=16,units='in',res=300)
rawplot
dev.off()

#Remove KB33 (below Kihansi, as only contains -999)
rufidat <- data.table(rufidat[rufidat$Gage.ID!='1KB33',])
nrow(rufidat[rufidat$Calculated.Flow..cms.==-999,])

########################################################################################################
#Format rufiji flow data to be used within the FlowScreen package following the USGS data import format
#FlowScreen package was developed to work with Water Survey of Canada (WSC) or the United States Geological Survey (USGS) data. 

#Download and import data from USGS using FlowScreen to see output of 'read.flows function'
testdat <- importDVs('06135000', code = "00060", stat = "00003", sdate = "1851-01-01",
                     edate = as.Date(Sys.Date(), format = "%Y-%m-%d"))
testtab=file.path(getwd(),paste('test_',as.character(format(Sys.Date(),'%Y%m%d')),'.csv',sep=""))
write.csv(testdat, testtab)
test<-read.flows(testtab)
colnames(test)

#Reproduce data structure from read.flows function
str(rufidat)
rufidat_screenform <- data.table(rufidat[,c('Gage.ID','Date.Time','Calculated.Flow..cms.','Record.Quality','Rating.Curve.Source')])
rufidat_screenform$Date.Time <- as.Date(rufidat_screenform$Date.Time)
rufidat_screenform <- rufidat_screenform[,list(Calculated.flow..cms.daily=mean(Calculated.Flow..cms., na.rm=T)), .(Gage.ID, Date.Time,Record.Quality,Rating.Curve.Source)]
rufidat_screenform <- rufidat_screenform[,c('Gage.ID','Date.Time','Calculated.flow..cms.daily','Record.Quality','Rating.Curve.Source')]
colnames(rufidat_screenform) <- colnames(test)
which(duplicated(rufidat_screenform[,c('Date','ID')]))

#Inspect data prior to cleaning
outdir=file.path(getwd(),paste('rufiji_hydrodatainspect',as.character(format(Sys.Date(),'%Y%m%d')),sep='_'))
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}
for (gage in unique(rufidat_screenform$ID)) {
  print(gage)
  #Fit Savitzky-Golay 1st order derivative
  # p = polynomial order w = window size (must be odd) m = m-th derivative (0 = smoothing) 
  gts<- create.ts(rufidat_screenform[rufidat_screenform$ID==gage,])
  d1 <- as.data.frame(savitzkyGolay(gts$Flow, p = 3, w = 5, m = 1))
  gts[3:(nrow(gts)-2),'sg.d1'] <- d1
  gts[gts$sg.d1<0.0001 & gts$sg.d1>-0.0001 & !is.na(gts$sg.d1),'Flag'] <- 'Y'
  
  #Make raw time series plot
  rawplot <-ggplot(gts, aes(x=Date, y=Flow)) + 
    geom_point(color='#045a8d', size=1) + 
    geom_point(data=gts[gts$Flow==0,],aes(x=Date, y=Flow), color='#e31a1c', size=1) +
    theme_bw() + 
    scale_y_sqrt()+
    labs(y='Discharge (m3/s)')
  png(file.path(outdir,paste(gage,'raw.png',sep="_")),width = 8, height=8,units='in',res=300)
  print(rawplot)
  dev.off()
  
  sgplot <-ggplot(gts, aes(x=Date, y=sg.d1)) + geom_point() +
    geom_point(data=gts[gts$Flag=='Y',], color='red') +
    theme_bw()
  png(file.path(outdir,paste(gage,'sg.png',sep="_")),width = 8, height=8,units='in',res=300)
  print(sgplot)
  dev.off()
}

#Notes on gages





##########################################
#Assess general record characteristics
rufidat$year <- format(rufidat$Date.Time, "%Y")
rufidat$month <- as.numeric(format(rufidat$Date.Time, "%m"))







station='1KA15A'
ts<- create.ts(rufidat_screenform[rufidat_screenform$ID==station,]) 
str(ts)
#Initial error because of NA values in the Dates. Had to correct formatting of Kilombero flow dates in hydrodataprep.R
#Then error because tried to use ts for all gages at the same time. Need to first subset by gage, then run ts.
res <- metrics.all(ts)
screen.summary(res, type="b")
screen.summary(res, type="h")
#####







#Visualize time series and qualitative assessment of data by ZTE
#Compute double derivative
#Check maximum diel and daily variation in discharge
#Look for non-stationarity in flow magnitude, timing, and variability, consider de-trending
#Get summary statistics on 
# grain of data
# length of record
# completeness in terms of frequency, length, and time periods of gaps
# overlap in terms of period and length

#Evaluate consistency of discharge with drainage area and precipitation + HydroSHEDS modeled data

#Test range of acceptance criteria (15 years, etc.)
