#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/22/2018
#Date last updated: 03/23/2018

#Purpose: visualize and assess quality and quantity of hydrological data in the Rufiji river basin, provided by CDMSmith Zach T. Eichenwald

library(ggplot2)
library(data.table)
library(FlowScreen)
library(waterData)
setwd("F:/Tanzania/Tanzania/results") #UPDATE
datadir = file.path(getwd(),paste('rufiji_hydrodataraw','20180323',sep='_')) #UPDATE
origdatadir = "F:/Tanzania/Tanzania/data"

setClass('myDate')
setAs("character","myDate", function(from)  as.POSIXlt(from, format= "%Y-%m-%d %H:%M:%S"))
rufidat <- read.csv(file.path(datadir,'ZTE_rufidat.csv'),colClasses = c('character','factor','myDate','numeric','numeric','factor','factor','factor'))

##########################################
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
rufidat_screenform <- rufidatsub <- rufidat_screenform[,list(Calculated.flow..cms.daily=mean(Calculated.Flow..cms.),Record.Quality,Rating.Curve.Source), .(Gage.ID, Date.Time)]
colnames(rufidat_screenform) <- colnames(test)
rufidat_screenform.ts <- create.ts(rufidat_screenform, hyrstart=1) 
#Initial error because of NA values in the Dates. 



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
