#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/22/2018
#Date last updated: 03/22/2018

#Purpose: compare data provided by CDMSmith Zachary T. Eichenwald (ZTE) and Japhet Kashaigili (JK)

library(ggplot2)
setwd("F:/Tanzania/Tanzania/results") #UPDATE
datadir = file.path(getwd(),paste('rufiji_hydrodataraw','20180322',sep='_')) #UPDATE

setClass('myDate')
setAs("character","myDate", function(from)  as.POSIXct(from, format= "%m/%d/%Y %H:%M"))
rufidat <- read.csv(file.path(datadir,'ZTE_rufidat.csv'),colClasses = c('factor','factor','myDate','numeric','numeric','factor','factor','factor'))
str(rufidat)
Japhet_dailydat <- read.csv(file.path(datadir,'JK_dailydat.csv'),colClasses = c('Date','numeric','factor'))
Japhet_Unimpaired <- read.csv(file.path(datadir,'JK_monthly_unimpaired.csv'),colClasses = c('Date','factor','numeric'))
Japhet_Rukwa <- read.csv(file.path(datadir,'JK_monthly_rukwa.csv'),colClasses = c('Date','factor','numeric'))
