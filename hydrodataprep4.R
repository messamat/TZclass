#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/26/2018
#Date last updated: 03/28/2018

#Purpose: impute/infill missing discharge data for all stream gauges

library(ggplot2)
library(data.table)
library(FlowScreen)
library(hydroTSM)
library(plyr)
library(zoo)
library(lemon)
library(foreign)
library(forecast)
library(imputeTS)
setwd("F:/Tanzania/Tanzania/results") #UPDATE
datadir = file.path(getwd(),paste('rufiji_hydrodatainspect','20180326',sep='_')) #UPDATE
origdatadir = "F:/Tanzania/Tanzania/data"
outdir=file.path(getwd(),paste('rufiji_hydrodataimpute',as.character(format(Sys.Date(),'%Y%m%d')),sep='_'))
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}
rufidat_clean <- read.csv(file.path(datadir,'rufidat_clean.csv'), colClasses=c('factor','Date','numeric','character','character'))

################################################
#Infilling/Interpolate data
################################################
rufidat_clean$Flowlog <- log(rufidat_clean$Flow+0.01)
rufidat_cast <- dcast(rufidat_clean, Date~ID, value.var='Flowlog')

#Try auto.arima from forecast package and KalmanSmoother, inspired from  https://stats.stackexchange.com/questions/104565/how-to-use-auto-arima-to-impute-missing-values
#Function to find the index, for each row, of the previous row with a non-NA value
#Takes in a univariate time series with NAs
na.lomf <- function(x) {
  if (length(x) > 0) {
    non.na.idx <- which(!is.na(x))
    rep.int(non.na.idx, diff(c(non.na.idx, length(x) + 1))) #Repeat index of previous row with non-NA as many times gap until next non-NA values
  }
}
#Function to find the index, for each row, of the next row with a non-NA value
#Takes in a univariate time series with NAs
na.lomb <- function(x) {
  if (length(x) > 0) {
    non.na.idx <- which(!is.na(x))
    y<-rep(NA,length(x))
    y[non.na.idx] <- non.na.idx
    zoo::na.locf(y, fromLast = TRUE)
  }
}
#Function to fill gaps for one stream gauge based on an ARIMAX model using time series characteristics from that gauge and the other 
#gages as exogenous regressors.
#tscast data frame where the first is called 'Date' and contains dates and all the other columns are time flow series
#sn: column of the time series to be infilled
#maxgap: maximum gap length to be imputed
CustomImpute <- function(tscast, sn, maxgap) {
  name<- colnames(tscast)[sn]
  print(name)
  #Restrict analysis to min and max year of records
  mindate <- tscast[min(which(!is.na(tscast[,sn]))),'Date']
  maxdate <- tscast[max(which(!is.na(tscast[,sn]))),'Date']
  tscastsub <- tscast[tscast$Date>mindate & tscast$Date<maxdate,]
  #Compute size of gap a record belongs to
  tscastsub$prevdate <- tscastsub[na.lomf(tscastsub[,sn]),'Date']
  tscastsub$nextdate <- tscastsub[na.lomb(tscastsub[,sn]),'Date']
  tscastsub$gap <- as.numeric(as.Date(tscastsub$nextdate)-as.Date(tscastsub$prevdate))
  tscastsub <- as.data.frame(tscastsub)
  #Fit an ARIMAX model
  pred <- tscastsub[,sn]
  fit <- auto.arima(tscastsub[,sn],xreg=tscastsub[,-c(1,sn,ncol(tscastsub)-c(0,1,2))])
  #summary(fit)
  id.na <- which((tscastsub$prevgap<=maxgap &
                    is.na(tscastsub[,sn])))
  #Kalman smoother
  kr <- KalmanSmooth(tscastsub[,sn], fit$model)
  for (i in id.na)
    pred[i] <- fit$model$Z %*% kr$smooth[i,]
  #Output observed and predicted predicted data
  df_pred <- data.frame(Date=tscastsub[id.na,'Date'],pred=pred[id.na])
}

CustomImpute

ggplot(rufidat_castsub, aes(x=Date, y=get(colnames(rufidat_castsub)[sn])))+
  geom_point(data=df_pred,aes(y=pred), color='red') + 
  geom_point(color='black') +
  scale_x_date(limits=as.Date(c('1965-01-01','1975-01-01')))


# #Try imputeTS package Kalman filter
# sn='1KA8A'
# mindate <- rufidat_castsub[min(which(!is.na(rufidat_castsub[,sn]))),'Date']
# x <- na.kalman(rufidat_castsub[rufidat_castsub$Date>mindate,sn], model = "auto.arima")
# rufidat_castsub[rufidat_castsub$Date>mindate,'pred1KA8A'] <- x
# ggplot(rufidat_castsub[rufidat_castsub$Date>mindate,], aes(x=Date, y=get(sn)))+
#   geom_point(aes(y=pred1KA8A), color='red') +
#   geom_point(color='black')+
#   scale_x_date(limits=as.Date(c('1965-01-01','1975-01-01')))