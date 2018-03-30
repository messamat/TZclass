#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/28/2018
#Date last updated: 03/29/2018

#Purpose: impute/infill missing discharge data for all stream gauges

library(ggplot2)
library(data.table)
library(FlowScreen)
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
rufidat <- read.csv(file.path(getwd(),'rufiji_hydrodataraw_20180324','ZTE_rufidat.csv'),colClasses = c('character','character','Date','numeric','numeric','factor','factor','factor'))
rufidat_clean <- read.csv(file.path(datadir,'rufidat_clean.csv'), colClasses=c('factor','Date','numeric','character','character'))
rufidat_deleted <- read.csv(file.path(datadir,'rufidat_deleted.csv'), colClasses=c('character','Date','numeric','character','character'))
gagesenv <- read.dbf(file.path(getwd(),'gages_netjoin.dbf'))

################################################
#Infilling/Interpolate data
################################################
#rufidat_clean$Flowlog <- log(rufidat_clean$Flow+0.01)
rufidat_cast <- dcast(rufidat_clean, Date~ID, value.var='Flow')

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
#Function to fill gaps for one stream gauge based on an ARIMAX model and Kalman smoother using time series characteristics from that gauge and the other 
#gages as exogenous regressors.
#tscast data frame where the first is called 'Date' and contains dates and all the other columns are time flow series
#sn: column of the time series to be infilled
#maxgap: maximum gap length to be imputed
#Example values to test function:
    #tscast=rufidat_cast
    #sn=2
    #maxgap=30
CustomImpute <- function(tscast, sn, maxgap) {
  if (sn > 1) {
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
    #Only use data from stream gauges that have at least 50% of overlapping non-NA data with that streamgauge
    overlap <- adply(tscastsub[!is.na(tscastsub[,sn]),-c(1,sn,ncol(tscastsub)-c(0,1,2))],2,function(x) {
      length(which(!is.na(x)))/length(which(!is.na(tscastsub[,sn])))
      })
    covar <- which(colnames(tscastsub)%in%overlap[overlap$V1>=0.5,'X1'])
    #Fit an ARIMAX model
    pred <- tscastsub[,sn]
    fit <- auto.arima(tscastsub[,sn],xreg=tscastsub[,covar]) 
    #ARIMAX without forcing seasonality onto the model does not pan out, but ARIMA takes way too long to apply to all gauges
    #Keep very simple interpolation for now, and increase complexity later
    #fit <- auto.arima(ts(tscastsub[,sn],frequency=365),D=1, approximation=F)
    #xreg=ts(tscastsub[,covar], frequency=365)
    #summary(fit)
    id.na <- which((tscastsub$gap<=maxgap &
                      is.na(tscastsub[,sn])))
    #Kalman smoother
    kr <- KalmanSmooth(tscastsub[,sn], fit$model)
    for (i in id.na)
      pred[i] <- fit$model$Z %*% kr$smooth[i,]
    #Output observed and predicted predicted data
    tscastsub[id.na, sn] <- pred[id.na]
    return(tscastsub[,c(1,sn)])
  } else {
    warning('sn must be >1 as 1st column must be "Date"')
  }
}

#Fill in every gauge
impute_preds <- data.frame(Date=rufidat_cast[,"Date"])
for (i in 2:(ncol(rufidat_cast))) {
  try({
  impute_preds<-merge(impute_preds, CustomImpute(rufidat_cast,i,maxgap=30), by='Date', all.x=T)
  })
}
#Append the two gauges for which there was not enough info
impute_preds <- cbind(impute_preds,rufidat_cast[,which(!(colnames(rufidat_cast) %in% colnames(impute_preds)))])
write.csv(impute_preds, file.path(outdir, 'rufidat_interp.csv'), row.names=F)

#############################################################################################################
#Plot clean and interpolated data
predsmelt <-melt(setDT(impute_preds),id.vars = 'Date',value.name='Flow',variable.name='ID')
predsmelt <- predsmelt[,c(2,1,3)]
predsmelt$SYM <- NA
predsmelt$Agency <- NA

for (gage in unique(predsmelt$ID)) {
  print(gage)
  genv <- gagesenv[gagesenv$RGS_No==gage,]
  gname <- paste(genv$RGS_Loc,genv$RGS_Name,sep=" at ")
  #Generate FlowScreen time series
  gts<- create.ts(predsmelt[predsmelt$ID==gage,])  #Cannot run ts on multiple gages. Need to first subset by gage, then run ts.
  #Compute and output flowScreen metrics and plots
  try({
    res <- metrics.all(gts)
    ginfo <- data.frame(StationID=genv$RGS_No, StnName=gname, ProvState='Rufiji Basin',Country='Tanzania',
                        Lat=genv$POINT_Y, Long=genv$POINT_X, Area=genv$WsArea, RHN='RBWB')
    png(file.path(outdir,paste(gage,'screenb.png',sep="_")),width = 20, height=12,units='in',res=300)
    screen.summary(res, type="b", StnInfo=ginfo)
    dev.off()
    png(file.path(outdir,paste(gage,'screenl.png',sep="_")),width = 20, height=12,units='in',res=300)
    screen.summary(res, type="l", StnInfo=ginfo)
    dev.off()
    png(file.path(outdir,paste(gage,'screenh.png',sep="_")),width = 20, height=12,units='in',res=300)
    screen.summary(res, type="h", StnInfo=ginfo)
    dev.off()
  })
  #Make raw time series plot
  rawsgplot <-ggplot(gts, aes(x=Date, y=Flow+0.01)) + 
    geom_point(color='#bf812d', size=1.5) + 
    geom_point(data=rufidat_clean[rufidat_clean$ID==gage,],aes(x=Date, y=Flow+0.01), color='#045a8d', size=1.5) +
    geom_point(data=rufidat_deleted[rufidat_deleted$ID==gage,],aes(x=Date, y=Flow+0.01), color='#e31a1c', size=1.5) +
    scale_y_sqrt(breaks=trans_breaks("sqrt", function(x) x ^ 2))+
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") + 
    labs(y='Discharge (m3/s)', title=paste(gage, gname,sep=" - ")) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  png(file.path(outdir,paste(gage,'raw_sg.png',sep="_")),width = 20, height=12,units='in',res=300)
  print(rawsgplot)
  dev.off()
}


######################################## EXTRA ######################################################
# ggplot(tscastsub, aes(x=Date, y=get(colnames(tscastsub)[sn])))+
#   geom_point(data=try,aes(y=pred), color='red') + 
#   geom_point(color='black') 
# #Try imputeTS package Kalman filter
# sn='1KA8A'
# mindate <- rufidat_castsub[min(which(!is.na(rufidat_castsub[,sn]))),'Date']
# x <- na.kalman(rufidat_castsub[rufidat_castsub$Date>mindate,sn], model = "auto.arima")
# rufidat_castsub[rufidat_castsub$Date>mindate,'pred1KA8A'] <- x
# ggplot(rufidat_castsub[rufidat_castsub$Date>mindate,], aes(x=Date, y=get(sn)))+
#   geom_point(aes(y=pred1KA8A), color='red') +
#   geom_point(color='black')+
#   scale_x_date(limits=as.Date(c('1965-01-01','1975-01-01')))