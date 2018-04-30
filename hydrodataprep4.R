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
library(hydroTSM)
library(plyr)
library(zoo)
library(lemon)
library(foreign)
library(forecast)
library(imputeTS)

setwd("F:/Tanzania/Tanzania/results") #UPDATE
datadir = file.path(getwd(),'rufiji_hydrodatainspect') #UPDATE
origdatadir = "F:/Tanzania/Tanzania/data"
outdir=file.path(getwd(),'rufiji_hydrodataimpute')
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}
rufidat <- read.csv(file.path(datadir,'rufidat_all.csv'), colClasses=c('character','Date','numeric','character','character'))
rufidat_clean <- read.csv(file.path(datadir,'rufidat_clean.csv'), colClasses=c('factor','Date','numeric','character','character'))
rufidat_deleted <- read.csv(file.path(datadir,'rufidat_deleted.csv'), colClasses=c('character','Date','numeric','character','character'))
gagesenv <- read.dbf(file.path(getwd(),'gages_netjoin.dbf'))

################################################
#Infilling/Interpolate data
################################################
#rufidat_clean$Flowlog <- log(rufidat_clean$Flow+0.01)
rufidat_cast <- as.data.frame(dcast(setDT(rufidat_clean), Date~ID, value.var='Flow'))

#Try auto.arima from forecast package and KalmanSmoother, inspired from  https://stats.stackexchange.com/questions/104565/how-to-use-auto-arima-to-impute-missing-values
#Function to find the index, for each row, of the previous row with a non-NA value
#Takes in a univariate time series with NAs
na.lomf <- function(x) {
  if (length(x) > 0) {
    non.na.idx <- which(!is.na(x))
    if(is.na(x[1]))
      non.na.idx=c(1,non.na.idx)
    rep.int(non.na.idx, diff(c(non.na.idx, length(x) + 1))) #Repeat index of previous row with non-NA as many times gap until next non-NA values
  }
}
#Function to find the index, for each row, of the next row with a non-NA value
#Takes in a univariate time series with NAs
na.lomb <- function(x) {
  if (length(x) > 0) {
    non.na.idx <- which(!is.na(x))
    if(is.na(x[length(x)]))
      non.na.idx=c(non.na.idx,length(x))
    y<-rep(NA,length(x))
    y[non.na.idx] <- non.na.idx
    zoo::na.locf(y, fromLast = TRUE)
  }
}

#################################Visualize correlation among gages########
rufidat_clean$Flowlog <- log(rufidat_clean$Flow+0.01)
setDT(rufidat_clean)[,ndays:=length(unique(Date)),.(ID)]
rufidat_castlog <- dcast(rufidat_clean[rufidat_clean$ndays>=10000,], Date~ID, value.var='Flowlog')
png(file.path(outdir,'hydropairs.png'),width=40, height=40,units='in',res=300)
hydropairs(as.data.frame(rufidat_castlog), dec=3, use="pairwise.complete.obs", method="pearson")
dev.off()
row.names(rufidat_castlog) <- rufidat_castlog$Date
rufidat_castlog <- rufidat_castlog[,-1]
library(Hmisc)
corrtab <- rcorr(as.matrix(rufidat_castlog), type="pearson")
corrtab_r <- corrtab$r
corrtab_p <- corrtab$P


########################## Interpolate data with na.interp ########################
#Use forecast na.interp for two gages that for which interpolation didn't work well
#Try out forecast package na.interp for seasonal series
CustomImpute_nainterp <- function(tscast, sn, maxgap,pplot=F) {
  if (sn > 1) {
    mindate <- tscast[min(which(!is.na(tscast[,sn]))),'Date']-37
    maxdate <- tscast[max(which(!is.na(tscast[,sn]))),'Date']+37
    tscastsub <- tscast[tscast$Date>mindate & tscast$Date<maxdate,]
    tscastsub$prevdate <- tscastsub[na.lomf(tscastsub[,sn]),'Date']
    tscastsub$nextdate <- tscastsub[na.lomb(tscastsub[,sn]),'Date']
    tscastsub$gap <- as.numeric(as.Date(tscastsub$nextdate)-as.Date(tscastsub$prevdate))
    tscastsub <- as.data.frame(tscastsub)
    id.na <- which((tscastsub$gap<=maxgap &
                      is.na(tscastsub[,sn])))
    bc <- BoxCox.lambda(tscastsub[,sn]+0.01,method='loglik',lower=0)
    pred_try <- data.frame(pred=na.interp(ts(tscastsub[,sn],frequency=365)), lambda=0)
    pred_try[pred_try$pred<min(tscastsub[,sn],na.rm=T),'pred'] <- min(tscastsub[,sn],na.rm=T) #if interpolated value is outside bound, bound it
    pred_try[pred_try$pred>max(tscastsub[,sn],na.rm=T),'pred'] <- max(tscastsub[,sn],na.rm=T)
    pred_try$Date <- tscastsub[,'Date']
    if (pplot==T){
      print(ggplot(tscastsub, aes(x=Date, y=get(colnames(tscastsub)[sn])))+
              geom_point(color='black')+
              geom_point(data=pred_try[id.na,],aes(y=pred), color='red') +
              ggtitle(colnames(tscastsub)[sn]) +
              scale_y_sqrt(name='Discharge (cms)')+
              scale_x_date(limits=c(as.Date(tscastsub[min(id.na),'Date']),
                                    as.Date(tscastsub[max(id.na),'Date'])))) #Only plot the area with NAs
    }
    out<-data.frame(pred_try$Date,pred_try$pred)
    colnames(out) <- c("Date",colnames(tscastsub)[sn])
    return(out)
  }
}

#Fill in every gauge
impute_preds <- data.frame(Date=rufidat_cast[,"Date"])
for (i in 2:(ncol(rufidat_cast))) {
  try({
    impute_preds<-merge(impute_preds, CustomImpute_nainterp(tscast=rufidat_cast,sn=i,maxgap=37, pplot=T), by='Date', all.x=T)
  })
}

#Notes on more data cleaning: 1KB31 (2001 peak, 2014 Jan-Feb), 1KB27 (1972-0973), 1KB24 (random peak 1973-74), 1KB18B (interp 1989, 2000-2001, 2004, low values few years before 2010),
# 1KB14A (interp low value 2003-2004), 1KA9 (interp 2014 low value, check value in 2012), 1KA8A (low value in 73-74 interp), 1KA42A (first records interp),
# 1KA31 (93-94 interp), 1KA15A (2000 interp), 1KA11A (1981 interp peak), 


write.csv(impute_preds, file.path(outdir, 'rufidat_interp.csv'), row.names=F)
impute_preds <- read.csv(file.path('rufiji_hydrodataimpute', 'rufidat_interp.csv'), colClasses=c('Date',rep('numeric',34)))
colnames(impute_preds)[2:(ncol(impute_preds))] <- substr(colnames(impute_preds),2,10)[2:(ncol(impute_preds))]

#############################################################################################################
#Plot clean and interpolated data
predsmelt <-melt(setDT(impute_preds),id.vars = 'Date',value.name='Flow',variable.name='ID')
predsmelt <- predsmelt[,c(2,1,3)]
predsmelt$SYM <- NA
predsmelt$Agency <- NA

for (gage in unique(predsmelt$ID)) {
  print(gage)
  genv <- gagesenv[gagesenv$RGS_No==gage,]
  gname <- paste(genv$RGS_Loc,"river at",genv$RGS_Name,sep=" ")
  #Generate FlowScreen time series
  gts<- create.ts(predsmelt[predsmelt$ID==gage,])  #Cannot run ts on multiple gages. Need to first subset by gage, then run ts.
  #Make raw time series plot
  rawsgplot <-ggplot(gts, aes(x=Date, y=Flow)) +
    #geom_rect(aes(xmin=as.Date('2001-10-01'), xmax=as.Date('2016-10-01'), ymin=min(gts$Flow-1), ymax=max(gts$Flow+1)), fill='#ffffbf', alpha=0.1) +
    geom_point(color='#bf812d', size=1.5)+ 
    geom_point(data=rufidat_deleted[rufidat_deleted$ID==gage,],aes(x=Date, y=Flow), color='#e31a1c', size=1.5) +
    geom_point(data=rufidat_clean[rufidat_clean$ID==gage,], color='#9ecae1', size=1.5) +
    geom_point(data=rufidat_clean[rufidat_clean$ID==gage & rufidat_clean$Date>'2001-10-01',], color='#045a8d', size=1.5) +
    scale_y_sqrt(expand=c(0,0))+
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") + 
    labs(y='Discharge (m3/s)', title=paste(gage, gname,sep=" - ")) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  #png(file.path(outdir,paste(gage,'raw_sg.png',sep="_")),width = 20, height=12,units='in',res=300)
  print(rawsgplot)
  #dev.off()
}

######################################## EXTRA ######################################################
# ggplot(tscastsub, aes(x=Date, y=get(colnames(tscastsub)[sn])))+
#   geom_point(data=try,aes(y=pred), color='red') + 
#   geom_point(color='black') 
########################## Try na.StructTS ########################
# tscast=as.data.frame(rufidat_cast)
# sn=2
# maxgap=180
# 
# #Restrict analysis to min and max year of records
# mindate <- tscast[min(which(!is.na(tscast[,sn]))),'Date']
# maxdate <- tscast[max(which(!is.na(tscast[,sn]))),'Date']+37
# tscastsub <- tscast[tscast$Date>mindate & tscast$Date<maxdate,]
# tscastsub <- ts(tscastsub[,sn],frequency=365)
# pred_try <- na.StructTS(tscastsub, maxgap=180) #Started at 5:25pm - ran for 24h for one station. Did not solve.
# 


############################################################
# #Try imputeTS package Kalman filter
# sn='1KA8A'
# mindate <- rufidat_castsub[min(which(!is.na(rufidat_castsub[,sn]))),'Date']
# x <- na.kalman(rufidat_castsub[rufidat_castsub$Date>mindate,sn], model = "auto.arima")
# rufidat_castsub[rufidat_castsub$Date>mindate,'pred1KA8A'] <- x
# ggplot(rufidat_castsub[rufidat_castsub$Date>mindate,], aes(x=Date, y=get(sn)))+
#   geom_point(aes(y=pred1KA8A), color='red') +
#   geom_point(color='black')+
#   scale_x_date(limits=as.Date(c('1965-01-01','1975-01-01')))
##########################################################
#Try out zoo na.approx - neither linear nor spline method is better
# tscastsubzoo <- zoo(tscastsub[,sn])
# coredata(tscastsubzoo)
# #pred_try <- data.frame(pred=na.approx(tscastsubzoo,x=index(tscastsubzoo), na.rm=F, rule=2,method='linear',maxgap=37))
# pred_try <- data.frame(pred=na.spline(tscastsubzoo,x=index(tscastsubzoo), na.rm=F,method='fmm'))
# pred_try$Date <- tscastsub[,'Date']
# ggplot(tscastsub, aes(x=Date, y=get(colnames(tscastsub)[sn])))+
#   geom_point(data=pred_try,aes(y=pred), color='red') +
#   geom_point(color='black')
#############################################################
#Function to fill gaps for one stream gauge based on an ARIMAX model and Kalman smoother using time series characteristics from that gauge and the other 
#gages as exogenous regressors.
#tscast data frame where the first is called 'Date' and contains dates and all the other columns are time flow series
#sn: column of the time series to be infilled
#maxgap: maximum gap length to be imputed
#Example values to test function:
#1KB8B, 1KB24, 1KB14A, 1KB50B, 1KA31, 1KA21A still have missing data
# tscast=rufidat_cast
# sn=34
# maxgap=180
# CustomImpute <- function(tscast, sn, maxgap) {
#   if (sn > 1) {
#     name<- colnames(tscast)[sn]
#     print(name)
#     #Restrict analysis to min and max year of records
#     mindate <- tscast[min(which(!is.na(tscast[,sn]))),'Date']-37
#     maxdate <- tscast[max(which(!is.na(tscast[,sn]))),'Date']+37
#     tscastsub <- tscast[tscast$Date>mindate & tscast$Date<maxdate,]
#     #Compute size of gap a record belongs to
#     tscastsub$prevdate <- tscastsub[na.lomf(tscastsub[,sn]),'Date']
#     tscastsub$nextdate <- tscastsub[na.lomb(tscastsub[,sn]),'Date']
#     tscastsub$prevgap <- as.numeric(tscastsub$Date-as.Date(tscastsub$prevdate))
#     tscastsub$nextgap <- as.numeric(as.Date(tscastsub$nextdate)-tscastsub$Date)
#     tscastsub$gap <- as.numeric(as.Date(tscastsub$nextdate)-as.Date(tscastsub$prevdate))
#     tscastsub <- as.data.frame(tscastsub)
#     #Only use data from stream gauges that have at least 50% of overlapping non-NA data with that streamgauge
#     overlap <- adply(tscastsub[!is.na(tscastsub[,sn]),-c(1,sn,ncol(tscastsub)-c(0,1,2,3,4))],2,function(x) {
#       length(which(!is.na(x)))/length(which(!is.na(tscastsub[,sn])))
#     })
#     covar <- which(colnames(tscastsub)%in%overlap[overlap$V1>=0.5,'X1'])
#     #Fit an ARIMAX model
#     pred <- tscastsub[,sn]
#     fit <- auto.arima(tscastsub[,sn],xreg=tscastsub[,covar]) 
#     #ARIMAX without forcing seasonality onto the model does not pan out, but ARIMA takes way too long to apply to all gauges
#     #Keep very simple interpolation for now, and increase complexity later
#     #fit <- auto.arima(ts(tscastsub[,sn],frequency=365),D=1, approximation=F)
#     #xreg=ts(tscastsub[,covar], frequency=365)
#     #summary(fit)
#     id.na <- which((tscastsub$gap<=maxgap | tscastsub$prevgap<=maxgap/2 | tscastsub$nextgap<=maxgap/2) &
#                      is.na(tscastsub[,sn]))
#     #Kalman smoother
#     kr <- KalmanSmooth(tscastsub[,sn], fit$model)
#     for (i in id.na)
#       pred[i] <- fit$model$Z %*% kr$smooth[i,]
#     #Output observed and predicted predicted data
#     tscastsub[id.na, sn] <- pred[id.na]
#     tscastsub[tscastsub[,sn]<0 & !is.na(tscastsub[,sn]),sn] <- 0
#     return(tscastsub[,c(1,sn)])
#   } else {
#     warning('sn must be >1 as 1st column must be "Date"')
#   }
# }
# 
# #Fill in every gauge
# impute_preds <- data.frame(Date=rufidat_cast[,"Date"])
# for (i in 2:(ncol(rufidat_cast))) {
#   try({
#     impute_preds<-merge(impute_preds, CustomImpute(rufidat_cast,i,maxgap=365), by='Date', all.x=T)
#   })
# }
# #Append the two gauges for which there was not enough info 1KA4A and 1KA33B
# impute_preds <- cbind(impute_preds,rufidat_cast[,which(!(colnames(rufidat_cast) %in% colnames(impute_preds)))])
