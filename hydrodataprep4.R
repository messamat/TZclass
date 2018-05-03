#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/28/2018
#Date last updated: 03/29/2018

#Purpose: impute/infill missing discharge data for all stream gauges

library(xlsx)
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
# rufidat_clean$Flowlog <- log(rufidat_clean$Flow+0.01)
# setDT(rufidat_clean)[,ndays:=length(unique(Date)),.(ID)]
# rufidat_castlog <- dcast(rufidat_clean[rufidat_clean$ndays>=10000,], Date~ID, value.var='Flowlog')
# png(file.path(outdir,'hydropairs.png'),width=40, height=40,units='in',res=300)
# hydropairs(as.data.frame(rufidat_castlog), dec=3, use="pairwise.complete.obs", method="pearson")
# dev.off()
# row.names(rufidat_castlog) <- rufidat_castlog$Date
# rufidat_castlog <- rufidat_castlog[,-1]
# library(Hmisc)
# corrtab <- rcorr(as.matrix(rufidat_castlog), type="pearson")
# corrtab_r <- corrtab$r
# corrtab_p <- corrtab$P


########################## Interpolate data with na.interp ########################
#Use forecast na.interp for two gages that for which interpolation didn't work well
#Try out forecast package na.interp for seasonal series
CustomImpute_nainterp <- function(tscast, sn, maxgap,pplot=F) {
  if (sn > 1) {
    print(colnames(tscast)[sn])
    #Only select records within 37 days of the first and last date of records for the gauge
    mindate <- tscast[min(which(!is.na(tscast[,sn]))),'Date']-37
    maxdate <- tscast[max(which(!is.na(tscast[,sn]))),'Date']+37
    tscastsub <- tscast[tscast$Date>mindate & tscast$Date<maxdate,]
    tscastsub$prevdate <- tscastsub[na.lomf(tscastsub[,sn]),'Date']
    tscastsub$nextdate <- tscastsub[na.lomb(tscastsub[,sn]),'Date']
    tscastsub$gap <- as.numeric(as.Date(tscastsub$nextdate)-as.Date(tscastsub$prevdate))
    tscastsub <- as.data.frame(tscastsub)
    
    #Determine Box-Cox lamba for tsl decomposition
    bc <- BoxCox.lambda(tscastsub[,sn]+10,method='loglik',lower=0)
    #Interpolate data
    pred_try <- data.frame(pred=na.interp(ts(tscastsub[,sn]+10,frequency=365), lambda=bc)-10)
    #Bound interpolation to min and max if outside bound of record
    pred_try[pred_try$pred<min(tscastsub[,sn],na.rm=T),'pred'] <- min(tscastsub[,sn],na.rm=T) #if interpolated value is outside bound, bound it
    pred_try[pred_try$pred>max(tscastsub[,sn],na.rm=T),'pred'] <- max(tscastsub[,sn],na.rm=T)
  
    #Create output dataframe
    out<-data.frame(tscastsub[,'Date'],pred_try$pred, tscastsub$gap)
    gapcol <- paste('gap',colnames(tscastsub)[sn],sep='_') 
    colnames(out) <- c("Date",colnames(tscastsub)[sn],gapcol)
    
    #Re-convert records for gaps longer than 'maxgap' to NA
    out[out[[gapcol]]>=maxgap, colnames(tscastsub)[sn]] <- NA
    
    #Flag interpolated vs. original data
    srccol <-  paste('source',colnames(tscastsub)[sn],sep='_') 
    out[out[[gapcol]]>0,srccol] <- 'interpolated' 
    out[out[[gapcol]]==0, srccol] <- 'original'
    
    #Plot original and interpolated data for subsequent QA/QC
    if (pplot==T){
      print(ggplot(out, aes(x=Date, y=get(colnames(tscastsub)[sn]), color=get(srccol))) +
              geom_point()+
              ggtitle(colnames(tscastsub)[sn]) +
              scale_y_sqrt(name='Discharge (cms)')+
              scale_x_date(limits=c(as.Date(min(rufidat_cast[!is.na(tscast[,sn]),'Date'])),
                                     as.Date(max(rufidat_cast[!is.na(tscast[,sn]),'Date']))))) #Only plot the area with NAs
    }
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

############QA/QC interpolated values/linearly interpolate some outliers
#Function to linearly interpolate for a set of dates.
interpcor <- function(df, station, cordates) {
  df[df$Date %in% cordates, station] <- NA
  df[,station] <- na.interpolation(df[,station], option='linear')
  gapcol <- paste('gap',station,sep='_') 
  df[df[[gapcol]]>37 | is.na(df[[gapcol]]),station] <- NA
  return(df)
}

ip1KB18B<- impute_preds[,c('Date','1KB18B','source_1KB18B')]
ggplot(ip1KB18B, aes(x=Date, y=`1KB18B`, color=source_1KB18B)) +
  geom_point()+
  scale_y_sqrt(name='Discharge (cms)')

#1KB18B: some abnormally low values
impute_preds <- interpcor(impute_preds, '1KB18B', as.Date(c('1987-05-17', '1987-06-11', '1989-05-27', '1989-06-24', '1990-05-12', '1990-07-01', '1999-11-08', '2004-11-06', 
                                            '2004-12-01', '2005-10-22', '2005-12-11', '2007-11-11', '2007-12-06', '2007-12-31', '2009-11-05',
                                            '2009-12-25', '2009-11-30','2009-10-11','2007-12-31','2007-12-06','2007-10-17','2006-12-21',
                                            '2006-11-26', '2006-11-01', '2006-10-07', '2006-09-12')))


#1KB17A: some abnormally high and low values
impute_preds <- interpcor(impute_preds, '1KB17A', as.Date(c('1972-08-13', '1973-06-02', '1973-06-27', '1974-06-15', '1976-06-06', '1977-06-16', '1977-08-12', '1978-06-01',
                                            '1978-06-26', '1979-05-30')))

#
bugval <- impute_preds[impute_preds$Date=='1995-09-06','1KA41']
impute_preds[impute_preds$`1KA41`<0.01 & !is.na(impute_preds$`1KA41`), '1KA41'] <- 0

write.csv(impute_preds, file.path(outdir, 'rufidat_interp.csv'), row.names=F)
impute_preds <- read.csv(file.path('rufiji_hydrodataimpute', 'rufidat_interp.csv'), colClasses=c('Date',rep(c('numeric','numeric','character'),39)))
colnames(impute_preds)[seq(2,ncol(impute_preds),3)] <- substr(colnames(impute_preds),2,10)[seq(2,ncol(impute_preds),3)]

################################Interpolate Kizigo data using precipitation data#####################
#Import and format data
dodoprecip <- openxlsx::read.xlsx(file.path(origdatadir,"RBWB_David20180430/Dodoma_Airport.xlsx"),sheet=1,startRow=1, detectDates=F)
dodoprecip$Date <- gsub(" ","", dodoprecip$Date)
dodoprecip$Date <- as.Date(as.character(dodoprecip$Date), format="%d/%m/%Y")
dodoprecip[dodoprecip$mm==-9.9,'mm'] <- NA
ggplot(dodoprecip[dodoprecip$Date<'1940-01-01',], aes(x=Date, y=mm)) + geom_line() + scale_y_sqrt()
#Need to change from 0.01 to 0

##############For 1KA41
precip1KA41 <- merge(impute_preds[,c('Date','1KA41')], dodoprecip, by='Date')
colnames(precip1KA41) <- c('Date','Flow','Precip')

#Try lm
# fit <- lm(Flow~Precip,data=precip1KA41)
# summary(fit)
# ggplot(precip1KA41, aes(x=Precip, y=Flow)) + geom_point() + geom_smooth()

#Create time series
mindate <- precip1KA41[min(which(!is.na(precip1KA41[,'Flow']))),'Date']
maxdate <- precip1KA41[max(which(!is.na(precip1KA41[,'Flow']))),'Date']
precip1KA41sub <- precip1KA41[precip1KA41$Date>mindate & precip1KA41$Date<maxdate,]
#Compute size of gap a record belongs to
precip1KA41sub$prevdate <- precip1KA41sub[na.lomf(precip1KA41sub[,'Flow']),'Date']
precip1KA41sub$nextdate <- precip1KA41sub[na.lomb(precip1KA41sub[,'Flow']),'Date']
precip1KA41sub$prevgap <- as.numeric(precip1KA41sub$Date-as.Date(precip1KA41sub$prevdate))
precip1KA41sub$nextgap <- as.numeric(as.Date(precip1KA41sub$nextdate)-precip1KA41sub$Date)
precip1KA41sub$gap <- as.numeric(as.Date(precip1KA41sub$nextdate)-as.Date(precip1KA41sub$prevdate))
precip1KA41sub <- as.data.frame(precip1KA41sub)
#Try with Fourier in response to https://robjhyndman.com/hyndsight/longseasonality/
#and https://robjhyndman.com/hyndsight/forecasting-weekly-data/
ts1KA41 <- msts(precip1KA41sub[, 'Flow'],seasonal.periods=365.25)
bestfit <- list(aicc=Inf)
for(i in 1:50)
{
  fit <- auto.arima(ts1KA41, seasonal=FALSE, xreg=cbind(fourier(ts1KA41, K=i), precip1KA41sub[, 'Precip']))
  if(fit$aicc < bestfit$aicc){
    print(i)
    bestfit <- fit
  }else {break};
}
kr <- KalmanSmooth(ts1KA41, bestfit$model)
id.na <- which(is.na(ts1KA41))
pred <- ts1KA41
for (i in id.na)
  pred[i] <- bestfit$model$Z %*% kr$smooth[i,]
precip1KA41sub[id.na,'Flow'] <- pred[id.na]
precip1KA41sub[precip1KA41sub$Flow<0,'Flow'] <- 0
ggplot(precip1KA41sub[,], aes(x=Date, y=Flow)) + geom_point() +
  geom_point(data=precip1KA41sub[id.na,], color='red')

##############For 1KA42



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
#Try auto.arima from forecast package and KalmanSmoother, inspired from  https://stats.stackexchange.com/questions/104565/how-to-use-auto-arima-to-impute-missing-values
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

###########################Try tbats: https://robjhyndman.com/publications/complex-seasonality/
#Gives error
# try <- tbats(ts1KA41)
# precip1KA41sub[1:7300, 'Flowinterp'] <- try
# ggplot(precip1KA41sub, aes(x=Date, y=Flowinterp)) + geom_point(color='red') +
#   geom_point(aes(y=Flow), color='black')
