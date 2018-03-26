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
library(FlowScreen) #to inspect data
library(waterData) #to import USGS data
library(prospectr) #for sg derivative
setwd("F:/Tanzania/Tanzania/results") #UPDATE
datadir = file.path(getwd(),paste('rufiji_hydrodataraw','20180324',sep='_')) #UPDATE
origdatadir = "F:/Tanzania/Tanzania/data"

setClass('myDate')
setAs("character","myDate", function(from)  as.POSIXlt(from, format= "%Y-%m-%d %H:%M:%S"))
rufidat <- read.csv(file.path(datadir,'ZTE_rufidat.csv'),colClasses = c('character','character','myDate','numeric','numeric','factor','factor','factor'))
str(rufidat)

#General plotting of time series
rawplot <- ggplot(rufidat, aes(x=Date.Time, y=Calculated.Flow..cms.)) + 
  geom_point(color='#045a8d', size=1) + 
  geom_point(data=rufidat[rufidat$Calculated.Flow..cms.==0,],aes(x=Date.Time, y=Calculated.Flow..cms.), color='#e31a1c', size=1) +
  facet_wrap(~Gage.ID+Station.Name, scale='free') +
  theme_bw() + 
  labs(y='Discharge (m3/s)')
rawplot
#png('rufidat_rawts.png',width = 32, height=16,units='in',res=300)
#rawplot
#dev.off()

#Remove KB33 (below Kihansi, as only contains -999)
rufidat <- data.table(rufidat[rufidat$Gage.ID!='1KB33',])
nrow(rufidat[rufidat$Calculated.Flow..cms.==-999,])

########################################################################################################
#Format rufiji flow data to be used within the FlowScreen package following the USGS data import format
#FlowScreen package was developed to work with Water Survey of Canada (WSC) or the United States Geological Survey (USGS) data. 

#Download and import data from USGS using FlowScreen to see output of 'read.flows function'
# testdat <- importDVs('06135000', code = "00060", stat = "00003", sdate = "1851-01-01",
#                      edate = as.Date(Sys.Date(), format = "%Y-%m-%d"))
# testtab=file.path(getwd(),paste('test_',as.character(format(Sys.Date(),'%Y%m%d')),'.csv',sep=""))
# write.csv(testdat, testtab)
# test<-read.flows(testtab)
# colnames(test)

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
  gname <- as.character(unique(rufidat[rufidat$Gage.ID==gage,'Station.Name']))
  #Generate FlowScreen time series
  gts<- create.ts(rufidat_screenform[rufidat_screenform$ID==gage,]) #Cannot run ts on multiple gages. Need to first subset by gage, then run ts.
  #Compute and output flowScreen metrics and plots
  try({
    res <- metrics.all(gts)
    ginfo <- data.frame(StationID=gage, StnName=gname, ProvState='Rufiji Basin',Country='Tanzania',Lat=0, Long=0, Area=0, RHN='RBWB')
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
  #Fit Savitzky-Golay 1st order derivative
  # p = polynomial order w = window size (must be odd) m = m-th derivative (0 = smoothing) 
  d1 <- as.data.frame(savitzkyGolay(gts$Flow, p = 3, w = 11, m = 1))
  gts[6:(nrow(gts)-5),'sg.d1'] <- d1
  gts[gts$sg.d1<(10^-10) & gts$sg.d1>-(10^-10) & !is.na(gts$sg.d1),'Flag'] <- 'Y'
  
  #Make raw time series plot
  rawsgplot <-ggplot(gts, aes(x=Date, y=Flow+0.01)) + 
    geom_point(color='#045a8d', size=1) + 
    geom_point(data=gts[gts$Flag=='Y',],aes(x=Date, y=Flow), color='#d01c8b', size=1.5) +
    geom_point(data=gts[gts$Flow==0,],aes(x=Date, y=Flow), color='#e31a1c', size=1.5) +
    theme_bw() + 
    scale_y_log10(limits=c(0,max(gts$Flow)))+
    labs(y='Discharge (m3/s)', title=paste(gage, gname,sep=" - "))
  png(file.path(outdir,paste(gage,'raw_sg.png',sep="_")),width = 20, height=12,units='in',res=300)
  rawsgplot
  
  dev.off()
}

#########################################################
# Clean out spurious data based on preliminary observation
#########################################################





#Notes on gages


##########################################
#Assess general record characteristics
rufidat$year <- format(rufidat$Date.Time, "%Y")
rufidat$month <- as.numeric(format(rufidat$Date.Time, "%m"))





res <- metrics.all(ts, Qmax = 0.95, Dur = 5, Qdr = 0.2, WinSize = 30,season = c(4:9), NAthresh = 0.5, language = "English")

#####


#To do:
#Get summary statistics on 
# grain of data
# length of record
# completeness in terms of frequency, length, and time periods of gaps
# overlap in terms of period and length
#Test range of acceptance criteria (15 years, etc.)
#Evaluate consistency of discharge with drainage area and precipitation + HydroSHEDS modeled data
#Add StnInfo to screen.summary