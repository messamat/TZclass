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
#setAs("character","myDate", function(from)  as.POSIXlt(from, format= "%Y-%m-%d %H:%M:%S"))
#rufidat <- read.csv(file.path(datadir,'ZTE_rufidat.csv'),colClasses = c('character','character','myDate','numeric','numeric','factor','factor','factor'))
rufidat <- read.csv(file.path(datadir,'ZTE_rufidat.csv'),colClasses = c('character','character','Date','numeric','numeric','factor','factor','factor'))
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
#rufidat_screenform$Date.Time <- as.Date(rufidat_screenform$Date.Time)
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
  d1 <- as.data.frame(savitzkyGolay(gts$Flow, p = 3, w = 21, m = 1)) 
  #A shorter period than 21 days would be used but it seems like granularity of stage measurements at low flow makes them appear very constant
  gts[11:(nrow(gts)-10),'sg.d1'] <- d1
  gts[gts$sg.d1<(10^-10) & gts$sg.d1>-(10^-10) & !is.na(gts$sg.d1),'Flag'] <- 'Y'
  
  #Make raw time series plot
  rawsgplot <-ggplot(gts, aes(x=Date, y=Flow)) + 
    geom_point(color='#045a8d', size=1) + 
    geom_point(data=gts[gts$Flag=='Y',],aes(x=Date, y=Flow), color='#d01c8b', size=1.5) +
    geom_point(data=gts[gts$Flow==0,],aes(x=Date, y=Flow), color='#e31a1c', size=1.5) +
    scale_y_sqrt()+
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") + 
    labs(y='Discharge (m3/s)', title=paste(gage, gname,sep=" - ")) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  png(file.path(outdir,paste(gage,'raw_sg.png',sep="_")),width = 20, height=12,units='in',res=300)
  print(rawsgplot)
  dev.off()
}

####################################################################################
# Clean out obviously spurious data based on preliminary observation
####################################################################################
rufidat_clean <- rufidat_screenform
###1KA2A	LITTLE RUAHA AT NDIUKA: period 1995 to late 2011 seems suspect: remove
#                                 Lots of missing data, appears like low-flow part of the year was cut-off or everything was shifted up?
#g1KA2A<-rufidat_clean[rufidat_clean$ID=='1KA2A',]
rufidat_deleted <- rufidat_clean[(rufidat_clean$ID=='1KA2A' & rufidat_clean$Date>'1994-12-19' & rufidat_clean$Date<'2001-09-05'),] 
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA2A' & rufidat_clean$Date>'1994-12-19' & rufidat_clean$Date<'2001-09-05'),] 
###1KA4A	GREAT RUAHA AT MSOSA: seems like stage measurements corresponding to a discharge right over 100cms are spurious as they plateau, 
#                               but not picked up by derivative. Will be kept.
###1KA9	KIMANI RIVER AT OLD GN: very large peak in 1988 seems spurious as 02/22:131,02/23:>20,000,02/24:350 cms
#g1KA9<-rufidat_clean[rufidat_clean$ID=='1KA9',]
rufidat_deleted <- rbind(rufidat_deleted,rufidat_clean[(rufidat_clean$ID=='1KA9' & rufidat_clean$Date=='1987-02-23'),]) 
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA9' & rufidat_clean$Date=='1987-02-23'),] 
###1KA15A	NDEMBERA AT ILONGO: data prior to 1970 are wonky. 
#ZTE posits that ' Old data are the outliers on the curve - possibility of old station (pre 1970) at different location?.  
g1KA9<-rufidat_clean[rufidat_clean$ID=='1KA9',]
rufidat_deleted <- rbind(rufidat_deleted,rufidat_clean[(rufidat_clean$ID=='1KA9' & rufidat_clean$Date=='1987-02-23'),]) 
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA9' & rufidat_clean$Date=='1987-02-23'),] 



###1KA21A	LITTLE RUAHA AT IHIMBU
###1KA22	MTITU AT MTITU
###1KA27	GREAT RUAHA AT MKUPULE

###1KA31	LITTLE RUAHA AT MAWENDE
###1KA32A	LITTLE RUAHA AT MAKALALA
###1KA33B	NDEMBERA AT MADIBILA
###1KA37A	LUKOSI AT MTANDIKA
###1KA38A	YOVI AT YOVI
###1KA41	KIZIGO AT ILANGALI
###1KA42A	KIZIGO RIVER AT CHINUGULU
##
###1KA50B	MSWISWI AT WILIMA
###1KA51A	UMROBO AT GNR
###1KA57A	MWEGA AT MALOLO
###1KA59	GREAT RUAHA AT MSEMBE FERRY
###1KA66	MLOWO AT ILONGO
###1KA71	GREAT RUAHA AT NYALUHANGA



###1KB14A	LUMEMO AT KIBURUBUTU:
###1KB15A	MNGETA RIVER AT MCHOMBE
###1KB16C	FURUA AT MALINYI
###1KB17A	KILOMBERO RIVER AT SWERO
###1KB18B	RUHUDJI BELOW KIFUNG'A FALLS
###1KB19A	HAGAFIRO AT HAGAFIRO
###1KB24	SANJE RIVER AT SANJE
###1KB33	KIHANSI BELOW KIHANSI DAM: only contains -999: remove
rufidat <- rufidat[rufidat$Gage.ID!='1KB33',]

###1KB36	MGUGWE RIVER AT MGUGWE
###1KB36	MGUGWE RIVER AT MGUGWE
###1KB4A	KILOMBERO RIVER AT IFWEMA
###1KB8B	MPANGA RIVER AT MPANGA MISSION
###1KB9	MNYERA RIVER AT TAWETA
R






#Notes on gages


##########################################
#Assess general record characteristics
rufidat$year <- format(rufidat$Date.Time, "%Y")
rufidat$month <- as.numeric(format(rufidat$Date.Time, "%m"))

#Test range of maximum gap length per year and total percentage missing data
#



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