#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/22/2018
#Date last updated: 03/24/2018

#Purpose: visualize, assess quality, and perform preliminary clean-up of hydrological data in the Rufiji river basin, provided by CDMSmith Zach T. Eichenwald

library(ggplot2)
library(data.table)
library(FlowScreen) #to inspect data
library(waterData) #to import USGS data
library(prospectr) #for sg derivative
library(compare)
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
g1KA2A<-rufidat_clean[rufidat_clean$ID=='1KA2A',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA2A' & rufidat_clean$Date>'1994-12-19' & rufidat_clean$Date<'2001-09-05'),] 

###1KA4A	GREAT RUAHA AT MSOSA: seems like stage measurements corresponding to a discharge right over 100cms are spurious as they plateau, 
#                               but not picked up by derivative. Will be kept (seem like the discharge at which there are spot measurements).

###1KA9	KIMANI RIVER AT OLD GN: very large peak in 1988 seems spurious as 02/22:131,02/23:>20,000,02/24:350 cms
g1KA9<-rufidat_clean[rufidat_clean$ID=='1KA9',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA9' & rufidat_clean$Date=='1987-02-23'),] 

###1KA15A	NDEMBERA AT ILONGO: data prior to 1968 are wonky. 
#ZTE posits that ' Old data are the outliers on the curve - possibility of old station (pre 1970) at different location?.  
g1KA15A<-rufidat_clean[rufidat_clean$ID=='1KA15A',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA15A' & rufidat_clean$Date<'1968-01-01'),]

###1KA21A	LITTLE RUAHA AT IHIMBU: data prior to 1967-11-23 area wonky 
g1KA21A<-rufidat_clean[rufidat_clean$ID=='1KA21A',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA21A' & rufidat_clean$Date<'1967-11-23'),]

###1KA22	MTITU AT MTITU: first few observations are spurious
g1KA22<-rufidat_clean[rufidat_clean$ID=='1KA22',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA22' & rufidat_clean$Date<'1957-05-27'),]

###1KA32A	LITTLE RUAHA AT MAKALALA: remove all 0 values prior to 1970
rufidat_deleted <- rbind(rufidat_deleted,rufidat_clean[(rufidat_clean$ID=='1KA32A' & rufidat_clean$Date<'1970-01-01' & rufidat_clean$Flow==0),]) 
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA32A' & rufidat_clean$Date<'1970-01-01' & rufidat_clean$Flow==0),]

###1KA37A	LUKOSI AT MTANDIKA: look wonky prior to 1968, also appears to be a flow duration curve break around 10 cms?
g1KA37A<-rufidat_clean[rufidat_clean$ID=='1KA37A',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA37A' & rufidat_clean$Date<'1968-01-01'),]

###1KA41	KIZIGO AT ILANGALI: seems like there was a shift in bed morphology post-2000. Overestimated 0-flows lead to long period of constant
#                             positive flow. Could also be some odd interpolation? 
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA41' & rufidat_clean$Date>'2000-01-01'),]

###1KA42A	KIZIGO RIVER AT CHINUGULU: Many instances of long zero flow period post 1993, could be due either to overwithdrawal or 
#                                     shift in bed morphology leading to spurious zero flows

###1KA50B	MSWISWI AT WILIMA: according to ZTE: 'Records poor pre 1968, post 1999. Possible station moved in 1999?' consider removing pre-1968

###1KA59	GREAT RUAHA AT MSEMBE FERRY: definite shift pre- and post-1987-1990 gap. Most likely different 
#ZTE: "No RC info for prior to 1994. Suspect poor records prior to 1990 when new rating curve (and perhaps station) were established. Likely decent records post 1990."

###1KA66	MLOWO AT ILONGO: anomalous period of constant values in late 2005-early 2006 + wonky 1993 data
g1KA66<-rufidat_clean[rufidat_clean$ID=='1KA66',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA66' & (rufidat_clean$Date<'1994-01-01' | 
                                                                (rufidat_clean$Date>'2005-07-23' & rufidat_clean$Date<'2006-02-04'))),] 

###1KA71	GREAT RUAHA AT NYALUHANGA: odd systematic break on descending limb in 2003, 2004, and 2005, but kept.

###1KB33	KIHANSI BELOW KIHANSI DAM: only contains -999: remove entirely
rufidat_clean <- rufidat_clean[rufidat_clean$ID!='1KB33',]

###1KB4A	KILOMBERO RIVER AT IFWEMA: long periods of spurious values
g1KB4A<-rufidat_clean[rufidat_clean$ID=='1KB4A',]
#Compute length of repetition and remove all records associated with periods of constant values of at least 20 days
remove_constant <- function(gage_data, gap_n=20) {
  gage_data[1,'Flag2'] <- 0
  delist <- list()
  for (i in 2:nrow(gage_data)){
    if (gage_data[i,'Flow'] == gage_data[i-1,'Flow']){
      gage_data[i,'Flag2'] <- gage_data[i-1,'Flag2']+1
    } else {
      gage_data[i,'Flag2'] <- 0
    }
    if (gage_data[i,'Flag2'] == gap_n) {
      delist <- c(delist, gage_data[(i-19):i,'Date'])
    } 
    if (gage_data[i,'Flag2'] > gap_n) {
      delist <- c(delist, gage_data[i,'Date'])
    }
  }
  return(delist)
}
g1KB4A_delist <- remove_constant(g1KB4A)
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KB4A' & rufidat_clean$Date %in% g1KB4A_delist),] 

###1KB8B	MPANGA RIVER AT MPANGA MISSION: two sets of seemingly low-flow spurious data in 1964 and 1976 
#                                       (also happen from one day to the next and encompass exactly one month, so remove)
g1KB8B<-rufidat_clean[rufidat_clean$ID=='1KB8B',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KB8B' & ((rufidat_clean$Date>'1963-12-01' & rufidat_clean$Date<'1964-01-30') | 
                                                                (rufidat_clean$Date>'1976-02-01' & rufidat_clean$Date<'1976-02-29'))),] 

###1KB9	MNYERA RIVER AT TAWETA: long periods of constant values
g1KB9<-rufidat_clean[rufidat_clean$ID=='1KB9',]
g1KB9_delist <- remove_constant(g1KB9)
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KB9' & rufidat_clean$Date %in% g1KB9_delist),] 

###1KB14A	LUMEMO AT KIBURUBUTU: quite a few spurious constant values post 1995. In addition pre mid-1967 data is off
g1KB14A<-rufidat_clean[rufidat_clean$ID=='1KB14A',]
g1KB14A_delist <- remove_constant(g1KB14A)
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KB14A' & (rufidat_clean$Date<'1966-07-01' | 
                                                                 (rufidat_clean$Date %in% g1KB14A_delist))),] 

###1KB15A	MNGETA RIVER AT MCHOMBE: a few spurious low flow values. Less than 20 records.
g1KB15A<-rufidat_clean[rufidat_clean$ID=='1KB15A',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KB15A' & rufidat_clean$Flow < 5),] 

###1KB16C	FURUA AT MALINYI: seemingly spurious 0 values
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KB16C' & rufidat_clean$Flow == 0),] 

###1KB17A	KILOMBERO RIVER AT SWERO: seemingly spurious 0 values
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KB17A' & rufidat_clean$Flow == 0),] 

###1KB18B	RUHUDJI BELOW KIFUNG'A FALLS: seemingly spurious 0 values and wondy data after April 2010
g1KB18B<-rufidat_clean[rufidat_clean$ID=='1KB18B',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KB18B' & (rufidat_clean$Flow == 0 | rufidat_clean$Date>'2010-06-05')),] 

###1KB19A	HAGAFIRO AT HAGAFIRO: spurious reccurring constant value
g1KB19A<-rufidat_clean[rufidat_clean$ID=='1KB19A',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KB19A' & rufidat_clean$Flow <=0.1812),] 

#################################
#Create dataset of deleted values
rufidat_deleted <- anti_join(rufidat_screenform, rufidat_clean, by=c("ID","Date"))

##########################################
#Assess general record characteristics
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