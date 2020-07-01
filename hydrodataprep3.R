#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/22/2018
#Date last updated: 05/10/2018

#Purpose: visualize, assess quality, and perform preliminary clean-up of hydrological data in the Rufiji river basin, provided by CDMSmith Zach T. Eichenwald

library(ggplot2)
library(data.table)
library(FlowScreen) #to inspect data
library(waterData) #to import USGS data
library(prospectr) #for sg derivative
library(dplyr)
library(compare)
setwd("F:/Tanzania/Tanzania/results") #####UPDATE THIS TO MATCH YOUR PROJECT FOLDER #######
datadir = file.path(getwd(),'rufiji_hydrodataraw') #####UPDATE THIS TO MATCH YOUR PROJECT FOLDER #######
origdatadir = "F:/Tanzania/Tanzania/data" #####UPDATE THIS TO MATCH YOUR PROJECT FOLDER #######

#Define output directory
outdir=file.path(getwd(),'rufiji_hydrodatainspect')
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}

#Function to remove spurious constant values
remove_constant <- function(gage_data, gap_n=20) {
  gage_data[1,'Flag2'] <- 0
  delist <- list()
  for (i in 2:nrow(gage_data)){
    #Accumulate number of days with constant values
    if (gage_data[i,'Flow'] == gage_data[i-1,'Flow']){
      gage_data[i,'Flag2'] <- gage_data[i-1,'Flag2']+1 
    } else {
      gage_data[i,'Flag2'] <- 0
    }
    #Flag period with at least gap_n days of constant values
    if (gage_data[i,'Flag2'] == gap_n) {
      delist <- c(delist, gage_data[(i-19):i,'Date'])
    } 
    if (gage_data[i,'Flag2'] > gap_n) {
      delist <- c(delist, gage_data[i,'Date'])
    }
  }
  return(delist)
}

setClass('myDate')
#setAs("character","myDate", function(from)  as.POSIXlt(from, format= "%Y-%m-%d %H:%M:%S"))
#rufidat <- read.csv(file.path(datadir,'ZTE_rufidat.csv'),colClasses = c('character','character','myDate','numeric','numeric','factor','factor','factor'))
rufidat <- read.csv(file.path(datadir,'ZTE_rufidat.csv'),colClasses = c('character','character','Date','numeric','numeric','factor','factor','factor'))
RBWBflow <- read.csv(file.path(datadir, 'RBWBflowdat.csv'), colClasses = c('Date','numeric','character'))
JKdat <- read.csv(file.path(datadir, 'JK_dailydat.csv'), colClasses= c('Date', 'numeric', 'character'))
str(rufidat)

#General plotting of time series
rawplot <- ggplot(rufidat, aes(x=Date.Time, y=Calculated.Flow..cms.)) + 
  geom_point(color='#045a8d', size=1) + 
  geom_point(data=rufidat[rufidat$Calculated.Flow..cms.==0,],aes(x=Date.Time, y=Calculated.Flow..cms.), color='#e31a1c', size=1) +
  facet_wrap(~Gage.ID+Station.Name, scale='free') +
  theme_bw() + 
  labs(y='Discharge (m3/s)')
#rawplot
#png('rufidat_rawts.png',width = 32, height=16,units='in',res=300)
#rawplot
#dev.off()
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
#rufidat_screenform$Date.Time <- as.Date(rufidat_screenform$Date.Time)
rufidat_screenform <- rufidat_screenform[,list(Calculated.flow..cms.daily=mean(Calculated.Flow..cms., na.rm=T)), .(Gage.ID, Date.Time,Record.Quality,Rating.Curve.Source)]
rufidat_screenform <- rufidat_screenform[,c('Gage.ID','Date.Time','Calculated.flow..cms.daily','Record.Quality','Rating.Curve.Source')]
colnames(rufidat_screenform) <- colnames(test)
which(duplicated(rufidat_screenform[,c('Date','ID')]))

RBWBflow <- RBWBflow[,c(3,1,2)]
colnames(RBWBflow) <- c('ID','Date','Flow')
RBWBflow$SYM <- NA
RBWBflow$Agency <- 'RBWB_DavidMunkyala'

JKdat <- JKdat[,c(3,1,2)]
colnames(JKdat) <- c('ID','Date','Flow')
JKdat$SYM <- 'simulated'
JKdat$Agency <- 'JaphetKashaigili'

####RUN THE FOLLOWING COMMENTED-OUT BLOCK IN ORDER TO GET DIAGNOSTIC PLOTS USED IN DATA QA/QCing####
# for (gage in unique(rufidat_screenform$ID)) {
#   print(gage)
#   gname <- as.character(unique(rufidat[rufidat$Gage.ID==gage,'Station.Name']))
#   #Generate FlowScreen time series
#   gts<- create.ts(rufidat_screenform[rufidat_screenform$ID==gage,]) #Cannot run ts on multiple gages. Need to first subset by gage, then run ts.
#   #Compute and output flowScreen metrics and plots
#   try({
#     res <- metrics.all(gts)
#     ginfo <- data.frame(StationID=gage, StnName=gname, ProvState='Rufiji Basin',Country='Tanzania',Lat=0, Long=0, Area=0, RHN='RBWB')
#     png(file.path(outdir,paste(gage,'screenb.png',sep="_")),width = 20, height=12,units='in',res=300)
#     screen.summary(res, type="b", StnInfo=ginfo)
#     dev.off()
#     png(file.path(outdir,paste(gage,'screenl.png',sep="_")),width = 20, height=12,units='in',res=300)
#     screen.summary(res, type="l", StnInfo=ginfo)
#     dev.off()
#     png(file.path(outdir,paste(gage,'screenh.png',sep="_")),width = 20, height=12,units='in',res=300)
#     screen.summary(res, type="h", StnInfo=ginfo)
#     dev.off()
#   })
#   #Fit Savitzky-Golay 1st order derivative
#   # p = polynomial order w = window size (must be odd) m = m-th derivative (0 = smoothing) 
#   d1 <- as.data.frame(savitzkyGolay(gts$Flow, p = 3, w = 21, m = 1)) 
#   #A shorter period than 21 days would be used but it seems like granularity of stage measurements at low flow makes them appear very constant
#   gts[11:(nrow(gts)-10),'sg.d1'] <- d1
#   gts[gts$sg.d1<(10^-10) & gts$sg.d1>-(10^-10) & !is.na(gts$sg.d1),'Flag'] <- 'Y'
#   
#   #Make raw time series plot
#   rawsgplot <-ggplot(gts, aes(x=Date, y=Flow)) + 
#     geom_point(color='#045a8d', size=1) + 
#     geom_point(data=gts[gts$Flag=='Y',],aes(x=Date, y=Flow), color='#d01c8b', size=1.5) +
#     geom_point(data=gts[gts$Flow==0,],aes(x=Date, y=Flow), color='#e31a1c', size=1.5) +
#     scale_y_sqrt()+
#     scale_x_date(date_breaks = "2 years", date_labels = "%Y") + 
#     labs(y='Discharge (m3/s)', title=paste(gage, gname,sep=" - ")) +
#     theme_bw() + 
#     theme(axis.text.x = element_text(angle = 45, hjust=1))
#   png(file.path(outdir,paste(gage,'raw_sg.png',sep="_")),width = 20, height=12,units='in',res=300)
#   print(rawsgplot)
#   dev.off()
# }

#################################Clean out obviously spurious data in data from David Munkyala RBWB ###########################
RBWBflow_clean <- RBWBflow[!is.na(RBWBflow$Flow),]
###1KA9 KIMANI RIVER AT OLD GN: amomalous peak flood nealry at 2x10^7 cms
dm1KA9<- RBWBflow_clean[RBWBflow_clean$ID=='1KA9',]
ggplot(dm1KA9, aes(x=Date,y=Flow)) + geom_point() + scale_y_sqrt()
RBWBflow_clean <- RBWBflow_clean[!(RBWBflow_clean$ID=='1KA9' & RBWBflow_clean$Flow>10^6),]

###1KA21A: LITTLE RUAHA AT IHIMBU: a few erroneous values
dm1KA21A<- RBWBflow_clean[RBWBflow_clean$ID=='1KA21A',]
#ggplot(dm1KA21A[dm1KA21A$Date>'1957-01-01' & dm1KA21A$Date<'1960-01-01',], aes(x=Date,y=Flow)) + geom_point() + scale_x_date(date_breaks='1 year') + theme(axis.text.x=element_text(angle=90))
RBWBflow_clean <- RBWBflow_clean[!(RBWBflow_clean$ID=='1KA21A' & RBWBflow_clean$Date %in% as.Date(c('1968-05-29','1968-11-09', '1969-12-11','1959-07-31'))),]

###1KA31: LITTLE RUAHA AT MAWANDE:
dm1KA31<- RBWBflow_clean[RBWBflow_clean$ID=='1KA31',]
#ggplot(dm1KA31[dm1KA31$Date>'2011-01-01' & dm1KA31$Date<'2019-01-01',], aes(x=Date,y=Flow)) + geom_point() + scale_x_date(date_breaks='1 year') + theme(axis.text.x=element_text(angle=90))

###1KA39A LITTLE RUAHA @ IWAWA: suspiciously low values in December 1964 (jump back to 10 cms on January 1st)
dm1KA39A<- RBWBflow_clean[RBWBflow_clean$ID=='1KA39A',]
ggplot(dm1KA39A, aes(x=Date,y=Flow)) + geom_point()
RBWBflow_clean <- RBWBflow_clean[!(RBWBflow_clean$ID=='1KA39A' & RBWBflow_clean$Flow<0.5),]

###1KB31 Kigogo-Ruaha @ Lugema (Mshongo): anomalously low values during peak periods
dm1KB31 <- RBWBflow_clean[RBWBflow_clean$ID=='1KB31',]
setDT(dm1KB31)[, delta := shift(Flow, 1L, type="lead")-Flow,]
#ggplot(dm1KB31, aes(x=Date,y=Flow)) + geom_point() + scale_y_sqrt()
#ggplot(dm1KB31, aes(x=Date,y=delta)) + geom_point()

RBWBflow_clean <- RBWBflow_clean[!(RBWBflow_clean$ID=='1KB31' & 
                                     RBWBflow_clean$Date %in% as.Date(c('2001-03-11', '2001-04-30', '2002-02-15', '2002-03-09', '2002-03-22', '2002-03-23', '2002-03-25',
                                                                '2002-04-16', '2003-04-07', '2004-04-15', '2010-03-16', '2014-01-26','2014-04-09'))),]

###1KB32 Kihansi @ Lutaki: anomalous fall in 2013, 2015
dm1KB32 <- RBWBflow[RBWBflow$ID=='1KB32',]
setDT(dm1KB32)[, delta := shift(Flow, 1L, type="lead")-Flow,]
#ggplot(dm1KB32[dm1KB32$Date>'2007-01-01',], aes(x=Date,y=Flow)) + geom_line() + scale_y_sqrt()
RBWBflow_clean <- RBWBflow_clean[!(RBWBflow_clean$ID=='1KB32' & 
                                     RBWBflow_clean$Date %in% as.Date(c('2013-04-20', '2015-04-30'))),]

###1K3A Rufiji @ Stiegler's Gorge: unrealistic rises and falls (e.g. see March 1978) - cannot use data
dm1K3A <- RBWBflow[RBWBflow$ID=='1K3A',]
#ggplot(dm1K3A, aes(x=Date,y=Flow)) + geom_point() + scale_y_sqrt()

#################################Clean out obviously spurious data in data from Japhet Kashaigili UDSM ###########################
JKdat_clean <- JKdat[!is.na(JKdat$Flow),]
#1KB27 Luipa @ Mbingu: some discontinuities in record (sudden falls) + some suddenly low values near zero (1964)
jk1KB27 <- JKdat[JKdat$ID=='1KB27',]
ggplot(jk1KB27[jk1KB27$Date>'1970-01-01',], aes(x=Date,y=Flow)) + geom_point() + scale_y_sqrt()
JKdat_clean <- JKdat_clean[!(JKdat_clean$ID=='1KB27' & ((JKdat_clean$Date>='1963-12-01' & JKdat_clean$Date<='1964-01-31') | JKdat_clean$Flow<2)),]

#1KB28: Kihansi @ Ludoga
jk1KB28 <- JKdat[JKdat$ID=='1KB28',]
ggplot(jk1KB28, aes(x=Date,y=Flow)) + geom_point() + scale_y_sqrt()

#BBM2: Udajaji River @ Bridge (E 0816580, N 9047210, UTM Zone N 36)
#NOTE: In the current version of the classification, this station was discarded as having insufficient quality data due to a coding typo.
#If future classifications include simulated data, re-include the full dataset for this station. 
JKdat_clean[which(JKdat_clean$ID=='Udagaji' & !is.na(JKdat_clean$ID)),'ID'] <- 'BBM2'
BBM2 <- JKdat_clean[JKdat_clean$ID=='BBM2',]
ggplot(BBM2, aes(x=Date,y=Flow)) + geom_point()

#################################Clean out obviously spurious data in cleaned out data from ZTE###########################
rufidat_clean <- rufidat_screenform
###1KA2A	LITTLE RUAHA AT NDIUKA: period 1995 to late 2001 seems suspect: lots of missing data, appears like low-flow part of the year was cut-off or 
#                                                                         everything was shifted up?
g1KA2A<-rufidat_clean[rufidat_clean$ID=='1KA2A',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA2A' & rufidat_clean$Date>='1994-12-19' & rufidat_clean$Date<'2001-09-05'),] 


###1KA4A	GREAT RUAHA AT MSOSA: seems like stage measurements corresponding to a discharge right over 100cms are spurious as they plateau, 
#                               but not picked up by derivative. Will be kept (seem like the discharge at which there are spot measurements).

###1KA9	KIMANI RIVER AT OLD GN: overall variability seems unrealistic. Replace with RBWB David Munkyala data
g1KA9<-rufidat_clean[rufidat_clean$ID=='1KA9',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA9'),] 

###1KA11A seems like data prior to 165 are shifted up, could either be due to change in channel morphology, rating curve, or serious change inf low, remove
g1KA11A<-rufidat_clean[rufidat_clean$ID=='1KA11A',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA11A' & rufidat_clean$Date<'1964-11-15'),] 

###1KA15A	NDEMBERA AT ILONGO: data prior to 1968 are wonky. 
#ZTE posits that ' Old data are the outliers on the curve - possibility of old station (pre 1970) at different location?.  
g1KA15A<-rufidat_clean[rufidat_clean$ID=='1KA15A',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA15A' & rufidat_clean$Date<'1968-01-01'),]

###1KA21A	LITTLE RUAHA AT IHIMBU: data prior to 1967-11-23 area wonky, use David Munkyala's data instead 
g1KA21A<-rufidat_clean[rufidat_clean$ID=='1KA21A',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA21A'),]

###1KA22	MTITU AT MTITU: first few observations are spurious
g1KA22<-rufidat_clean[rufidat_clean$ID=='1KA22',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA22' & rufidat_clean$Date<'1957-05-27'),]

###1KA31: LITTLE RUAHA AT MAWANDE: Date pre-1970 are wonky, use David Munkyala's data instead
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA31'),]

###1KA32A	LITTLE RUAHA AT MAKALALA: remove all 0 values prior to 1970
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA32A' & rufidat_clean$Date<'1970-01-01' & rufidat_clean$Flow==0),]

###1KA37A	LUKOSI AT MTANDIKA: look wonky prior to 1968, also appears to be a flow duration curve break around 10 cms?
g1KA37A<-rufidat_clean[rufidat_clean$ID=='1KA37A',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA37A' & rufidat_clean$Date<'1968-01-01'),]

###1KA41	KIZIGO AT ILANGALI: seems like there was a shift in bed morphology post-2000. Overestimated 0-flows lead to long period of constant
#                             positive flow. Could also be some odd interpolation? Overestimated floods at the end of 1984
#                             During early period, seems like 0 flows are counted as 0.01. Replace with 0s.
g1KA41<-rufidat_clean[rufidat_clean$ID=='1KA41',]
#g1KA41_delist <- remove_constant(g1KA41)
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA41' & (rufidat_clean$Date>'2000-01-01' | rufidat_clean$Flow>500 |
                                                                (rufidat_clean$Date %in% g1KA41_delist & rufidat_clean$Flow > 1))),]
rufidat_clean[rufidat_clean$ID=='1KA41' & rufidat_clean$Flow<=0.01 & !is.na(rufidat_clean$Flow), 'Flow'] <- 0 
ggplot(g1KA41, aes(x=Date,y=Flow)) + geom_point() + scale_y_sqrt()

###1KA42A	KIZIGO RIVER AT CHINUGULU: Many instances of long zero flow period post 1993. Not artefacts (according to conversation with RBWB) 
g1KA42A<-rufidat_clean[rufidat_clean$ID=='1KA42A',]
ggplot(g1KA42A, aes(x=Date,y=Flow)) + geom_point() + scale_y_sqrt()
rufidat_clean[rufidat_clean$ID=='1KA42A' & rufidat_clean$Flow<=0.01 & !is.na(rufidat_clean$Flow), 'Flow'] <- 0 

###1KA50B	MSWISWI AT WILIMA: according to ZTE: 'Records poor pre 1968, post 1999. Possible station moved in 1999?' consider removing pre-1968

###1KA59	GREAT RUAHA AT MSEMBE FERRY: apparent shift pre- and post-1987-1990 gap due to water abstraction according to RBWB.

###1KA66	MLOWO AT ILONGO: anomalous period of constant values in late 2005-early 2006 + wonky 1993 data
# It in fact seems like the entire rating curve is off. 1KA51A which is right upstream has a discharge >20 times smaller even after cleaning
g1KA66<-rufidat_clean[rufidat_clean$ID=='1KA66',]
mean(g1KA66$Flow)
g1KA51A <-rufidat_clean[rufidat_clean$ID=='1KA51A',]
mean(g1KA51A$Flow)
#Gage has a watershed of 82km2, which means that if the mean annual discharge were to be 16cms, then it would require 6000mm of precipitation
#to fall on the watershed and 100% of it to accumulate at that location — very unlikely.
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KA66' & (rufidat_clean$Date<'1994-01-01' | 
                                                                (rufidat_clean$Date>'2005-07-23' & rufidat_clean$Date<'2006-02-04') |
                                                                rufidat_clean$Flow>1000)),] 
###1KA71	GREAT RUAHA AT NYALUHANGA: odd systematic break on descending limb in 2003, 2004, and 2005, but kept.

###1KB33	KIHANSI BELOW KIHANSI DAM: only contains -999: remove entirely
rufidat_clean <- rufidat_clean[rufidat_clean$ID!='1KB33',]

###1KB4A	KILOMBERO RIVER AT IFWEMA: long periods of spurious values
g1KB4A<-rufidat_clean[rufidat_clean$ID=='1KB4A',]
#Compute length of repetition and remove all records associated with periods of constant values of at least 20 days
g1KB4A_delist <- remove_constant(g1KB4A)
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KB4A' & rufidat_clean$Date %in% g1KB4A_delist),] 

###1KB8B	MPANGA RIVER AT MPANGA MISSION: two sets of seemingly low-flow spurious data in 1964 and 1976 
#                                       (also happen from one day to the next and encompass exactly one month, so remove)
g1KB8B<-rufidat_clean[rufidat_clean$ID=='1KB8B',]
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KB8B' & ((rufidat_clean$Date>'1963-12-01' & rufidat_clean$Date<'1964-01-30') | 
                                                                (rufidat_clean$Date>'1976-02-01' & rufidat_clean$Date<'1976-02-29') |
                                                                rufidat_clean$Flow<7)),] 
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
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KB15A' & (rufidat_clean$Flow < 5 | rufidat_clean$Date == '1970-04-13')),] 

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

###1KB24A SANJE AT SANJE: Very high peak at the beginning of the record. Way outside of the rating curve's range
#(max spot measurement at 6cms, here flood at nearly 2500 cms, corresponding to 18m stage increase — implausible)
#Some excessive peaks over 20cms during the short rains season pre-1985
g1KB24<- rufidat_clean[rufidat_clean$ID=='1KB24',]
ggplot(g1KB24[g1KB24$Date>'1972-01-01' & g1KB24$Date<'1975-01-01'], aes(x=Date,y=Flow)) + geom_point()

###1KB24A SANJE AT SANJE: Very high peak at the beginning of the record. Way outside of the rating curve's range
rufidat_clean <- rufidat_clean[!(rufidat_clean$ID=='1KB24' & rufidat_clean$Date<'1985-01-01'& rufidat_clean$Flow>20),] 

###1KB36 MGUGWE at MGUGWE
g1KB36<- rufidat_clean[rufidat_clean$ID=='1KB36',]
#ggplot(g1KB36, aes(x=Date,y=Flow)) + geom_point()

#################################Merge datasets#############################################
rufidat_all <- rbind(rufidat_screenform, 
                        RBWBflow[RBWBflow$ID %in% c('1KA9','1KA21A', '1KA31', '1KA39A','1KB31','1KB32'),], 
                        JKdat[JKdat$ID == '1KB27',]) #Does not include 1K3A, BBM2, or 1KB28

rufidat_cleanall <- rbind(rufidat_clean[rufidat_clean$ID!='1KA66',], 
                          RBWBflow_clean[RBWBflow_clean$ID %in% c('1KA9','1KA21A', '1KA31','1KA39A','1KB31','1KB32'),], 
                          JKdat_clean[JKdat_clean$ID == '1KB27',])

#################################Create dataset of deleted values###########################
rufidat_deleted <- anti_join(rufidat_screenform[!(rufidat_screenform$ID %in% c('1KA9','1KA21A','1KA31','1KA66')),], rufidat_clean, by=c("ID","Date"))
RBWBflow_deleted <- anti_join(RBWBflow, RBWBflow_clean, by=c("ID","Date"))
JKdat_deleted <- anti_join(JKdat, JKdat_clean, by=c("ID","Date"))
rufidat_deletedall <- rbind(rufidat_deleted, RBWBflow_deleted, JKdat_deleted)

#################################Export data to directory########################################
write.csv(rufidat_all, file.path(outdir, 'rufidat_all.csv'), row.names = F)
write.csv(rufidat_cleanall, file.path(outdir,'rufidat_clean.csv'),row.names = F)
write.csv(rufidat_deletedall, file.path(outdir,'rufidat_deleted.csv'),row.names = F)

################################# REPORT ##########################################
#Number of stations with some data (exlude 1KB33 as no actual data):
unique(rufidat_cleanall[!(rufidat_cleanall$ID %in% c('1KB33','1KB27','1KB28')),'ID'])
nrow(unique(rufidat_all[!(rufidat_all$ID %in% c('1KB33','1KB27','1KB28')),'ID'])) + 1 #For 1K3A 
#Of which 37 were kept for subsequent analysis
# 1KB27, 1KB28, and BBM2 that were not used in the classification due to insufficient data. 

#Total number of mean daily discharge records pre-cleaning
precleaning <- nrow(unique(as.data.frame(rufidat_all)[!(rufidat_all$ID %in% c('1KB33','1KB27','1KB28')),c('ID','Date')])) + nrow(dm1K3A)
precleaning

#Number of stations for which new rating curves were developed
unique(rufidat_all[,c('ID','Agency')])
nrow(unique(rufidat_cleanall[rufidat_cleanall$Agency=='New','ID'])) 
nrow(unique(rufidat_all[rufidat_all$Agency %in% c('RBWB_DavidMunkyala','RBWB',''),'ID']))

#Number of deleted records
postcleaning <- nrow(unique(as.data.frame(rufidat_cleanall)[!(rufidat_cleanall$ID %in% c('1KB33','1KB27','1KB28')),c('ID','Date')]))
(precleaning-postcleaning)/precleaning

