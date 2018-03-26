#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/22/2018
#Date last updated: 03/23/2018

#Purpose: compare data provided by CDMSmith Zachary T. Eichenwald (ZTE) and Japhet Kashaigili (JK)

library(ggplot2)
library(data.table)
setwd("F:/Tanzania/Tanzania/results") #UPDATE
datadir = file.path(getwd(),paste('rufiji_hydrodataraw','20180323',sep='_')) #UPDATE
origdatadir = "F:/Tanzania/Tanzania/data"

#setClass('myDate')
#setAs("character","myDate", function(from)  as.POSIXct(from, format= "%m/%d/%Y %H:%M"))
rufidat <- read.csv(file.path(datadir,'ZTE_rufidat.csv'),colClasses = c('character','factor','Date','numeric','numeric','factor','factor','factor'))
str(rufidat)
JK_dailydat <- read.csv(file.path(datadir,'JK_dailydat.csv'),colClasses = c('Date','numeric','character'))
JK_dailydat <- JK_dailydat[!is.na(JK_dailydat$Station),]
JK_Unimpaired <- read.csv(file.path(datadir,'JK_monthly_unimpaired.csv'),colClasses = c('Date','factor','numeric'))
JK_Rukwa <- read.csv(file.path(datadir,'JK_monthly_rukwa.csv'),colClasses = c('Date','factor','numeric'))
rufidat[rufidat$Gage.ID=='1KA11A' & !is.na(rufidat$Gage.ID),'Gage.ID'] <- '1KA11'
rufidat[rufidat$Gage.ID=='1KB8B' & !is.na(rufidat$Gage.ID),'Gage.ID'] <- '1KB8'
rufidat[rufidat$Gage.ID=='1KB17A' & !is.na(rufidat$Gage.ID),'Gage.ID'] <- '1KB17'

#Check which stations are both in JK daily dataset and ZKE dataset and merge them
unique(JK_dailydat$Station)[(which(unique(JK_dailydat$Station) %in% unique(rufidat$Gage.ID)))]
rufidatsubraw <- data.table(rufidat[rufidat$Gage.ID %in% unique(JK_dailydat$Station),])
JK_dailydatsub <- JK_dailydat[JK_dailydat$Station %in% unique(rufidatsubraw$Gage.ID),]
rufidatsub <- rufidatsubraw[,list(Calculated.flow..cms.daily=mean(Calculated.Flow..cms.)), .(Gage.ID, Date.Time)] #Get daily average discharge from ZTE data for comparison with JK
JK_ZTE_dat <- merge(rufidatsub,JK_dailydatsub, by.x=c('Gage.ID','Date.Time'), by.y=c('Station','Date'), all.x=T, all.y=T)
JK_ZTE_dat <- merge(JK_ZTE_dat, unique(rufidatsubraw[,c("Gage.ID",'Station.Name')]),by='Gage.ID',all.x=T)
colnames(JK_ZTE_dat) <- c("Gage.ID", "Date", "Flow_ZTE", "Flow_JK","Station.Name")

#Compare data graphically
comparison1 <- ggplot(JK_ZTE_dat, aes(x=Flow_ZTE, y=Flow_JK, color=Gage.ID)) + 
  geom_point(alpha=0.5) + 
  facet_wrap(~Gage.ID+Station.Name, scale='free') + 
  geom_abline(slope=1, intercept=0, size=1) + 
  labs(x='Daily mean discharge from Zach E. (cms)', y='Daily mean discharge from Japhet K. (cms)')+
  theme_bw()
comparison1
pdf('ZTE_JK_datacomparison1.pdf',width = 11.5, height=8)
comparison1
dev.off()

comparison2 <- ggplot(JK_ZTE_dat, aes(x=Date, y=Flow_JK)) + 
  geom_line(color='blue') + 
  geom_line(aes(y=Flow_ZTE), color='red', alpha=0.75) +
  facet_wrap(~Gage.ID+Station.Name, scale='free') +
  scale_y_sqrt() +
  theme_bw() + 
  labs(y='Discharge (m3/s) - Blue: JK, Red:ZTE')
comparison2
pdf('ZTE_JK_datacomparison2.pdf',width = 11.5, height=8.5)
comparison2
dev.off()

##################################################
#Inspect data from 1KA32A Little Ruaha at Makalala
##################################################
gage1KA32A <- JK_ZTE_dat[JK_ZTE_dat$Gage.ID =='1KA32A',] 
gage1KA32A[gage1KA32A$Flow_ZTE>0,][order(Date),'Date'][1] #No data over 0 cms prior to February 1st 1965
#Check whether this issue is present in original data from ZTE
ruahaflow <- read.csv(file.path(origdatadir,"sharepoint20180316/flow/2018-03-06 Corrected Ruaha Stage and Flow_data.csv"),,colClasses = c('factor','factor','myDate','numeric','numeric')) 
ruahaflow <- ruahaflow[,!names(ruahaflow) %in% c("X","X.1")]
orig1KA32A <- ruahaflow[ruahaflow$Gage.ID =='1KA32A',] #Definitely in the original data from Zach — needs to be deleted

##################################################
#Inspect data from 1KB9 Mynera river at Taweta
##################################################
ZE1KB9 <- rufidat[rufidat$Gage.ID =='1KB9',]#It appears like there is a unit mismatch between JK and 
gage1KB9 <- JK_ZTE_dat[JK_ZTE_dat$Gage.ID =='1KB9',] 
gage1KB9$ratio <- gage1KB9$Flow_JK/gage1KB9$Flow_ZTE
gage1KB9$Flow_ZTE_cor <- gage1KB9$Flow_ZTE*mean(gage1KB9$ratio,na.rm=T)
comparison1KB9 <- ggplot(gage1KB9, aes(x=Flow_ZTE_cor, y=Flow_JK)) + 
  geom_point(alpha=0.5) + 
  geom_abline(slope=1, intercept=0, size=1) + 
  labs(x='Daily mean discharge from Zach E. (cms)', y='Daily mean discharge from Japhet K. (cms)')+
  theme_bw()
comparison1KB9

comparison1KB9_2 <- ggplot(gage1KB9, aes(x=Date, y=Flow_JK)) + 
  geom_line(color='blue') + 
  geom_line(aes(y=Flow_ZTE_cor), color='red', alpha=0.75) +
  scale_y_sqrt() +
  theme_bw() + 
  labs(y='Discharge (m3/s) - Blue: JK, Red:ZTE')
comparison1KB9_2

#It definitely seems like a unit difference + obviously different rating curves. But there also seems to be constant filled in data in ZTE data?
