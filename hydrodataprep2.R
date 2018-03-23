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
library(data.table)
setwd("F:/Tanzania/Tanzania/results") #UPDATE
datadir = file.path(getwd(),paste('rufiji_hydrodataraw','20180322',sep='_')) #UPDATE

setClass('myDate')
setAs("character","myDate", function(from)  as.Date(from, format= "%m/%d/%Y"))
rufidat <- read.csv(file.path(datadir,'ZTE_rufidat.csv'),colClasses = c('factor','factor','myDate','numeric','numeric','factor','factor','factor'))
str(rufidat)
JK_dailydat <- read.csv(file.path(datadir,'JK_dailydat.csv'),colClasses = c('Date','numeric','factor'))
JK_Unimpaired <- read.csv(file.path(datadir,'JK_monthly_unimpaired.csv'),colClasses = c('Date','factor','numeric'))
JK_Rukwa <- read.csv(file.path(datadir,'JK_monthly_rukwa.csv'),colClasses = c('Date','factor','numeric'))

#Check which stations are both in JK daily dataset and ZKE dataset and merge them
unique(JK_dailydat$Station)[(which(unique(JK_dailydat$Station) %in% unique(rufidat$Gage.ID)))]
rufidatsub <- data.table(rufidat[rufidat$Gage.ID %in% unique(JK_dailydat$Station),])
JK_dailydatsub <- JK_dailydat[JK_dailydat$Station %in% unique(rufidatsub$Gage.ID),]
rufidatsub <- rufidatsub[,list(Station.Name, Calculated.flow..cms.daily=mean(Calculated.Flow..cms.)), .(Gage.ID, Date.Time)] #Get daily average discharge from ZTE data for comparison with JK
JK_ZTE_dat <- merge(rufidatsub,JK_dailydatsub, by.x=c('Gage.ID','Date.Time'), by.y=c('Station','Date'), all.x=T, all.y=T)
colnames(JK_ZTE_dat) <- c("Gage.ID", "Date","Station.Name", "Flow_ZTE", "Flow_JK")

#Compare data graphically
comparison1 <- ggplot(JK_ZTE_dat, aes(x=Flow_ZTE, y=Flow_JK, color=Gage.ID)) + 
  geom_point(alpha=0.5) + 
  facet_wrap(~Gage.ID, scale='free') + 
  geom_abline(slope=1, intercept=0, size=1) + 
  theme_bw()
pdf('ZTE_JK_datacomparison1.pdf',width = 8, height=8)
comparison1
dev.off()

comparison2 <- ggplot(JK_ZTE_dat, aes(x=Date, y=Flow_JK)) + 
  geom_line(color='lightblue') + 
  geom_line(aes(y=Flow_ZTE), color='red', alpha=0.5) +
  facet_wrap(~Gage.ID, scale='free') +
  scale_y_sqrt() +
  theme_bw() + 
  labs(y='Discharge (m3/s) - Blue: JK, Red:ZTE')
comparison2
pdf('ZTE_JK_datacomparison2.pdf',width = 11.5, height=8.5)
comparison2
dev.off()






