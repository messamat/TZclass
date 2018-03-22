library(plyr)
library(dplyr)
library(ggplot2)

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Date created: 03/15/2018
#Date last updated: 03/22/2018




setwd("F:/Tanzania/Tanzania/results")
#Import and merge ruaha data from Zach 
ruahaflow <- read.csv("F:/Tanzania/Tanzania/data/sharepoint20180316/flow/2018-03-06 Corrected Ruaha Stage and Flow_data.csv")
length(table(ruahaflow$Gage.ID))
ruahaflow <- ruahaflow[,!names(ruahaflow) %in% c("X","X.1")]
str(ruahaflow)
ruahastations <- read.csv("F:/Tanzania/Tanzania/data/sharepoint20180316/flow/2018-03-06 Corrected Ruaha Stage and Flow_stations.csv")
str(ruahastations)
ruahadat <- merge(ruahaflow, ruahastations[,c('Station_ID','Station_Name','Rating.Curve.Source','Record.Quality','Rating.Curve.Quality')], 
                  by.x='Gage.ID', by.y='Station_ID')
#Import and merge Kilombero data
kilomflow <- read.csv("F:/Tanzania/Tanzania/data/sharepoint20180316/flow/2018-03-13 Corrected Kilombero Stage and Flow_data.csv")
length(table(kilomflow$Gage.ID))



rufihydro <- read.csv("2018-01-23_HydroDatabase_stations.csv")

ruaha <- read.csv("F:/Tanzania/Tanzania/data/sharepoint20180311/flow/2018-03-06 Corrected Ruaha Stage and Flow_data.csv")
length(table(ruaha$Gage.ID))
str(ruaha)
ruaha$dateformat <- as.POSIXct(ruaha$Date.Time, format= "%m/%d/%Y %H:%M")
ruaha$date <- as.POSIXct(ruaha$Date.Time, format= "%m/%d/%Y")

summary<-ddply(ruaha, .(Gage.ID), summarize, mind=min(date), max=max(date), 
               inter=as.numeric((max(date)-min(date)))/365.25, rec=length(unique(date)),
               recper=length(unique(date))/as.numeric(max(date)-min(date)))
summary(summary)

ggplot(ruaha, aes(x=date, y=Calculated.Flow..cms., color=Gage.ID)) + geom_line() +
  scale_y_log10()


hydromerge <- merge(rufihydro, ruaha, by.x="Station_ID", by.y="GageID", all.y=F, all.x=T)
hydromerge <- hydromerge[!duplicated(hydromerge$Station_ID),]
str(hydromerge)
table(hydromerge$RATING) 
nrow(hydromerge[hydromerge$Q.MEASURE.COUNT >= 9,])


#Screening using hydroTSM (management, analysis, interpolation and plotting of flow time series)
#Assessment of daily streamflow time series quality with FlowScreen and/or EGRET packages
#Hydrologic metrics with hydrostats

#Prediction using caret




