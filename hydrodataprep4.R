#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/26/2018
#Date last updated: 03/26/2018

#Purpose: filter our streamgages based on length of record, gaps, overlap, and non-stationarity

library(ggplot2)
library(data.table)
library(FlowScreen) #to inspect data
library(compare)
setwd("F:/Tanzania/Tanzania/results") #UPDATE
datadir = file.path(getwd(),paste('rufiji_hydrodatainspect','20180326',sep='_')) #UPDATE
origdatadir = "F:/Tanzania/Tanzania/data"


rufidat_clean <- read.csv(file.path(datadir,'rufidat_clean.csv'), colClasses=c('factor','Date','numeric','character','character'))
rufidat_deleted <- read.csv(file.path(datadir,'rufidat_deleted.csv'), colClasses=c('character','Date','numeric','character','character'))

##########################################
#BUild summary tables
rufidat_clean$year <- as.numeric(format(rufidat_clean$Date, "%Y"))
rufidat_clean$month <- as.numeric(format(rufidat_clean$Date, "%m"))
rufidat_clean$prevdate <- NA #Get previous date of record for that station for each date
for (i in 2:nrow(rufidat_clean)){
  if (rufidat_clean[i,'ID']==rufidat_clean[i-1,'ID']) {
    rufidat_clean[i,'prevdate'] < rufidat_clean[i-1,'Date']
  }
}

rufidat_dt <- data.table(rufidat_clean)
rufidat_datesummary <- rufidat_dt[,list(length(unique(ID))), .(Date)]

ggplot(data=rufidat_clean, aes(x=Date, y=ID)) + 
  geom_bar(data=rufidat_datesummary, aes(x=Date,y='1KB9',color=V1), stat='identity') +
  geom_point() +
  scale_x_date(date_breaks ="5 years") +
  scale_colour_distiller('spectral') +
  theme_bw()

rufidat_summary <- rufidat_dt[,list(min_year=min(year), max_year=max(year), max_len=max(year)-min(year), 
                                    max_lend=max(Date)-min(Date), gap_d=as.numeric((max(Date)-min(Date)))-length(unique(Date)), gap_per=1-(length(unique(Date))/as.numeric((max(Date)-min(Date)))),
                                    max_gap = max(Date-prevdate),
                                       ), .(ID)]

rufidat_gapsummary <- rufidat_dt[,list(gap_d=as.numeric(format(as.Date(paste(year, "12", "31", sep="-")), "%j"))-len(unique(Date)),
                                          gap_per=length(unique(Date))/as.numeric(format(as.Date(paste(year, "12", "31", sep="-")), "%j")),
                                          ), .(ID,year)]


#Test range of maximum gap length per year and total percentage missing data

#####

#To do:
#Get summary statistics on 
# grain of data
# length of record
# completeness in terms of frequency, length, and time periods of gaps
# overlap in terms of period and length
#Test range of acceptance criteria (15 years, etc.)
#Evaluate consistency of discharge with drainage area and precipitation + HydroSHEDS modeled data
