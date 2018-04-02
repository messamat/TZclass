#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/30/2018
#Date last updated: 04/02/2018

#Purpose: compute hydrological metrics and classify stream gauges

library(plyr)
library(dplyr)
library(hydrostats)
library(data.table)
devtools::install_github("USGS-R/EflowStats") #Intro Vignette: https://cdn.rawgit.com/USGS-R/EflowStats/9507f714/inst/doc/intro.html
library(EflowStats)

rootdir="F:/Tanzania/Tanzania" #UPDATE
setwd(file.path(rootdir,"results")) 
datadir = file.path(getwd(),paste('rufiji_hydrodatafilter','20180327',sep='_')) #UPDATE
origdatadir = file.path(rootdir,"data") 
outdir=file.path(getwd(),'rufiji_hydrodatastats_20180331')
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}

gagesenv <- read.dbf(file.path(getwd(),'gages_netjoinclean.dbf'))
gagesenvrec <- merge(gagesenv, unique(rufidat_clean[,c('ID','SYM')]), by.x='RGS_No', by.y='ID', all.x=F)

rufidat_clean <- read.csv(file.path('rufiji_hydrodatainspect_20180326','rufidat_clean.csv'), colClasses=c('factor','Date','numeric','character','character'))
rufidat_impute <- read.csv(file.path('rufiji_hydrodataimpute_20180329', 'rufidat_interp.csv'), colClasses=c('Date',rep('numeric',34)))
colnames(rufidat_impute)[2:(ncol(rufidat_impute))] <- substr(colnames(rufidat_impute),2,10)[2:(ncol(rufidat_impute))]
rufidat_gapsummary <- read.csv(file.path(datadir, 'rufidat_gapsummary.csv'))
rufidat_post1991<-read.csv(file.path(datadir, 'gageselect_post1991comp90.csv'))

predsmelt <-melt(setDT(rufidat_impute),id.vars = 'Date',value.name='Flow',variable.name='ID')
predsmelt <- predsmelt[,c(2,1,3)]
predsmelt$year <- as.numeric(format(predsmelt$Date, "%Y"))
predsmelt$month <- as.numeric(format(predsmelt$Date, "%m"))
predsmelt<-hyear.internal(predsmelt,hyrstart=10) #Ignore hdoy
predsmelt <- merge(predsmelt, rufidat_gapsummary, by=c('ID','hyear'),all.x=T)
predsmelt <- merge(predsmelt, rufidat_post1991, by='ID',all.x=T)

rufidat_select <- predsmelt[predsmelt$hyear>=1991 & predsmelt$hyear<2017 & predsmelt$ycount>=10 & predsmelt$gap_per<=0.1,]
str(rufidat_select)

#####################################################################################
# Compute hydrologic metrics
HITcomp <- function(dfhydro, dfenv, gageID, hstats="all", floodquantile=0.95) {
  #Check data for completeness
  dailyQClean <- validate_data(dfhydro[dfhydro$ID==gageID,c("Date", "Flow")], yearType="water")
  #Calculate all hit stats
  calc_allHITout <- calc_allHIT(dailyQClean, yearType="water", stats=hstats, digits=3, pref="mean",
                                drainArea=dfenv[dfenv$RGS_No==gageID,'WsArea'], floodThreshold = quantile(dailyQClean$discharge, floodquantile))
  return(calc_allHITout)
}
#Get template
HITall_template <- HITcomp(rufidat_select, gagesenvrec, '1KA9')
colnames(HITall_template)[2] <- '1KA9'
HITall <- data.frame(indice=HITall_template$indice) 
#Compute metrics for all gages
for (gage in unique(rufidat_select$ID)) {
  print(gage)
  try({
    calc_allHITout <- HITcomp(rufidat_select, gagesenvrec, gage)
    colnames(calc_allHITout)[2] <- gage
    HITall <- merge(HITall, calc_allHITout, by='indice')
  })
}

gageID <- '1KA51A'
dailyQClean <- validate_data(rufidat_select[rufidat_select$ID==gageID,c("Date", "Flow")], yearType="water")
#Calculate all hit stats
calc_allHITout <- calc_allHIT(dailyQClean, yearType="water", stats='all', digits=3, pref="mean",
                              drainArea=dfenv[dfenv$RGS_No==gageID,'WsArea'], floodThreshold = quantile(dailyQClean$discharge, 0.95))

colnames(calc_allHITout)[2] <- gage
HITall <- merge(HITall, calc_allHITout, by='indice')

#Troubleshoot error message (check out https://github.com/USGS-R/EflowStats/blob/23a680d8a4c657674828bc5d32f2d6a887f14fbb/R/calc_timingAverage.R):
#Error in findInterval(x$log_discharge, break_pts) : 'vec' must be sorted non-decreasingly and not contain NAs"""
gageID <- '1KA51A'
dailyQClean <- validate_data(rufidat_select[rufidat_select$ID==gageID,c("Date", "Flow")], yearType="water")
#Calculate all hit stats



#Calculate mag7 stats
#magnifStatsOut <- calc_magnifSeven(dailyQClean,yearType="water",digits=3)


######################################################################
#Classify gauges

#Standardize metrics
 #First magnitude metrics by drainage area


hydrostats_tra <- data.stand(envdata,method='standardize',margin='column',plot=T) #Then columnwise z-standardize
gauge_eucd<- vegdist(hydrostats_tra, method='euclidean') #Compute Euclidean distance
gaugecla_ward <-hclust(gauge_eucd, method='ward.D') #Classify using hierarchical agglomerative using Ward's minimum variance method

hclus.table(gaugecla_ward)
plot(gaugecla_ward)

