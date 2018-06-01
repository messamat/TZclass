#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created:05/31/2018
#Date last updated: 05/31/2018

#Purpose: classify Tanzanian river network using a deductive approach

# library(foreign)
# library(plyr)
# library(dplyr)
# library(hydrostats)
# library(data.table)
# #devtools::install_github("messamat/EflowStats") #Intro Vignette: https://cdn.rawgit.com/USGS-R/EflowStats/9507f714/inst/doc/intro.html
# #Corrected a glitch in package, need to re-change package download to USGS-R/EflowStats
# library(EflowStats)
# library(vegan) 
# library(pastecs)
# library(FD)
# library(cluster)
# library(pvclust)
# library(clusteval)
# library(adabag)
# library(rpart)
# library(rpart.plot)
# library(ggplot2)
# library(grid)
# library(gridExtra)
# library(ggdendro)
# library(dendextend)
# library(dendroextras)
# library(stringr)
# library(zoo)
# library(stargazer)

rootdir="F:/Tanzania/Tanzania" #####UPDATE THIS TO MATCH YOUR ROOT PROJECT FOLDER #######

source(file.path(rootdir,"bin/outside_src/Biostats.R"))

#Set folder structure
setwd(file.path(rootdir,"results")) 
origdatadir = file.path(rootdir,"data") 
outdir=file.path(getwd(),'TZ_deductive')
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}

#Import data
tzenv <- read.csv('streamnet118_final.csv')

##################### Format data ###################################
####Make subset of data/remove uneeded columns
colnames(tzenv)
outcols <- c(1:5,7,8,10,46:87,110:204,250:291,311:405, which(colnames(tzenv) %in% c('CatFlowAcc','wsFlowAcc','CatElvMin','CatDen','CatDamDen','CatFlowAcc','CatLCMaj',
                                                                     'WsPAPer','WsDamDen','WsGeolMaj','WsLCMaj','ReaElvMin',
                                                                     'ReaElvMax','SUM_LENGTH_GEO','Shape_Length','CatSoilg','CatPAPer','WsPAPer') |
                                              !is.na(str_match(colnames(tzenv),'DirSum*'))))
tzenvsub <- tzenv[,-outcols]
tzenvsub$ReaDirMaj <- as.factor(tzenvsub$ReaDirMaj)
#colnames(tzenvsub)

####Data transformation
str(tzenvsub)
colnames(tzenvsub)[factcol]

factcol <- c('GridID','ReaOrd','CatDirMaj','CatGeolMaj','CatSoilMaj','WsDirMaj','WsSoilMaj','ReaDirMaj') #Columns that should be considered as factor
#Make factor colums factors
tzenvsub[,factcol] <- sapply(tzenvsub[,factcol], as.factor) #Factorize 
#hist.plots(tzenvsub[,-factcol]) #Inspect data
logcols <- c('CatPopDen','ReaSloAvg','WsArea','WsPopDen','WsAIAvg') #Columns to be log-transformed
tzenvsub[,logcols] <- data.trans(data.frame(tzenvsub[,logcols]), method = 'log', plot = T)
sqrtcols <- c('CatAIAvg', 'CatBio13Av','CatBio14Av','CatBio16Av','CatBio17Av','CatBio18Av','CatBio19Av','CatElvMax','CatDRocAvg', 'CatElvAvg','CatSloAvg','CatSloStd',
              'CatEroAvg','CatPAPer','CatRoadDen','CatWatcha','CatMineDen','CatWatOcc','ReaPAPer','ReaElvAvg','WsBio13Av','WsBio14Av','WsBio17Av','WsBio19Av','WsElvMax',
              'WsElvAvg','WsEroAvg','WsSloAvg','WsSloStd','WsDen','WsRoadDen','WsWatcha','WsMineDen','WsWatOcc','WsWatSea','WsDRocAvg','ReaElvAvg') #Columns to be sqrt transform
tzenvsub[,sqrtcols][tzenvsub[,sqrtcols]<0] <- 0 #A few precipitation values are negative, correct back to 0
tzenvsub[,sqrtcols] <- data.trans(tzenvsub[,sqrtcols], method = 'power',exp=.5, plot = T)

###Standardization to mean of 0 and unit variance by variable
tzenvsub_std <- tzenvsub[,!which(factcol %in% colnames(tzenvsub))]
tzenvsub_std <- cbind(data.stand(tzenvsub_std, method = "standardize", margin = "column", plot = F),
                        tzenvsub[,factcol])
###Join standardized columns to gages
gagesenv_format <- gagesenv[,c('RGS_No','GridID')]
gagesenv_format <- merge(gagesenv_format,  tzenvsub_std, by='GridID')





