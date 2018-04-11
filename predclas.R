#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 04/11/2018
#Date last updated: 04/02/2018

#Purpose: interpolate hydrological metrics to entire network (predict) then classify network

library(foreign)
library(plyr)
library(dplyr)
library(hydrostats)
library(data.table)
#devtools::install_github("messamat/EflowStats") #Intro Vignette: https://cdn.rawgit.com/USGS-R/EflowStats/9507f714/inst/doc/intro.html
#Corrected a glitch in package, need to re-change package download to USGS-R/EflowStats
library(EflowStats)
source(file.path(rootdir,"bin/outside_src/Biostats.R"))
library(vegan) 
library(pastecs)
library(FD)
library(cluster)
library(pvclust)
library(clusteval)
library(adabag)
library(rpart)
library(rpart.plot)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggdendro)
library(dendextend)
library(dendroextras)

rootdir="F:/Tanzania/Tanzania" #UPDATE
source(file.path(rootdir,"bin/outside_src/Biostats.R"))
source(file.path(rootdir,"bin/outside_src/Flowscreen.hyear.internal.R"))

setwd(file.path(rootdir,"results")) 
datadir = file.path(getwd(),'rufiji_hydrodatastats') #UPDATE
origdatadir = file.path(rootdir,"data") 
outdir=file.path(getwd(),'rufiji_hydrodatapredclas')
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}

gagesenv <- read.csv(file.path(getwd(),'gages_netjoinclean.csv'))
rufienv <- read.csv(file.path(getwd(),'streamnet118_rufiji_finaltabclean.csv'))
HITo15y <- read.csv(file.path(datadir, 'HITo15y.csv'))
HITo15y <- dcast(setDT(HITo15y), ID~indice)

