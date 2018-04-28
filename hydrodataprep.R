#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/15/2018
#Date last updated: 03/24/2018

#Purpose: import and merge hydrological data for the Rufiji basin of Tanzania — provided by CDMSmith Zachary T. Eichenwald (ZTE) and Japhet Kashaigili (JK)

devtools::install_github("awalker89/openxlsx")
library(openxlsx)
library(reshape2)
library(DescTools)
library(stringr)
rootdir= "F:/Tanzania/Tanzania" #UPDATE
setwd(file.path(rootdir, "results")) 
datadir = file.path(rootdir, "data")
########################################
#Import and merge hydrological data
########################################
###Import and merge ruaha data from Zach 
setClass('myDate')
setAs("character","myDate", function(from)  as.POSIXct(from, format= "%m/%d/%Y %H:%M")) #Create format to directly import data column in appropriate format
ruahaflow <- read.csv(file.path(datadir,"sharepoint20180316/flow/2018-03-06 Corrected Ruaha Stage and Flow_data.csv"),
                      colClasses = c('factor','factor','myDate','numeric','numeric','character','character'))
length(table(ruahaflow$Gage.ID))
ruahaflow <- ruahaflow[,!names(ruahaflow) %in% c("X","X.1")]
#str(ruahaflow)
ruahastations <- read.csv(file.path(datadir,"sharepoint20180316/flow/2018-03-06 Corrected Ruaha Stage and Flow_stations.csv"))
#Check whether comments on record and rating curve quality are unique for each gage
nrow(ruahastations)
#str(ruahastations)
ruahadat <- merge(ruahaflow, ruahastations[,c('Station_ID','Rating.Curve.Source','Record.Quality','Rating.Curve.Quality')], 
                  by.x='Gage.ID', by.y='Station_ID')
str(ruahadat)

###Import and merge Kilombero data (some data are in "%m/%d/%Y %H:%M" and others are in "%m/%d/%Y")
kilomflow <- read.csv(file.path(datadir,"sharepoint20180316/flow/2018-03-13 Corrected Kilombero Stage and Flow_data.csv"),colClasses = c('factor','factor','myDate','numeric','numeric','character','character','character','character','character'))
kilomflow2 <- read.csv(file.path(datadir,"sharepoint20180316/flow/2018-03-13 Corrected Kilombero Stage and Flow_data.csv"),colClasses = c('factor','factor','character','numeric','numeric','character','character','character','character','character'))
kilomflow[is.na(kilomflow$Date.Time),'Date.Time'] <- as.POSIXct(paste(kilomflow2[is.na(kilomflow$Date.Time),'Date.Time'], "12:00:00"), format= "%m/%d/%Y %H:%M:%S") 
kilomflow <- kilomflow[,!names(kilomflow) %in% c("X","X.1","X.2","X.3","X.4")]
str(kilomflow)
colnames(kilomflow) <- colnames(ruahaflow)
length(table(kilomflow$Gage.ID))
kilomstations <- read.csv(file.path(datadir,"sharepoint20180316/flow/2018-03-13 Corrected Kilombero Stage and Flow_stations.csv"))
#Check whether comments on record and rating curve quality are unique for each gage
nrow(kilomstations)
#Comments are not unique. Go back to original data: it seems like comments for first 1KB36 in stations
#"1973 RC good. 2015-2017 RC poor, with little data." in fact correspond to comments on station 1KB27 Luipa/Ruipa river
#Observed on OneNote file from ZTE. However, it does not appear that there are data for 1KB27 in the dataset? To inquire with ZTE.
kilomstations <- kilomstations[!(kilomstations$Station_ID=='1KB36' & kilomstations$Record.Quality=='Poor, many data gaps'),]
str(kilomstations)
kilomdat <- merge(kilomflow, kilomstations[,c('Station_ID','Rating.Curve.Source','Record.Quality','Rating.Curve.Quality')], by.x='Gage.ID', by.y='Station_ID')

###Merge great Ruaha data and Kilombero data
rufidat <- rbind(ruahadat, kilomdat)

##########################################################################################################################################
###Import raw stage data from RBWB given by ZTE
setAs("character","myDate", function(from)  as.POSIXct(from, format= "%m/%d/%Y")) #Create format to directly import data column in appropriate format
RBWBhdtb <- read.csv(file.path(datadir,"sharepoint20180316/flow/2018-03-06_HydroDatabase_StageMeasurementDatabase.csv"),
                     colClasses = c('factor','factor','character','character','myDate',rep('character',13)))

RBWBhdtb[,"Date.Time"] <- as.POSIXct(paste(RBWBhdtb$Date.Time, RBWBhdtb$Time), format="%Y-%m-%d %H:%M:%S")
RBWBhdtb <- RBWBhdtb[,!(names(RBWBhdtb) %like% "X%")]
RBWBhdtb <- RBWBhdtb[!(RBWBhdtb$Stage..m.=='m'),] 
RBWBhdtb$Stage..m. <- as.numeric(RBWBhdtb$Stage..m.)

##########################################################################################################################################
###Import flow data from RBWB David
RBWBflow <- data.frame(Date=NULL, `Flow`=NULL, Gage.ID=NULL)
for (sheet in 1:25) {
  flowdf <- read.xlsx(file.path(datadir,"RBWB_David20180415/Sheet 2_ flows_edit.xlsx"),sheet=sheet,startRow=1, detectDates=F)
  gage <- colnames(flowdf)[2]
  print(gage)
  flowdf <- read.xlsx(file.path(datadir,"RBWB_David20180415/Sheet 2_ flows_edit.xlsx"),sheet=sheet,startRow=3, detectDates=F)
  colnames(flowdf) <- c('Date', 'Flow')
  flowdf$Gage.ID <- gage
  flowdf$Date <- as.Date(flowdf$Date, origin='1899-12-30')
  RBWBflow <- rbind(RBWBflow, flowdf)
}


##########################################################################################################################################
###Import daily data from Japhet K.
Japhet_GRuaha <- read.xlsx(file.path(datadir,"JaphetK_20180316/Rufiji_L.RukwaBasins_Flow Data_formatMM20180322.xlsx"),sheet=1,startRow=4, detectDates=T)
Japhet_1KA59 <- Japhet_GRuaha[,c(1,2)]
Japhet_1KA59[,'Station'] <- '1KA59'
Japhet_1KA11 <- Japhet_GRuaha[,c(3,4)]
Japhet_1KA11[,'Station'] <- '1KA11'
Japhet_GRuaha <- rbind(Japhet_1KA59, Japhet_1KA11)
Japhet_GRuaha <- Japhet_GRuaha[!is.na(Japhet_GRuaha$Date) & Japhet_GRuaha$`Flow.(m3/s)`!='m' & Japhet_GRuaha$`Flow.(m3/s)`!='-',]
Japhet_GRuaha$`Flow.(m3/s)` <- as.numeric(Japhet_GRuaha$`Flow.(m3/s)`)
str(Japhet_GRuaha)

Japhet_Kilom <- read.xlsx(file.path(datadir,"JaphetK_20180316/Rufiji_L.RukwaBasins_Flow Data_formatMM20180322.xlsx"),sheet=2,startRow=3, detectDates=T)
Japhet_Kilom <- Japhet_Kilom[-1,]
colnames(Japhet_Kilom)[c(1,6)] <- c('Date','Udagaji')
Japhet_Kilom <- melt(Japhet_Kilom, id.vars = 'Date',value.name='Flow.(m3/s)',variable.name='Station')
str(Japhet_Kilom)
colnames(Japhet_Kilom) <- colnames(Japhet_GRuaha)[c(1,3,2)]
Japhet_Kilom$`Flow.(m3/s)` <- as.numeric(Japhet_Kilom$`Flow.(m3/s)`)
Japhet_Kilom$Date <- as.Date(Japhet_Kilom$Date)
Japhet_Kilom <- Japhet_Kilom[!is.na(Japhet_Kilom$`Flow.(m3/s)`),]

Japhet_LRuaha <- read.xlsx(file.path(datadir,"JaphetK_20180316/Rufiji_L.RukwaBasins_Flow Data_formatMM20180322.xlsx"),sheet=3,startRow=3, detectDates=T)
Japhet_LRuaha <- Japhet_LRuaha[-1,]
colnames(Japhet_LRuaha)[1] <- 'Date'
colnames(Japhet_LRuaha)[4] <- '1KA31'
Japhet_LRuaha <- melt(Japhet_LRuaha, id.vars = 'Date',value.name='Flow.(m3/s)',variable.name='Station')
colnames(Japhet_LRuaha) <- colnames(Japhet_GRuaha)[c(1,3,2)]
Japhet_LRuaha$`Flow.(m3/s)` <- as.numeric(Japhet_LRuaha$`Flow.(m3/s)`)
Japhet_LRuaha$Date <- as.Date(Japhet_LRuaha$Date)
str(Japhet_LRuaha)
Japhet_LRuaha <- Japhet_LRuaha[!is.na(Japhet_LRuaha$`Flow.(m3/s)`),]

Japhet_dailydat <- rbind(Japhet_GRuaha, Japhet_Kilom, Japhet_LRuaha)

#Import monthly unimpaired data from Japhet K.
Japhet_Unimpaired <- read.xlsx(file.path(datadir,"JaphetK_20180316/Rufiji_L.RukwaBasins_Flow Data_formatMM20180322.xlsx"),sheet=4,startRow=2, detectDates=T)
Japhet_Unimpaired <- Japhet_Unimpaired[-1,]
Japhet_Unimpaired <- melt(Japhet_Unimpaired, id.vars = 'Date',value.name='Flow.(m3/s)',variable.name='Station')
str(Japhet_Unimpaired)
colnames(Japhet_Unimpaired) <- colnames(Japhet_GRuaha)[c(1,3,2)]
Japhet_Unimpaired$`Flow.(m3/s)` <- as.numeric(Japhet_Unimpaired$`Flow.(m3/s)`)
Japhet_Unimpaired <- Japhet_Unimpaired[!is.na(Japhet_Unimpaired$`Flow.(m3/s)`),] #No NA, must be modeled?

#Import monthly Rukwa data from Japhet K.
Japhet_Rukwa <- read.xlsx(file.path(datadir,"JaphetK_20180316/Rufiji_L.RukwaBasins_Flow Data_formatMM20180322.xlsx"),sheet=5,startRow=2, detectDates=T)
Japhet_Rukwa <- Japhet_Rukwa[-1,c(-3,-5)]
colnames(Japhet_Rukwa) <- c('Date','3CD2','3B2','3A17')
Japhet_Rukwa <- melt(Japhet_Rukwa, id.vars = 'Date',value.name='Flow.(m3/s)',variable.name='Station')
colnames(Japhet_Rukwa) <- colnames(Japhet_GRuaha)[c(1,3,2)]
Japhet_Rukwa$`Flow.(m3/s)` <- as.numeric(Japhet_Rukwa$`Flow.(m3/s)`)
Japhet_Rukwa <- Japhet_Rukwa[!is.na(Japhet_Rukwa$`Flow.(m3/s)`),] #No NA, must be modeled?

########################################
#Export data to directory
########################################
outdir=file.path(getwd(),'rufiji_hydrodataraw')
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}
write.csv(rufidat, file.path(outdir,'ZTE_rufidat.csv'),row.names = F)
write.csv(RBWBhdtb, file.path(outdir, 'ZTE_RBWBhydrodat.csv'), row.names=F)
write.csv(RBWBflow, file.path(outdir, 'RBWBflowdat.csv'), row.names=F)
write.csv(Japhet_dailydat, file.path(outdir,'JK_dailydat.csv'),row.names = F)
write.csv(Japhet_Unimpaired, file.path(outdir,'JK_monthly_unimpaired.csv'),row.names = F)
write.csv(Japhet_Rukwa, file.path(outdir,'JK_monthly_rukwa.csv'),row.names = F)