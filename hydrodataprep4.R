#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/26/2018
#Date last updated: 03/27/2018

#Purpose: filter our streamgages based on length of record, gaps, overlap, and non-stationarity

library(ggplot2)
library(data.table)
library(FlowScreen) #to inspect data
library(reshape)
library(plyr)
setwd("F:/Tanzania/Tanzania/results") #UPDATE
datadir = file.path(getwd(),paste('rufiji_hydrodatainspect','20180326',sep='_')) #UPDATE
origdatadir = "F:/Tanzania/Tanzania/data"
outdir=file.path(getwd(),paste('rufiji_hydrodatafilter',as.character(format(Sys.Date(),'%Y%m%d')),sep='_'))
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}

rufidat_clean <- read.csv(file.path(datadir,'rufidat_clean.csv'), colClasses=c('factor','Date','numeric','character','character'))
rufidat_deleted <- read.csv(file.path(datadir,'rufidat_deleted.csv'), colClasses=c('character','Date','numeric','character','character'))

##########################################
#BUild summary tables
rufidat_clean$year <- as.numeric(format(rufidat_clean$Date, "%Y"))
rufidat_clean$month <- as.numeric(format(rufidat_clean$Date, "%m"))
rufidat_dt <- data.table(rufidat_clean)

rufidat_datesummary <- rufidat_dt[,list(length(unique(ID))), .(Date)] #Compute number of gages with data for each date
#record_overview <-ggplot(data=rufidat_clean, aes(x=Date, y=ID)) + 
#  geom_bar(data=rufidat_datesummary, aes(x=Date,y='1KB9',color=V1), stat='identity') +
#  geom_point() +
#  scale_x_date(date_breaks ="5 years") +
#  scale_colour_distiller('spectral') +
#  theme_bw()
#png(file.path(outdir,paste(gage,'record_overview.png',sep="_")),width = 20, height=12,units='in',res=300)
#print(record_overview)
#dev.off()

#Compute length of gap between current daily record and previous record both within and between years for each date and gage
rufidat_dt[, prevdate := ifelse(shift(ID, 1L, type="lag")==ID,shift(as.character(Date), 1L, type="lag"), NA)]
rufidat_dt[, prevdateyr := ifelse((shift(ID, 1L, type="lag")==ID) & (shift(year, 1L, type="lag")==year),shift(as.character(Date), 1L, type="lag"), NA)]
rufidat_dt$prevgap <- as.numeric(rufidat_dt$Date-as.Date(rufidat_dt$prevdate))
rufidat_dt$prevgapyr <- as.numeric(rufidat_dt$Date-as.Date(rufidat_dt$prevdateyr))

#Compute summary statistics in terms of amount of record, percentage of gaps, max length of gaps both
#regarding entire record
rufidat_summary <- rufidat_dt[,list(min_year=min(year), max_year=max(year), max_len=max(year)-min(year), 
                                    max_lend=max(Date)-min(Date), gap_d=as.numeric((max(Date)-min(Date)))-length(unique(Date)), gap_per=1-(length(unique(Date))/as.numeric((max(Date)-min(Date)))),
                                    max_gap = max(prevgap,na.rm=T))
                              ,.(ID)]
#regardomg yearly record
rufidat_gapsummary <- rufidat_dt[,list(gap_d=as.numeric(format(as.Date(paste(year, "12", "31", sep="-")), "%j"))-length(unique(Date)),
                                       gap_per=1-(length(unique(Date))/as.numeric(format(as.Date(paste(year, "12", "31", sep="-")), "%j"))),
                                       max_gap = max(prevgapyr,na.rm=T))
                                 ,.(ID,year)]

#Compute number of valid years on record depending on the percentage of missing data tolerated to consider a year valid
rufidat_gapyear <- data.frame(ID=unique(rufidat_gapsummary$ID))
for (i in seq(0.5,0,-0.05)) {
  df<-ddply(rufidat_gapsummary[rufidat_gapsummary$gap_per<=i,], .(ID), summarise, gcount=length(year))
  rufidat_gapyear <- merge(rufidat_gapyear,df, by='ID',all.x=T)
}
colnames(rufidat_gapyear) <- c('ID',paste('gagecount_gap', seq(0.5,0,-0.05),sep="_"))

#Compute for a range of durations, the number of stations that have data for at least this duration, for each maximum percentage gap threshold
rufidat_gapplot <- ldply(seq(5,50,1), function(y) {adply(rufidat_gapyear[,2:ncol(rufidat_gapyear)], 2, function(x) length(which(x>y)))})
rufidat_gapplot$minyr <- sort(rep(seq(5,50,1), ncol(rufidat_gapyear[,2:ncol(rufidat_gapyear)])))
rufidat_gapplot$maxgap <- as.numeric(substr(rufidat_gapplot$X1,15,18))
  
ggplot(rufidat_gapplot, aes(x=minyr, y=V1, color=as.factor(maxgap))) + 
  scale_x_continuous(name='Record length', breaks=seq(5,50,5), expand=c(0,0)) +
  scale_y_continuous(name='Number of gages', limits=c(0,35), breaks=seq(0,35,5),expand=c(0,0))+
  scale_color_discrete(name="Maximum percentage of missing record per year") +
  guides(color=guide_legend(ncol=4)) +
  geom_line(size=1) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        legend.position=c(0.65,0.85),
        legend.background = element_blank())

#Now compute the number of gages as a function of record length and record overlap given that we only keep years with less than 10% of missing data
max_yrgap <- 0.1
dtoverlap <-rufidat_gapsummary[rufidat_gapsummary$gap_per<=0.1,]
dtoverlap <- merge(dtoverlap, rufidat_gapyear[,c("ID","gagecount_gap_0.1")], by="ID")

#period_len=5
#cyr = 1964
#completeness=1
#minyr=10
#Compute for a given year, length of period, degree of overlap, and minimum length of total record, the number of gages
overlap_comp <- function(period_len, cyr, completeness, minyr) {
  return(c(period_len, cyr, completeness, minyr,
           length(which(dtoverlap[dtoverlap$gagecount_gap_0.1>=minyr & dtoverlap$year>=cyr & dtoverlap$year<cyr+period_len,.N/period_len,.(ID)][,2]>=completeness)))
  )
}
rufidat_overlapplot <- ldply(seq(5,35,1), function(a) {
  ldply(seq(1955,2017,1), function(b) {
    ldply(seq(1,a)/a, function(c) {
      ldply(seq(5,50), function(d) {
      tryCatch(overlap_comp(period_len=a, cyr=b, completeness=c, minyr), error=function(e) NULL)
        })
      })
    })
  })
rufidat_overlapplot (bla) <- c("period_len","cyr","completeness","minyr","count")

#Get maximum number of overlapping gages for each percentage of overlap,length of overlap period, and total length of record 



#rufidat_overlapplot <- ldply(seq(5,50,1), function(y) {
#  rufidat_dtgap[ID %in% as.character(rufidat_gapyear[rufidat_gapyear$gagecount_gap_0.1>y,'ID']),list(length(unique(ID)), minyr=y), .(Date)]
#})


#####
#To do:
#Based on these requirements, assess number of gages based on non-stationarity
#Evaluate consistency of discharge with drainage area and precipitation + HydroSHEDS modeled data
#See hyoSTM Hydropairs
#Interpolate data
#