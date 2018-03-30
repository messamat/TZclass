#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/26/2018
#Date last updated: 03/29/2018

#Purpose: further assess flow record and filter out streamgages based on length of record, gaps, overlap, and non-stationarity
#N.B: here, filtering of gages based on environmental disturbance and non-stationarity is subsequent to analysis based on flow record length and overlap
# in order to be as inclusive as possible.

library(ggplot2)
library(data.table)
library(FlowScreen)
library(hydroTSM)
library(reshape)
library(plyr)
library(lemon)
library(foreign)
library(forecast)
library(imputeTS)
library(scales) # to access break formatting functions
library(grid)
library(gridExtra)
library(cowplot)
setwd("F:/Tanzania/Tanzania/results") #UPDATE
datadir = file.path(getwd(),paste('rufiji_hydrodatainspect','20180326',sep='_')) #UPDATE
origdatadir = "F:/Tanzania/Tanzania/data"
outdir=file.path(getwd(),'rufiji_hydrodatafilter_20180327')
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}

rufidat_clean <- read.csv(file.path(datadir,'rufidat_clean.csv'), colClasses=c('factor','Date','numeric','character','character'))
rufidat_deleted <- read.csv(file.path(datadir,'rufidat_deleted.csv'), colClasses=c('character','Date','numeric','character','character'))
rufienv <- read.dbf(file.path(getwd(),'streamnet118_rufiji_finaltab.dbf'))
gagesenv <- read.dbf(file.path(getwd(),'gages_netjoin.dbf'))
gagesenvrec <- merge(gagesenv, unique(rufidat_clean[,c('ID','SYM')]), by.x='RGS_No', by.y='ID', all.x=F)


##########################################
#BUild summary tables
rufidat_clean$year <- as.numeric(format(rufidat_clean$Date, "%Y"))
rufidat_clean$month <- as.numeric(format(rufidat_clean$Date, "%m"))
rufidat_dt <- data.table(rufidat_clean)

rufidat_datesummary <- rufidat_dt[,list(length(unique(ID))), .(Date)] #Compute number of gages with data for each date
rufidat_clean$ID <- factor(rufidat_clean$ID, levels = unique(rufidat_clean$ID[order(rufidat_clean$Date)]))

record_overview <-ggplot(data=rufidat_clean, aes(x=Date, y=ID)) +
  geom_point(size=2) +
  geom_bar(data=rufidat_datesummary, aes(x=Date,y='1KB36',color=V1), stat='identity') +
  geom_point(size=2) +
  scale_x_date(date_breaks ="5 years") +
  scale_colour_distiller(name='Number of gages',palette='Spectral') +
  theme_bw()
png(file.path(outdir,'record_overview.png'),width = 20, height=12,units='in',res=300)
print(record_overview)
dev.off()

#Compute length of gap between current daily record and previous record both within and between years for each date and gage
rufidat_dt[, prevdate := ifelse(shift(ID, 1L, type="lag")==ID,shift(as.character(Date), 1L, type="lag"), NA)]
rufidat_dt[, prevdateyr := ifelse((shift(ID, 1L, type="lag")==ID) & (shift(year, 1L, type="lag")==year),shift(as.character(Date), 1L, type="lag"), NA)]
rufidat_dt$prevgap <- as.numeric(rufidat_dt$Date-as.Date(rufidat_dt$prevdate))
rufidat_dt$prevgapyr <- as.numeric(rufidat_dt$Date-as.Date(rufidat_dt$prevdateyr))

###############################################
#Compute summary statistics in terms of amount of record, percentage of gaps, max length of gaps both
#regarding entire record
rufidat_summary <- rufidat_dt[,list(min_year=min(year), max_year=max(year), max_len=max(year)-min(year), 
                                    max_lend=max(Date)-min(Date), gap_d=as.numeric((max(Date)-min(Date)))-length(unique(Date)), gap_per=1-(length(unique(Date))/as.numeric((max(Date)-min(Date)))),
                                    max_gap = max(prevgap,na.rm=T))
                              ,.(ID)]
#regarding yearly record
rufidat_gapsummary <- rufidat_dt[,list(gap_d=as.numeric(format(as.Date(paste(year, "12", "31", sep="-")), "%j"))-length(unique(Date)),
                                       gap_per=1-(length(unique(Date))/as.numeric(format(as.Date(paste(year, "12", "31", sep="-")), "%j"))),
                                       max_gap = max(prevgapyr,na.rm=T))
                                 ,.(ID,year)]

#################################################
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

gapplot <- ggplot(rufidat_gapplot, aes(x=minyr, y=V1, color=as.factor(maxgap))) + 
  scale_x_continuous(name='Record length (years)', breaks=seq(5,50,5), expand=c(0,0)) +
  scale_y_continuous(name='Number of gages', limits=c(0,35), breaks=seq(0,35,5),expand=c(0,0))+
  scale_color_discrete(name="Maximum length of missing record per year (%)") +
  guides(color=guide_legend(ncol=4)) +
  geom_line(size=1) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        legend.position=c(0.65,0.85),
        legend.background = element_blank())
png(file.path(outdir,'gapplot.png'),width = 8, height=8,units='in',res=300)
print(gapplot)
dev.off()

############################################################
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
overlapplot <- ldply(seq(5,35,1), function(a) {
  ldply(seq(1955,2017,1), function(b) {
    ldply(seq(1,a)/a, function(c) {
      ldply(seq(5,50), function(d) {
        tryCatch(overlap_comp(period_len=a, cyr=b, completeness=c, minyr=d), error=function(e) NULL)
      })
    })
  })
})
colnames(overlapplot) <- c("period_len","cyr","completeness","minyr","count")
write.csv(overlapplot, file.path(outdir,'rufidat_overlapplot.csv'),row.names = F)
#Get maximum number of overlapping gages for each percentage of overlap,length of overlap period, and total length of record 
overlapplot <- as.data.table(overlapplot)
overlapplotmax <- overlapplot[, .SD[which.max(count)], .(period_len, completeness, minyr)]

#Plot
overlap_labels<-setNames(paste('Period of overlap:', seq(5,35,5),"years",sep=" "),seq(5,35,5))
overlapplot_out <-ggplot(overlapplotmax[overlapplotmax$period_len %in% seq(5,35,5),], aes(x=minyr, y=count, group=completeness, color=completeness)) +
  geom_line(size=1) +
  facet_wrap(~period_len, labeller=as_labeller(overlap_labels)) +
  scale_x_continuous(name='Record length (years)', breaks=seq(5,50,5), expand=c(0,0)) +
  scale_y_continuous(name='Number of gages', limits=c(0,35), breaks=seq(0,35,5),expand=c(0,0))+
  theme_bw() +
  scale_color_distiller(palette='Spectral', breaks=c(0,0.25,0.5,0.75,1),name='Minimum overlap (% of years)') +
  theme(legend.position="bottom",
        text=element_text(size=14)) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 20, barheight = 4))
png(file.path(outdir,'overlapplot.png'),width=20, height=12,units='in',res=300)
reposition_legend(overlapplot_out, 'left', panel='panel-3-3')
dev.off()

#Check best date of period based on a set of requirements
tryoverlap<-overlapplot[period_len>=15 & minyr>=15 & completeness>=0.3 & count>=15 & cyr>=2000 & cyr<=2018-15]
print(tryoverlap)

################################################
#Visualize correlation among gages
over5yr <- as.character(rufidat_gapyear[rufidat_gapyear$gagecount_gap_0.1>10,'ID']) #Only keep gages with at least 5 years of data (as otherwise correlation does not run)
rufidat_clean$Flowlog <- log(rufidat_clean$Flow+0.01)
rufidat_cast <- dcast(rufidat_clean[rufidat_clean$ID %in% over5yr,], Date~ID, value.var='Flowlog')
#png(file.path(outdir,'hydropairs.png'),width=40, height=40,units='in',res=300)
#hydropairs(rufidat_cast, dec=3, use="pairwise.complete.obs", method="pearson")
#dev.off()

#####################################################################
#Assess representativity of gages regarding environmental variables
#####################################################################

################################
#With individual variable plots
theme_env <- function () { 
  theme_bw(base_size=14) %+replace% 
    theme(
      panel.background  = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x=element_line(),
      axis.line.y=element_line()
    )
}
theme_envnoy <- function () { 
  theme_bw(base_size=14) %+replace% 
    theme(
      panel.background  = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x=element_line(),
      axis.line.y=element_line(),
      axis.title.y=element_blank()
    )
}

envplot_area <- ggplot(rufienv, aes(x=WsArea)) + 
  geom_vline(xintercept=gagesenvrec$WsArea, color='red',alpha=0.4, size=0.75) +
  geom_histogram(bins=50, fill='#878787', alpha=0.5) + 
  scale_x_log10(name=expression('Watershed area'~(km^3)),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand=c(0,0)) +
  scale_y_continuous(name='Number of reaches', expand=c(0,0)) + 
  theme_env()
#envplot_area
envplot_elv <- ggplot(rufienv, aes(x=ReaElvAvg)) + 
  geom_vline(xintercept=gagesenvrec$ReaElvAvg, color='red',alpha=0.4, size=0.75) +
  geom_histogram(bins=50,fill='#fdbf6f', alpha=0.65) + 
  scale_x_continuous(name=expression('River reach average elevation (m)'),
                     expand=c(0,0)) +
  scale_y_continuous(name='Number of reaches',expand=c(0,0)) + 
  theme_envnoy()
#envplot_elv
envplot_preci <- ggplot(rufienv, aes(x=WsBio12Av)) + 
  geom_vline(xintercept=gagesenvrec$WsBio12Av, color='red',alpha=0.4, size=0.75) +
  geom_histogram(bins=50,fill='#1f78b4', alpha=0.4) + 
  scale_x_continuous(name=expression('Upstream mean annual precipitation (mm)'),
                     expand=c(0,0)) +
  scale_y_continuous(name='Number of reaches',expand=c(0,0)) + 
  theme(axis.title.y=element_blank()) +
  theme_envnoy()
#envplot_preci


envplot_watext <- ggplot(rufienv, aes(x=100*CatWatExt)) + 
  geom_vline(xintercept=100*gagesenvrec$CatWatExt, color='red',alpha=0.4, size=0.75) +
  geom_histogram(bins=50,fill='#35978f', alpha=0.4) + 
  scale_x_continuous(name=expression('River catchment maximum water extent 1984-2015 (% of catchment area)'),
                     expand=c(0,0)) +
  scale_y_continuous(trans='sqrt',name='Number of reaches',expand=c(0,0)) + 
  theme_env()
#envplot_watext
rufienv$CatGeolMaj <- factor(rufienv$CatGeolMaj, levels = unique(rufienv$CatGeolMaj[order(as.numeric(as.character(rufienv$CatGeolMaj)))]))
envplot_geol <-  ggplot(rufienv, aes(x=CatGeolMaj)) + 
  geom_bar(fill='#8c510a', alpha=0.4, stat='count') +
  geom_bar(data=gagesenvrec, color='red',fill='white',alpha=0.4, size=0.75) +
  scale_x_discrete(name=expression('Main catchment lithology')) +
  scale_y_continuous(trans='log10',name='Number of reaches') + 
  theme_envnoy()
#envplot_geol
envplot_pop <- ggplot(rufienv, aes(x=CatPopDen)) + 
  geom_vline(xintercept=gagesenvrec$CatPopDen, color='red',alpha=0.4, size=0.75) +
  geom_histogram(bins=50,fill='#4d9221', alpha=0.4) + 
  scale_x_log10(name=expression('River catchment population density'~(persons/km^2)),
                breaks=c(1,5,10,50,100,500,1000,5000),
                expand=c(0,0)) +
  scale_y_continuous(name='Number of reaches',expand=c(0,0)) + 
  theme_envnoy()
#envplot_pop


envplot_forlos <- ggplot(rufienv, aes(x=100*WsFLosSum_)) + 
  geom_vline(xintercept=100*gagesenvrec$WsFLosSum_, color='red',alpha=0.4, size=0.75) +
  geom_histogram(bins=50,fill='#bf812d', alpha=0.6) + 
  scale_x_sqrt(name=expression('Upstream forest loss (% of watershed area)'),
               breaks=c(0,1,5,10,25,50,75,100),
               expand=c(0,0)) +
  scale_y_continuous(name='Number of reaches',expand=c(0,0)) + 
  theme_env()
#envplot_forlos
envplot_urb <- ggplot(rufienv, aes(x=100*LCSum_89)) + 
  geom_vline(xintercept=100*gagesenvrec$LCSum_89, color='red',alpha=0.4, size=0.75) +
  geom_histogram(bins=50,fill='#cab2d6', alpha=0.4) + 
  scale_x_continuous(trans='log10',
                     name='Upstream urban area (% of watershed area)',
                     breaks = c(0.001,0.01,0.1,1,10),
                     labels= c(0.001,0.01,0.1,1,10),
                     expand=c(0,0)) +
  scale_y_continuous(name='Number of reaches',expand=c(0,0)) + 
  theme_envnoy()
#envplot_urb
envplot_resind <- ggplot(rufienv, aes(x=WsResInd)) + 
  geom_vline(xintercept=gagesenvrec$WsResInd, color='red',alpha=0.4, size=0.75) +
  geom_histogram(bins=50,fill='#de77ae', alpha=0.4) + 
  scale_x_continuous(name=expression('Proportion of watershed draining from a reservoir (%)'),
                     expand=c(0,0)) +
  scale_y_log10(name='Number of reaches',expand=c(0,0),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  theme_envnoy()
#envplot_resind
#plot_grid(envplot_area, envplot_elv, envplot_preci, envplot_watext, envplot_geol, envplot_pop, envplot_forlos, envplot_urb, envplot_resind, align = "v", nrow = 3)

png(file.path(outdir,'gage_env.png'),width=20, height=12,units='in',res=300)
plot_grid(envplot_area, envplot_elv, envplot_preci, envplot_watext, envplot_geol, envplot_pop, envplot_forlos, envplot_urb, envplot_resind, align = "v", nrow = 3)
dev.off()

################################
#In multidimensional environment





#rufidat_overlapplot <- ldply(seq(5,50,1), function(y) {
#  rufidat_dtgap[ID %in% as.character(rufidat_gapyear[rufidat_gapyear$gagecount_gap_0.1>y,'ID']),list(length(unique(ID)), minyr=y), .(Date)]
#})


#####
#To do:
#Based on these requirements, assess number of gages based on non-stationarity
#Create plots of data showing deleted and interpolated data
#Evaluate consistency of discharge with drainage area and precipitation + HydroSHEDS modeled data
