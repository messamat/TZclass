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
library(reshape)
library(plyr)
library(lemon)
library(foreign) #For impot/export of dbf
library(scales) # to access break formatting functions
library(grid)
library(gridExtra) #For multi plot export
library(stringr) #For string matching
library(pastecs) #For data preparation
library(vegan) #For data preparation
library(FD) #For gower's dist
#library(fifer) #for stratified sampling

rootdir <- rprojroot::find_root(rprojroot::has_dir("src"))
source(file.path(rootdir,"bin/outside_src/Biostats.R"))
source(file.path(rootdir,"bin/outside_src/Flowscreen.hyear.internal.R"))
setwd(file.path(rootdir,"results")) #UPDATE
datadir = file.path(getwd(),'rufiji_hydrodatainspect') #UPDATE
origdatadir = "F:/Tanzania/Tanzania/data"
outdir=file.path(getwd(),'rufiji_hydrodatafilter')
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}

#Function to extract legend from graph as a grob to be re-inserted later
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#Import and format QA/QCed data
rufidat_clean <- read.csv(file.path(datadir,'rufidat_clean.csv'), colClasses=c('factor','Date','numeric','character','character'))
rufidat_clean$year <- as.numeric(format(rufidat_clean$Date, "%Y"))
rufidat_clean$month <- as.numeric(format(rufidat_clean$Date, "%m"))
rufidat_clean<-hyear.internal(rufidat_clean,hyrstart=10) #Ignore hdoy

#Import deleted records
rufidat_deleted <- read.csv(file.path(datadir,'rufidat_deleted.csv'), colClasses=c('character','Date','numeric','character','character'))

#Import rufiji river network environmental data
# rufienv <- read.csv(file.path(getwd(),'streamnet118_rufiji_finaltab.csv'))
# numcol <- colnames(rufienv)[which(sapply(rufienv, is.numeric))]
# sumcol <-adply(rufienv[,numcol],2,function(x) sum(x, na.rm=T))
# incol <- colnames(rufienv)[!(colnames(rufienv) %in% sumcol[sumcol$V1==0,'X1'])]
# rufienv <- rufienv[,incol] #Take out all 0 columns
# rufienv[,'WsVegPer'] <-  with(rufienv, LCSum_12+LCSum_23+LCSum_34) #Sum % trees, scrubs, grassland into a vegetation variable
# rufienv[,'WsAgriPer'] <-  with(rufienv, LCSum_45+LCSum_78) #Sum % cropland and bare ground into an agriculture variable
# rufienv[rufienv$WsVegPer>1,'WsVegPer']<- 1 #Some rounding must have led to insignificant exceedance of 1
# rufienv[rufienv$WsAgriPerr>1,'WsAgriPer'] <- 1
# write.csv(rufienv, file.path(getwd(),'streamnet118_rufiji_finaltabclean.csv'),row.names=F)
rufienv <- read.csv(file.path(getwd(),'streamnet118_rufiji_finaltabclean.csv'))

#Import gages environmental data
# gagesenv <- read.csv(file.path(getwd(),'gages_netjoin.csv'))
# incol<-colnames(gagesenv)[!(colnames(gagesenv) %in% sumcol[sumcol$V1==0,'X1'])] #Take out columns with only 0 values
# gagesenv <- gagesenv[,incol] #Take out all 0 columns
# gagesenv[,'WsVegPer'] <-  with(gagesenv, LCSum_12+LCSum_23+LCSum_34) #Sum % trees, scrubs, grassland into a vegetation variable
# gagesenv[,'WsAgriPer'] <-  with(gagesenv, LCSum_45+LCSum_78) #Sum % cropland and bare ground into an agriculture variable
# gagesenv[gagesenv$WsVegPer>1,'WsVegPer']<- 1 #Some rounding must have led to insignificant exceedance of 1
# gagesenv[gagesenv$WsAgriPerr>1,'WsAgriPer'] <- 1
# write.csv(gagesenv, file.path(getwd(),'gages_netjoinclean.csv'),row.names=F)
gagesenv <- read.csv(file.path(getwd(),'gages_netjoinclean.csv'))
gagesenvrec <- merge(gagesenv, unique(rufidat_clean[,c('ID','SYM')]), by.x='RGS_No', by.y='ID', all.x=F)
#write.csv(gagesenvrec, file.path(getwd(),'maps/gageenvrec_20180515.csv'),row.names=F)

#Import and format interpolated hydrological data
impute_preds <- read.csv(file.path('rufiji_hydrodataimpute', 'rufidat_interp.csv'), colClasses=c('Date',rep(c('numeric','numeric','character'),39)))
colnames(impute_preds)[seq(2,ncol(impute_preds),3)] <- substr(colnames(impute_preds),2,10)[seq(2,ncol(impute_preds),3)]

predsmelt <-melt.data.table(as.data.table(
  as.data.frame(impute_preds)[,
                              c(1,which(colnames(impute_preds) %like% "^1K*"))]),
  id.vars = 'Date',value.name='Flow',variable.name='ID') %>%
  .[, ID := gsub('^X', '', ID)]
  
predsmelt <- predsmelt[,c(2,1,3)]
predsmelt$SYM <- NA
predsmelt$Agency <- NA
predsmelt <- predsmelt[!is.na(predsmelt$Flow),]
predsmelt$year <- as.numeric(format(predsmelt$Date, "%Y"))
predsmelt$month <- as.numeric(format(predsmelt$Date, "%m"))
predsmelt<-hyear.internal(predsmelt,hyrstart=10) #hdoy doesn't work
#Compute hydrologic day
predsmelt$doy <- as.numeric(format(predsmelt$Date,"%j"))
predsmelt[,hdoy:=ifelse(month>=10,
                        doy-as.numeric(format(as.Date(paste(year,'-10-01',sep="")),"%j")),
                        doy+as.numeric(format(as.Date(paste(year-1,'-12-31',sep="")),"%j"))-as.numeric(format(as.Date(paste(year-1,'-10-01',sep="")),"%j")))] 

##########################################
#Build summary tables
rufidat_dt <- data.table(rufidat_clean)

rufidat_datesummary <- rufidat_dt[,list(length(unique(ID))), .(Date)] #Compute number of gages with data for each date
rufidat_clean$ID <- factor(rufidat_clean$ID, levels = unique(rufidat_clean$ID[order(rufidat_clean$Date)]))

#Compute length of gap between current daily record and previous record both within and between years for each date and gage
rufidat_dt[, prevdate := ifelse(data.table::shift(ID, 1L, type="lag")==ID, 
                                data.table::shift(as.character(Date), 1L, type="lag"), NA)]
rufidat_dt[, prevdateyr := ifelse((data.table::shift(ID, 1L, type="lag")==ID) & (data.table::shift(hyear, 1L, type="lag")==hyear),
                                  data.table::shift(as.character(Date), 1L, type="lag"), NA)]
rufidat_dt[, nextdate := ifelse(data.table::shift(ID, 1L, type="lead")==ID,
                                data.table::shift(as.character(Date), 1L, type="lead"), NA)]
rufidat_dt[, nextdateyr := ifelse((data.table::shift(ID, 1L, type="lead")==ID) & (data.table::shift(hyear, 1L, type="lead")==hyear),
                                  data.table::shift(as.character(Date), 1L, type="lead"), NA)]

rufidat_dt$gap <- as.numeric(as.Date(rufidat_dt$nextdate)-as.Date(rufidat_dt$prevdate))-2
rufidat_dt$gapyr <- as.numeric(as.Date(rufidat_dt$nextdateyr)-as.Date(rufidat_dt$prevdateyr))-2

############################################### summary statistics###################
#Compute summary statistics in terms of amount of record, percentage of gaps, max length of gaps both
#regarding entire record
rufidat_summary <- rufidat_dt[,list(min_hyear=min(hyear), max_hyear=max(hyear), max_len=max(hyear)-min(hyear), 
                                    max_lend=max(Date)-min(Date), gap_d=as.numeric((max(Date)-min(Date)))-length(unique(Date)), 
                                    gap_per=1-(length(unique(Date))/as.numeric((max(Date)-min(Date)))),
                                    max_gap = max(gap,na.rm=T))
                              ,.(ID)]
#regarding yearly record
rufidat_gapsummary <- rufidat_dt[,list(gap_d=as.numeric(format(as.Date(paste(hyear, "12", "31", sep="-")), "%j"))-length(unique(Date)), #Check the number of days in that year to account for leap years
                                       gap_per=1-(length(unique(Date))/as.numeric(format(as.Date(paste(hyear, "12", "31", sep="-")), "%j"))),
                                       max_gap = max(gap,na.rm=T), max_gapyr = max(gapyr,na.rm=T))
                                 ,.(ID,hyear)]
write.csv(rufidat_gapsummary, file.path(outdir, 'rufidat_gapsummary.csv'), row.names=F)

############################## FIGURE 4-3 ##############################################
############################## Overall figure of record (record_overview)###############
rufidat_cleanenv <- merge(rufidat_clean, gagesenv, by.x='ID', by.y='RGS_No')
rufidat_cleanenv$label <- as.factor(paste(rufidat_cleanenv$RGS_Loc,"R. at", rufidat_cleanenv$RGS_Name,"-", rufidat_cleanenv$ID, sep=" "))
rufidat_cleanenv$label <- factor(rufidat_cleanenv$label, levels = unique(rufidat_cleanenv$label[order(rufidat_cleanenv$Date)]))
rufidat_cleanenv$labelnum <-as.numeric(rufidat_cleanenv$label)
rufidat_datesummary$label <- factor('Mgugwe R. at Mgugwe - 1KB36', levels = unique(rufidat_cleanenv$label[order(rufidat_cleanenv$Date)]))
rufidat_datesummary$labelnum<- max(rufidat_cleanenv$labelnum)

record_overview_name <- ggplot(data=rufidat_cleanenv, aes(x=Date, y=labelnum)) +
  geom_bar(data=rufidat_datesummary, aes(color=V1), stat='identity') +
  geom_point(data=rufidat_cleanenv, size=2) +
  scale_x_date(breaks=as.Date(paste(c(1954,seq(1955,2015,5), 2017),'-01-01',sep="")), expand=c(0,0), date_labels = "%Y") +
  scale_y_continuous(name='River gauge name (format: River at Location - ID)',
                     breaks=1:max(rufidat_cleanenv$labelnum),
                     labels= unique(rufidat_cleanenv[order(rufidat_cleanenv$labelnum), 'label'])) + 
  scale_colour_distiller(name='No. of river gauges',palette='Spectral',breaks=c(5,10,15,20,25,max(rufidat_datesummary$V1)),
                         limits=c(min(rufidat_datesummary$V1),max(rufidat_datesummary$V1))) +
  theme_classic() +
  coord_cartesian(expand=FALSE, clip='off') +
  theme(text=element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.y = element_text(vjust=-10),
        legend.key.size = unit(3,"line"),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
png(file.path(outdir,'record_overview20200628_names.png'),width=20, height=12,units='in',res=600)
print(record_overview_name)
dev.off()

############################## FIGURE 4-4 ##############################################
############################## Gap plot ################################################
#Compute number of valid years on record depending on the percentage of missing data tolerated to consider a year valid and the corresponding
#maximum gap.
rufidat_gapyear <- data.frame(ID=unique(rufidat_gapsummary$ID))
for (i in seq(0.5,0,-0.1)) {
  df<-ddply(rufidat_gapsummary[rufidat_gapsummary$gap_per<=i & rufidat_gapsummary$max_gap<=(i*365.25),], .(ID), summarise, gcount=length(hyear))
  rufidat_gapyear <- merge(rufidat_gapyear,df, by='ID',all.x=T)
}
colnames(rufidat_gapyear) <- c('ID',paste('gagecount_gap', seq(0.5,0,-0.1),sep="_"))

#Compute for a range of durations, the number of stations that have data for at least this duration, for each percentage gap threshold
rufidat_gapplot <- ldply(seq(1,50,1), function(y) {
  adply(rufidat_gapyear[,2:ncol(rufidat_gapyear)], 2, function(x) length(which(x>y)))})
rufidat_gapplot$minyr <- sort(rep(seq(1,50,1), ncol(rufidat_gapyear[,2:ncol(rufidat_gapyear)])))
rufidat_gapplot$threshgap <- as.numeric(substr(rufidat_gapplot$X1,15,18))

gapplot <- ggplot(rufidat_gapplot, aes(x=minyr, y=V1, color=as.factor(100-100*threshgap))) + 
  scale_x_continuous(name='Record length (years)', breaks=c(1,seq(5,60,5)), expand=c(0,0)) +
  scale_y_continuous(name='Number of river gauges', limits=c(0,40), breaks=seq(0,40,5),expand=c(0,0))+
  scale_color_discrete(name="Minimum yearly record completeness (% of days)") +
  guides(color=guide_legend(ncol=4)) +
  geom_line(size=1.5) +
  theme_bw() + 
  theme(text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position=c(0.65,0.85),
        legend.background = element_blank(),
        panel.grid.minor=element_blank(),
        legend.key.size = unit(3,"line"),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
png(file.path(outdir,'gapplot20200628.png'),width = 12, height=12,units='in',res=600)
print(gapplot)
dev.off()

############REPORT
rufidat_gapplot[rufidat_gapplot$minyr==15 & rufidat_gapplot$threshgap==0,]
rufidat_gapplot[rufidat_gapplot$minyr==15 & rufidat_gapplot$threshgap==0.5,]

##################################################Overlapplot############
#Now compute the number of gages as a function of record length and record overlap given that we only keep years with less than 10% of missing data
max_yrgap <- 0.1
dtoverlap <-rufidat_gapsummary[rufidat_gapsummary$gap_per<=0.1,]
dtoverlap <- merge(dtoverlap, rufidat_gapyear[,c("ID","gagecount_gap_0.1")], by="ID")

period_len=5
cyr = 1964
completeness=1
minyr=10
#Compute for a given year, length of period, degree of overlap, and minimum length of total record, the number of gages
overlap_comp <- function(dtoverlap, period_len, cyr, completeness, minyr) {
  return(c(period_len, cyr, completeness, minyr,
           length(which(dtoverlap[dtoverlap$gagecount_gap_0.1>=minyr & dtoverlap$hyear>=cyr & dtoverlap$hyear<cyr+period_len,.N/period_len,.(ID)][,2]>=completeness)))
  )
}

dtoverlap[dtoverlap$gagecount_gap_0.1>=minyr & dtoverlap$hyear>=cyr & dtoverlap$hyear<cyr+period_len,.N/period_len,.(ID)]
which(dtoverlap[dtoverlap$gagecount_gap_0.1>=minyr & dtoverlap$hyear>=cyr & dtoverlap$hyear<cyr+period_len,.N/period_len,.(ID)][,2]>=completeness)

overlapplot <- ldply(seq(5,35,1), function(a) {
  ldply(seq(1955,2017,1), function(b) {
    ldply(seq(1,a)/a, function(c) {
      ldply(seq(5,50), function(d) {
        tryCatch(overlap_comp(dtoverlap, period_len=a, cyr=b, completeness=c, minyr=d), error=function(e) NULL)
      })
    })
  })
})
colnames(overlapplot) <- c("period_len","cyr","completeness","minyr","count")
write.csv(overlapplot, file.path(outdir,'rufidat_overlapplot.csv'),row.names = F)
overlapplot <- read.csv(file.path(outdir,'rufidat_overlapplot.csv'))
#Get maximum number of overlapping gages for each percentage of overlap,length of overlap period, and total length of record 
overlapplot <- as.data.table(overlapplot)
overlapplotmax <- overlapplot[, .SD[which.max(count)], .(period_len, completeness, minyr)]

#Plot
overlap_labels<-setNames(paste('Period of overlap:', seq(5,35,5),"hydrologic years",sep=" "),seq(5,35,5))
overlapplot_out <-ggplot(overlapplotmax[overlapplotmax$period_len %in% seq(5,35,5),], aes(x=minyr, y=count, group=completeness, color=completeness)) +
  geom_line(size=1) +
  facet_wrap(~period_len, labeller=as_labeller(overlap_labels)) +
  scale_x_continuous(name='Record length (hydrologic years)', breaks=seq(5,50,5), expand=c(0,0)) +
  scale_y_continuous(name='Number of gages', limits=c(0,35), breaks=seq(0,35,5),expand=c(0,0))+
  theme_bw() +
  scale_color_distiller(palette='Spectral', breaks=c(0,0.25,0.5,0.75,1),name='Minimum overlap (% of years)') +
  theme(legend.position="bottom",
        text=element_text(size=14)) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 20, barheight = 4))
png(file.path(outdir,'overlapplot.png'),width=20, height=12,units='in',res=600)
reposition_legend(overlapplot_out, 'left', panel='panel-3-3')
dev.off()

#Check best date of period based on a set of requirements
tryoverlap<-overlapplot[period_len>=30 & minyr>=15 & completeness>=0.25 & count>=15 & cyr>=1990 & cyr<=2018-15]
print(tryoverlap)

#####################################################################
#Compute number of years of data for later subsetting of gages (USED IN FIGURE 4-1) 
#####################################################################
#All
ycount_all <-  rufidat_gapsummary[,length(unique(hyear)),.(ID)]
#< 10% missing data, full length of record
ycount_full <- rufidat_gapsummary[gap_per<=0.1 & max_gap<37,length(unique(hyear)),.(ID)]
#< 10% missing data, pre-1983
ycount_pre83 <- rufidat_gapsummary[gap_per<=0.1 & max_gap<37 & hyear>1958 & hyear<=1983,length(unique(hyear)),.(ID)]
#< 10% missing data, post-2001
ycount_post01 <- rufidat_gapsummary[gap_per<=0.1 & max_gap<37 & hyear>1991 & hyear<=2016,length(unique(hyear)),.(ID)]
rufidat_ycount <- merge(ycount_full, ycount_pre83, by='ID',all.x=T)
rufidat_ycount <- merge(rufidat_ycount, ycount_post01, by='ID',all.x=T)
rufidat_ycount[is.na(rufidat_ycount)] <- 0
colnames(rufidat_ycount) <- c('ID','ycount_full', 'ycount_pre83', 'ycount_post01')
write.csv(rufidat_ycount, file.path(outdir, 'rufidat_ycount.csv'), row.names=F)

####
rufidat_clean <- merge(rufidat_clean, rufidat_gapsummary, by=c('ID','hyear'),all.x=T)
#rufidat_clean <- merge(rufidat_clean, rufidat_post1991, by='ID',all.x=T)

#Final set of gauges
#Select subset of gauges: include 1KB32 even if only 14 years of data + Kisigo stations but skip 1KB28 with simulated data
predsmelt <- merge(predsmelt, rufidat_ycount, by='ID',all.x=T)
predsmelt <- merge(predsmelt, rufidat_gapsummary, by=c('ID','hyear'),all.x=T)
rufidat_select_o15y <- predsmelt[predsmelt$gap_per<=0.1 & predsmelt$max_gap < 37 & predsmelt$hyear<2017 & predsmelt$ID !='1KB28' & 
                                   (predsmelt$ycount_full>=15 | predsmelt$ID=='1KB32') | 
                                   (predsmelt$ID=='1KA41' & predsmelt$max_gap < 273 & predsmelt$gap_per<0.75 & predsmelt$hyear<1996) |
                                   (predsmelt$ID=='1KA42A' & predsmelt$max_gap < 273 & predsmelt$gap_per<0.75 & predsmelt$hyear>1958 & predsmelt$hyear<2017),]


#######################################################################
#Plot clean and interpolated data (ANNEX B)
#######################################################################
##################################Get example legend ######
gage='1KA9'
print(gage)
genv <- gagesenv[gagesenv$RGS_No==gage,]
#Generate FlowScreen time series
gts<- create.ts(predsmelt[predsmelt$ID==gage & !is.na(predsmelt$Flow),])  #Cannot run ts on multiple gages. Need to first subset by gage, then run ts.
gts_sel <- merge(gts,  rufidat_gapyear[!is.na(rufidat_gapyear$gagecount_gap_0.1),c('ID', 'gagecount_gap_0.1')], by='ID', all.x=T)
gname <- paste(genv$RGS_Loc," river at ",genv$RGS_Name,". Selected data: ", unique(gts_sel$gagecount_gap_0.1)," years.",sep="")

#Make raw time series plot
template <-ggplot() +
  geom_point(data=gts_sel, aes(x=Date, y=Flow, color='brown'), size=1.5)+ 
  geom_point(data=rufidat_deleted[rufidat_deleted$ID==gage,], aes(x=Date, y=Flow,color='red'), size=1.5)+
  geom_point(data=rufidat_clean[rufidat_clean$ID==gage,], aes(x=Date, y=Flow,color='lightblue'), size=1.5)+
  geom_point(data=rufidat_clean[rufidat_clean$ID==gage & 
                                  rufidat_clean$hyear %in%  as.data.frame(rufidat_gapsummary)[rufidat_gapsummary$ID==gage & rufidat_gapsummary$gap_per<=0.1,'hyear'],],
             aes(x=Date, y=Flow,color='darkblue'), size=1.5)+
  scale_colour_manual(name='Hydrologic record',
                      values =c('brown'='#bf812d','red'='#e31a1c','lightblue'='#9ecae1','darkblue'='#045a8d'),
                      labels = c('Interpolated','< 10% missing data','> 10% missing data','Deleted')) +
  scale_y_sqrt(expand=c(0,0),limits=c(0,max(gts$Flow)+1))+
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", limits=as.Date(c('1954-01-01','2018-06-01')), expand=c(0,0)) + 
  labs(y='Discharge (m3/s)', title=paste(gage, gname,sep=" - "))+
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.text = element_text(size=8, lineheight=0.1),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.margin=margin(t = 0, unit='cm'),
        legend.key.height=unit(0.75,"line"),
        legend.title = element_text(size=8))
#print(template)
#Get legend as a grob
legendts <- get_legend(template)

##################################Plot 'em #############
plotseries <- function(gage){ #Make a graph of a time series highlighting deleted, interpolated, used and non-used data
  print(gage)
  genv <- gagesenv[gagesenv$RGS_No==gage,]
  #Generate FlowScreen time series
  gts<- create.ts(predsmelt[predsmelt$ID==gage,])  #Cannot run ts on multiple gages. Need to first subset by gage, then run ts.
  gts_sel <- merge(gts, rufidat_ycount, by='ID', all.x=T)
  gts_sel <- merge(gts_sel, rufidat_gapsummary[,c('ID','hyear','gap_per')], by=c('ID','hyear'))
  print(genv$RGC_Loc)
  #Make raw time series plot
  rawsgplot <-ggplot() +
    geom_point(data=gts_sel, aes(x=Date, y=Flow), color='#bf812d', size=0.5)+ 
    geom_point(data=rufidat_deleted[rufidat_deleted$ID==gage,], aes(x=Date, y=Flow),color='#e31a1c', size=0.5)+
    geom_point(data=rufidat_clean[rufidat_clean$ID==gage,], aes(x=Date, y=Flow),color='#9ecae1', size=0.5)+
    geom_point(data=rufidat_clean[rufidat_clean$ID==gage & 
                                    rufidat_clean$hyear %in%  as.data.frame(rufidat_gapsummary)[rufidat_gapsummary$ID==gage & 
                                                                                                  rufidat_gapsummary$gap_per<=0.1,'hyear'],],
               aes(x=Date, y=Flow),color='#045a8d', size=0.5)+
    scale_y_sqrt(expand=c(0,0),limits=c(0,max(gts$Flow)+1))+
    scale_x_date(date_breaks = "2 years", date_labels = "%Y", limits=as.Date(c('1954-01-01','2018-06-01')), expand=c(0,0)) + 
    labs(y='Discharge (m3/s)', 
         title=bquote(paste(.(gage)-.(as.character(genv$RGS_Loc)) ~'River at' ~ .(as.character(genv$RGS_Name)),'.')))+
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          legend.position = 'none',
          text= element_text(size=8,lineheight = 0),
          plot.title= element_text(size=8,margin = margin(t = 0, r = 0, b = 0, l = 0)),
          plot.margin=unit(c(-0,0,-0,-0), "cm"),
          axis.title.x= element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.y=element_text(vjust=-0.5),
          panel.grid.minor= element_blank(),
          panel.border = element_blank(),
          axis.line = element_line()
          )
  p <- ggplot_gtable(ggplot_build(rawsgplot))
  lay= t(c(rep(1,4),2))
  png(file.path(outdir,paste(gage,'raw_sg.png',sep="_")),width = 6.5, height=4.5,units='in',res=600)
  print(grid.arrange(p, legendts, ncol = 11, layout_matrix = lay))
  dev.off()
}

plotflowscreen <- function(gage, div,thrs){ #make graphs of FlowScreen package outputs
  print(gage)
  genv <- gagesenv[gagesenv$RGS_No==gage,]
  #Generate FlowScreen time series
  gts<- create.ts(predsmelt[predsmelt$ID==gage,])  #Cannot run ts on multiple gages. Need to first subset by gage, then run ts.
  gname <- paste(genv$RGS_Loc," river at ",genv$RGS_Name,sep="")
  gts$Flow <- gts$Flow/div

  #Compute and output flowScreen metrics and plots
  try({
    res <- metrics.all(gts,NAthresh=thrs)
    ginfo <- data.frame(StationID=genv$RGS_No, StnName=gname, ProvState='Rufiji Basin',Country='Tanzania',
                        Lat=genv$POINT_Y, Long=genv$POINT_X, Area=genv$WsArea, RHN='RBWB')
    png(file.path(outdir,paste(gage,'screenb.png',sep="_")),width = 20, height=12,units='in',res=600)
    screen.summary(res, type="b", StnInfo=ginfo)
    dev.off()
    png(file.path(outdir,paste(gage,'screenl.png',sep="_")),width = 20, height=12,units='in',res=600)
    screen.summary(res, type="l", StnInfo=ginfo)
    dev.off()
    png(file.path(outdir,paste(gage,'screenh.png',sep="_")),width = 20, height=12,units='in',res=600)
    screen.summary(res, type="h", StnInfo=ginfo)
    dev.off()
  })
}
for (g in unique(predsmelt$ID)) {
  plotseries(g)
  plotflowscreen(g,div=1,thrs=0.5)
}

##################################Report statistics#######################
mean(setDT(rufidat_gapsummary)[,length(hyear),ID]$V1)
min(setDT(rufidat_gapsummary)[,length(hyear),ID]$V1)
max(setDT(rufidat_gapsummary)[,length(hyear),ID]$V1)
"On average, gauges contained 40 years of daily discharge data (min = 2 years, max = 64 years)"

#Compare discharge record for Ruaha and Kilombero
kilombero_avg <- mean(setDT(rufidat_gapsummary)[substr(ID,1,3)=='1KB',length(hyear),ID]$V1)
ruaha_avg <- mean(setDT(rufidat_gapsummary)[substr(ID,1,3)=='1KA',length(hyear),ID]$V1)
ruaha_avg-kilombero_avg
#the discharge records of river gauges in the Kilombero River Basin were 16 years shorter than those in the Great Ruaha River Basin, on average 

mean(rufidat_gapsummary$gap_per)
min(setDT(rufidat_gapsummary)[,mean(gap_per),ID]$V1)
max(setDT(rufidat_gapsummary)[,mean(gap_per),ID]$V1)
"The average percentage of missing data per year across all gauges was 15% (min=1%, max=59%)." 

#Mean number of years with 90% completeness and 37 max gap criterion
mean(ycount_full$V1)

############################## FIGURE 4-2 ##############################
#####################################################################
#Assess representativity of gages regarding environmental variables
#####################################################################
##################################With individual variable plots#####
theme_env <- function () { 
  theme_classic(base_size=18) %+replace% 
    theme(
      panel.background  = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x=element_line(),
      axis.line.y=element_line(),
      legend.position = 'none'
    )
}
theme_envnoy <- function () { 
  theme_classic(base_size=18) %+replace% 
    theme(
      panel.background  = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x=element_line(),
      axis.line.y=element_line(),
      axis.title.y=element_blank(),
      legend.position='none'
    )
}

envplot <- function(selected_gages, plotname) {
  gagesenvrec[gagesenvrec$RGS_No %in% unique(selected_gages$ID),'select'] <- 'Y'
  gagesenvrec[!(gagesenvrec$RGS_No %in% unique(selected_gages$ID)),'select'] <- 'N'
  
  selRGB <- rgb(168,0,0,maxColorValue = 255)
  notselRGB <- rgb(150,150,150, maxColorValue = 255)
  
  envplot_area <- ggplot(rufienv, aes(x=WsArea)) + 
    geom_vline(data=gagesenvrec, aes(xintercept=WsArea, color=select),alpha=0.75, size=0.75) +
    geom_histogram(bins=50, fill='#878787', alpha=0.5) + 
    scale_x_log10(name=expression('Catchment area'~(km^2)),
                  breaks = c(1,10,100,1000,10000,100000),
                  labels = c(1,10,expression(10^2),expression(10^3),expression(10^4),expression(10^5)),
                  expand=c(0,0)) +
    scale_y_continuous(name='Number of river reaches', expand=c(0,0)) + 
    scale_color_manual(values=c(notselRGB,selRGB))+
    theme_env() + 
    labs(subtitle = "(a)")
  #envplot_area
  envplot_elv <- ggplot(rufienv, aes(x=ReaElvAvg)) +
    geom_vline(data=gagesenvrec, aes(xintercept=ReaElvAvg, color=select),alpha=0.75, size=0.75) +
    geom_histogram(bins=50,fill='#fdbf6f', alpha=0.65) + 
    scale_x_continuous(name=expression('River reach average elevation (m)'),
                       expand=c(0,0)) +
    scale_y_continuous(name='Number of river reaches',expand=c(0,0)) +
    scale_color_manual(values=c(notselRGB,selRGB))+
    theme_envnoy()+ 
    labs(subtitle = "(b)")
  #envplot_elv
  envplot_preci <- ggplot(rufienv, aes(x=WsBio12Av)) + 
    geom_vline(data=gagesenvrec, aes(xintercept=WsBio12Av, color=select),alpha=0.75, size=0.75) +
    geom_histogram(bins=50,fill='#1f78b4', alpha=0.4) + 
    scale_x_continuous(name=expression('Catchment mean annual precipitation (mm)'),
                       expand=c(0,0)) +
    scale_y_continuous(name='Number of river reaches',expand=c(0,0)) + 
    theme(axis.title.y=element_blank()) +
    scale_color_manual(values=c(notselRGB,selRGB))+
    theme_envnoy()+ 
    labs(subtitle = "(c)")
  #envplot_preci
  
  
  envplot_watext <- ggplot(rufienv, aes(x=100*CatWatExt)) + 
    geom_vline(data=gagesenvrec, aes(xintercept=100*CatWatExt, color=select),alpha=0.75, size=0.75) +
    geom_histogram(bins=50,fill='#35978f', alpha=0.4) + 
    scale_x_continuous(name=expression('Subcatchment maximum water extent 1984-2015 (% area)'),
                       expand=c(0,0)) +
    scale_y_continuous(trans='sqrt',name='Number of river reaches',expand=c(0,0)) + 
    scale_color_manual(values=c(notselRGB,selRGB))+
    theme_env()+
    theme(axis.title.x=element_text(hjust=0.90))+ 
    labs(subtitle = "(d)")
  #envplot_watext
  rufienv$CatGeolMaj <- factor(rufienv$CatGeolMaj, levels = unique(rufienv$CatGeolMaj[order(as.numeric(as.character(rufienv$CatGeolMaj)))]))
  envplot_geol <-  ggplot(rufienv, aes(x=CatGeolMaj)) + 
    geom_bar(fill='#8c510a', alpha=0.4, stat='count') +
    geom_bar(data=gagesenvrec, aes(x=as.factor(CatGeolMaj), color=select),fill='white',alpha=0.4, size=0.75) +
    scale_x_discrete(name=expression('Main subcatchment lithology (GMIS number)')) +
    scale_y_continuous(trans='log10',name='Number of river reaches') + 
    scale_color_manual(values=c(notselRGB,selRGB))+
    theme_envnoy() +
    theme(axis.text.x = element_text(angle = 45, hjust=1))+ 
    labs(subtitle = "(e)")
  #envplot_geol
  envplot_pop <- ggplot(rufienv, aes(x=CatPopDen)) + 
    geom_vline(data=gagesenvrec, aes(xintercept=CatPopDen, color=select),alpha=0.75, size=0.75) +
    geom_histogram(bins=50,fill='#4d9221', alpha=0.4) + 
    scale_x_log10(name=expression('Subcatchment population density'~(persons~km^-2)),
                  breaks=c(1,5,10,50,100,1000),
                  expand=c(0,0)) +
    scale_y_continuous(name='Number of river reaches',expand=c(0,0)) + 
    scale_color_manual(values=c(notselRGB,selRGB))+
    theme_envnoy()+ 
    labs(subtitle = "(f)")
  envplot_pop
  
  
  envplot_forlos <- ggplot(rufienv, aes(x=100*WsFLosSum_1)) + 
    geom_vline(data=gagesenvrec, aes(xintercept=100*WsFLosSum_1, color=select),alpha=0.75, size=0.75) +
    geom_histogram(bins=50,fill='#bf812d', alpha=0.6) + 
    scale_x_sqrt(name=expression('Catchment forest cover loss 2000-2016 (% area)'),
                 breaks=c(0,1,5,10,25,50,75,100),
                 expand=c(0,0)) +
    scale_y_continuous(name='Number of river reaches',expand=c(0,0)) + 
    scale_color_manual(values=c(notselRGB,selRGB))+
    theme_env() +
    theme(axis.title.x=element_text(hjust=1))+ 
    labs(subtitle = "(g)")
  #envplot_forlos
  envplot_urb <- ggplot(rufienv, aes(x=100*LCSum_89)) + 
    geom_vline(data=gagesenvrec, aes(xintercept=100*LCSum_89, color=select),alpha=0.75, size=0.75) +
    geom_histogram(bins=50,fill='#cab2d6', alpha=0.4) + 
    scale_x_continuous(trans='log10',
                       name='Catchment urban cover (% area)',
                       breaks = c(0.001,0.01,0.1,1,10),
                       labels= c(0.001,0.01,0.1,1,10),
                       expand=c(0,0)) +
    scale_y_continuous(name='Number of river reaches',expand=c(0,0)) + 
    scale_color_manual(values=c(notselRGB,selRGB))+
    theme_envnoy()+ 
    labs(subtitle = "(h)")
  #envplot_urb
  envplot_resind <- ggplot(rufienv, aes(x=100*WsResInd)) + 
    geom_vline(data=gagesenvrec, aes(xintercept=100*WsResInd, color=select),alpha=0.75, size=0.75) +
    geom_histogram(bins=50,fill='#de77ae', alpha=0.4) + 
    scale_x_continuous(name=expression('Catchment draining through nearest upstream reservoir (%)'),
                       expand=c(0,0)) +
    scale_y_log10(name='Number of river reaches',expand=c(0,0),
                  breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_color_manual(values=c(notselRGB,selRGB))+
    theme_envnoy() +
    theme(axis.title.x=element_text(hjust=0.90))+ 
    labs(subtitle = "(i)")
  #envplot_resind
  #plot_grid(envplot_area, envplot_elv, envplot_preci, envplot_watext, envplot_geol, envplot_pop, envplot_forlos, envplot_urb, envplot_resind, align = "v", nrow = 3)
  
  png(file.path(outdir,plotname),width=20, height=12,units='in',res=600)
  print(grid.arrange(envplot_area, envplot_elv, envplot_preci, envplot_watext,
                     envplot_geol, envplot_pop, envplot_forlos, envplot_urb, 
                     envplot_resind, nrow = 3))
  dev.off()
}
envplot(rufidat_select_o15y, 'gage_envo15y_20200628.png')

##################################In multidimensional environment (never fully implemented)####################
#Make subset of data
# colnames(rufienv)
# outcols <- c(1:4,6,7,9,45:70,92:112, which(colnames(rufienv) %in% c('CatFlowAcc','CatElvMin','CatDen','CatDamDen','CatFlowAcc','CatLCMaj',
#                                                                     'WsPAPer','WsDamDen','WsGeolMaj','WsLCMaj','ReaElvMin',
#                                                                     'ReaElvMax','SUM_LENGTH_GEO','Shape_Leng') |
#                                              !is.na(str_match(colnames(rufienv),'DirSum*'))))
# rufienvsub <- rufienv[,-outcols]
# rufienvsub$ReaDirMaj <- as.factor(rufienvsub$ReaDirMaj)
# colnames(rufienvsub)
# envsubcat <- rufienvsub[,c(1:55,157:161)]
# envsubws <-rufienvsub[,-c(1:55)]
# 
# #Data transformation
# #Transform catchments
# str(envsubcat)
# colnames(envsubcat)
# factcol <- c(1,2,54,55,57)
# #hist.plots(envsubcat[,-factcol]) #Inspect data
# logcols <- c('CatPopDen','ReaSloAvg')
# envsubcat[,logcols] <- data.trans(data.frame(envsubcat[,logcols]), method = 'log', plot = F)
# asincols <- c('CatFLosSum', paste('LCSum',c(1,2,3,4,5,6,7,8,10),sep='_'),'CatWatExt','CatResInd','CatLakInd')
# envsubcat[,asincols] <- data.trans(envsubcat[,asincols], method = 'asin', plot = F)
# sqrtcols <- c('CatAIAvg', 'CatBio14Av','CatBio17Av','CatBio19Av','CatElvMax', 'CatElvAvg','CatSloAvg','CatSloStd','CatLen_1','CatPAPer',
#               'CatRoadDen','CatWatcha','CatMineDen','CatWatOcc','ReaPAPer','ReaElvAvg')
# envsubcat[,sqrtcols] <- data.trans(envsubcat[,sqrtcols], method = 'power',exp=.5, plot = F)
# #Transform watersheds
# #hist.plots(envsubws) #Inspect data
# 
# #Then standardize to mean of 0 and unit variance
# envsubcat_std <- envsubcat[,-factcol]
# envsubcat_std <- cbind(data.stand(envsubcat_std, method = "standardize", margin = "column", plot = F),
#                        envsubcat[,factcol])
# #Transfer first column to row name
# rownames(envsubcat_std) <- envsubcat_std$GridID
# envsubcat_std <- envsubcat_std[,-which(colnames(envsubcat_std)=='GridID')]
# #Establish variable weights
# #envar <- colnames(rufienvsub)
# #weight_gow <- c(0.5, 0.5, 1,1,1, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.333, 0.333, 0.333)
# 
# #Computer Gower's dissimilarity on full dataset
# set.seed(1)
# envsubcat_samp<-stratified(envsubcat_std, "ReaOrd", 0.1, select=NULL)
# rufi_gowdis <- gowdis(envsubcat_samp, w=rep(1,ncol(envsubcat_samp)), asym.bin = NULL)
# attributes(rufi_gowdis)
# #Convert gowdis output to matrix
# rufi_gowdis<- as.matrix(rufi_gowdis)
# #Run PCoA without a constant added using the stats package
# rufi_pcoa <- cmdscale(rufi_gowdis, k = 3, eig = T, add = F)
# #Check eigenvalues
# rufi_pcoa$eig
# rufi_pcoa$GOF
# 
# #rufi_pcoa_eucldist <- vegdist(rufi_pcoa$points, method = "euclidean")
# #Shepard <- Shepard_diy(rufi_gowdis, rufi_pcoa_eucldist)
# #Calculate PC loadings using correlation analysis - misses 124 species, so not necessarily that representative
# vec_tr <- envfit(fish_pcoa,  Fish_traits_select, perm = 1000, na.rm = T)
# arrows <- as.data.frame(vec_tr$vectors$arrows)
# factors <- as.data.frame(vec_tr$vectors$factors)
# 
# #Visualize scores
# pcoa_scores <- as.data.frame(fish_pcoa$points)
# pcoa_scores$spe <- rownames(pcoa_scores)
# colnames(pcoa_scores) <- c("PC1", "PC2")
# ggplot(pcoa_scores, aes(x = PC1, y = PC2, label = rownames(pcoa_scores))) + geom_label() +
#   geom_segment(data = arrows, aes(x = rep(0, 13), y = rep(0, 13), xend = Dim1 , yend = Dim2, label = rownames(arrows))) + 
#   geom_text(aes(x = Dim1, y = Dim2, label = rownames(arrows)), data = arrows, color = "red") +
#   
#   
#   ggplot(arrows) + geom_segment(aes(x = rep(0, 13), y = rep(0, 13), xend = Dim1 , yend = Dim2))
# 
# ###########################################
# #Extra: 
# #rufidat_overlapplot <- ldply(seq(5,50,1), function(y) {
# #  rufidat_dtgap[ID %in% as.character(rufidat_gapyear[rufidat_gapyear$gagecount_gap_0.1>y,'ID']),list(length(unique(ID)), minyr=y), .(Date)]
# #})
