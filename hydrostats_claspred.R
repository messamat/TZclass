#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/30/2018
#Date last updated: 05/25/2018

#Purpose: compute hydrological metrics, classify stream gauges, and predict class membership for Rufiji network

library(foreign)
library(plyr)
library(dplyr)
library(hydrostats)
library(data.table)
#devtools::install_github("messamat/EflowStats") #Intro Vignette: https://cdn.rawgit.com/USGS-R/EflowStats/9507f714/inst/doc/intro.html
#Corrected a glitch in package, need to re-change package download to USGS-R/EflowStats
library(EflowStats)
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
library(zoo)

rootdir="F:/Tanzania/Tanzania" #####UPDATE THIS TO MATCH YOUR ROOT PROJECT FOLDER #######

source(file.path(rootdir,"bin/outside_src/Biostats.R"))
source(file.path(rootdir,"bin/outside_src/Flowscreen.hyear.internal.R"))

#Set folder structure
setwd(file.path(rootdir,"results")) 
datadir = file.path(getwd(),'rufiji_hydrodatafilter')
origdatadir = file.path(rootdir,"data") 
outdir=file.path(getwd(),'rufiji_hydrodatastats')
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}

#Import data
rufidat_clean <- read.csv(file.path('rufiji_hydrodatainspect','rufidat_clean.csv'), colClasses=c('factor','Date','numeric','character','character')) #Cleaned hydro data
rufidat_impute <- read.csv(file.path('rufiji_hydrodataimpute', 'rufidat_interp.csv'), colClasses=c('Date',rep(c('numeric','numeric','character'),39))) #Cleaned + interpolated hydro data
colnames(rufidat_impute)[seq(2,ncol(rufidat_impute),3)] <- substr(colnames(rufidat_impute),2,10)[seq(2,ncol(rufidat_impute),3)]
rufidat_gapsummary <- read.csv(file.path(datadir, 'rufidat_gapsummary.csv')) #Assessment of data availability 
#rufidat_post1991<-read.csv(file.path(datadir, 'gageselect_post1991comp90.csv')) #Selection of gauges that have at least 10 years  since 1991
rufidat_o15y<-read.csv(file.path(datadir, 'gageselect_o15comp90.csv')) #Number of years with >90% data for each gauge

gagesenv <- read.csv(file.path(getwd(),'gages_netjoinclean.csv')) #Import gauges' environmental data
gagesenvrec <- merge(gagesenv, unique(rufidat_clean[,c('ID','SYM')]), by.x='RGS_No', by.y='ID', all.x=F)
rufienv <- read.csv(file.path(getwd(),'streamnet118_rufiji_finaltabclean.csv')) #Import river network environmental data

#Format data for computing hydro metrics
predsmelt <-melt(setDT(as.data.frame(rufidat_impute)[,c(1,which(colnames(rufidat_impute) %like% "^1K*"))]),id.vars = 'Date',value.name='Flow',variable.name='ID')
predsmelt <- predsmelt[,c(2,1,3)]
predsmelt$year <- as.numeric(format(predsmelt$Date, "%Y"))
predsmelt$month <- as.numeric(format(predsmelt$Date, "%m"))
predsmelt<-hyear.internal(predsmelt,hyrstart=10) #hdoy doesn't work
#Compute hydrologic day
predsmelt$doy <- as.numeric(format(predsmelt$Date,"%j"))
predsmelt[,hdoy:=ifelse(month>=10,
                        doy-as.numeric(format(as.Date(paste(year,'-10-01',sep="")),"%j")),
                        doy+as.numeric(format(as.Date(paste(year-1,'-12-31',sep="")),"%j"))-as.numeric(format(as.Date(paste(year-1,'-10-01',sep="")),"%j")))] 
predsmelt <- merge(predsmelt, rufidat_gapsummary, by=c('ID','hyear'),all.x=T)
predsmelt <- merge(predsmelt, rufidat_post1991, by='ID',all.x=T)
predsmelt <- merge(predsmelt, rufidat_o15y, by='ID',all.x=T)

############################################### Compute hydrologic metrics##########################################
allHITcomp <- function(dfhydro, dfenv, gageID, templateID='1KA9',hstats="all", floodquantile=0.95) {
  ####Get template
  dailyQClean <- validate_data(dfhydro[dfhydro$ID==templateID,c("Date", "Flow")], yearType="water")
  #Calculate all hit stats
  HITall_template <- calc_allHIT(dailyQClean, yearType="water", stats=hstats, digits=10, pref="mean",
                                drainArea=dfenv[dfenv$RGS_No==templateID,'WsArea'], floodThreshold = quantile(dailyQClean$discharge, floodquantile))
  colnames(HITall_template)[2] <- templateID
  HITall <- data.frame(indice=HITall_template$indice) 
  ####Compute metrics for all gages
  for (gage in unique(dfhydro[,gageID])) {
    print(gage)
    try({
      #Check data for completeness
      dailyQClean <- validate_data(dfhydro[dfhydro$ID==gage,c("Date", "Flow")], yearType="water")
      #Calculate all hit stats
      calc_allHITout <- calc_allHIT(dailyQClean, yearType="water", stats=hstats, digits=10, pref="mean",
                                    drainArea=dfenv[dfenv$RGS_No==gage,'WsArea'], floodThreshold = quantile(dailyQClean$discharge, floodquantile))
      colnames(calc_allHITout)[2] <- gage
      HITall <- merge(HITall, calc_allHITout, by='indice')
    })
  }
  HITall_formatmelt <-melt(setDT(HITall), id.vars = "indice",variable.name = gageID) 
  HITall_formatmelt[is.infinite(HITall_formatmelt$value),'value'] <- NA
  HITall_formatmelt[is.nan(HITall_formatmelt$value),'value'] <- NA
  return(HITall_formatmelt)
}

#Select subset of gauges: include 1KB32 even if only 14 years of data + Kisigo stations
rufidat_select_o15y <- predsmelt[predsmelt$gap_per<=0.1 & predsmelt$max_gap < 37 & predsmelt$hyear<2017 & 
                                   (predsmelt$ycount_o15>=15 | predsmelt$ID=='1KB32') | 
                                   (predsmelt$ID=='1KA41' & predsmelt$max_gap < 273 & predsmelt$gap_per<0.75 & predsmelt$hyear<1996) |
                                   (predsmelt$ID=='1KA42A' & predsmelt$max_gap < 273 & predsmelt$gap_per<0.75 & predsmelt$hyear>1958 & predsmelt$hyear<2017),]
#Compute HITs
HITo15y <- allHITcomp(as.data.frame(rufidat_select_o15y), gagesenv, 'ID')
write.csv(HITo15y, file.path(outdir, 'HITo15y.csv'),row.names=F)

########################Inspect NA values in HITs
#1KA41
median(as.data.frame(rufidat_select_o15y)[rufidat_select_o15y$ID=='1KA41','Flow']) #Any metric that relies on dividing by median or some monthly flow is NA
#dl19 is NA for many stations that have no 0-flow days
#1KA2A
as.data.frame(rufidat_select_o15y)[rufidat_select_o15y$ID=='1KA2A' & 
                                     rufidat_select_o15y$Flow> 7*median(as.data.frame(rufidat_select_o15y)[rufidat_select_o15y$ID=='1KA2A','Flow']),'Flow']  #No discharge >7x median flow (mh23, mh26)
#1KA42: any metric that relies on dividing by some monthly flow is NA
#1KA50B
ggplot(as.data.frame(rufidat_select_o15y)[rufidat_select_o15y$ID=='1KA50B',], aes(x=Date,y=Flow))+geom_point() #ma31 and ma32 required dividing by monthly flows, yields NA
#1KA59:  any metric that relies on dividing by some monthly flow is NA
#All others with mh22, mh23, mh25, and mh26 don't have flows exceeding a given factor of median flow

############################################### Box plot and table of metrics##############################
HITboxplot <- function(HITdf, plotname) {
  HITallbox<- HITdf
  HITallbox$group1 <- as.factor(substr(HITdf$indice,1,1)) #Subset metric name into main category
  HITallbox$group2 <- as.factor(substr(HITdf$indice,2,2)) #and l, h, and a
  HITallbox$indice_sub <- substr(HITdf$indice,3,5) #and metric number
  HITallbox$indice_sub <- factor(HITallbox$indice_sub, levels = unique(HITallbox$indice_sub[order(as.numeric(as.character(HITallbox$indice_sub)))]))
  HITallbox$group1 <- factor(HITallbox$group1, levels = c('m','f','d','t','r'), 
                             labels = c("Magnitude-m", "Frequency-f", "Duration-d",'Timing-t',"Rate of change-r")) 
  HITallbox$group2 <- factor(HITallbox$group2, levels = c('h','a','l'), labels=c("High flow-h",'Average flow-a',"Low flow-l"))
  # lab <- ddply(HITallbox, .(indice),
  #              labels=median(value)-1.58*(quantile(value,0.75,na.rm=T)-quantile(value,0.25,na.rm=T))/sqrt(length(value)),
  #              labels2=min(value))
  #HITallbox <-  merge(HITallbox,lab, by=c('ID','indice'))
  HITallboxplot <-ggplot(HITallbox, aes(x=indice_sub, y=value, color=group1)) + 
    scale_y_log10(name='Metric value') +
    geom_boxplot() +
    facet_grid(group2~group1, scales = "free", space="free_x") + 
    scale_x_discrete(name='Metric number (Appendix 1)')+
    theme_classic() +
    theme(axis.title = element_text(size=16),
          axis.text.y = element_text(size=16),
          strip.text = element_text(size = 15.5),
          legend.position='none')
  png(file.path(outdir,plotname),width=20, height=12,units='in',res=300)
  print(HITallboxplot)
  dev.off()
}
#HITboxplot(HITpost1991,'HITallboxplotpost1991.png')
HITboxplot(HITo15y, 'HITallboxploto15y.png')

# Make table reporting mean (SD) for each metric for all gauges (appendix)
HITdf_cast <- dcast(HITo15y, ID ~ indice)

############################################### Compute redundancy in hydrologic metrics and subset them ################
#Compute redundancy in metrics
HITo15y_filled <- replace.missing(HITo15y, method='mean')
HITdf_cast <- as.data.frame(dcast(HITo15y_filled, ID ~ indice))
rownames(HITdf_cast) <- as.character(HITdf_cast$ID)
HITcor <- cor(HITdf_cast[,-1],use='pairwise.complete.obs')
HITcor[HITcor<0.9] <- NA
write.csv(HITcor, file.path(outdir, 'HITcor.csv'), row.names=T)

#Subset 1
sub1 <- c('dh1','dh2','dh4','dh6','dh9','dh10','dh11','dh12','dh13','dl1','dl3','dl5','dl6','dl7','dl8',
          'dl11','dl12','ma2','ma7','ma10','ma38','ma40','mh22','mh25','mh26','ml15','ml16','ml19','ml21')
HITo15ysub1 <- droplevels(HITo15y[!(HITo15y$indice %in% sub1),]) 
HITcorsub1 <- cor(dcast(HITo15ysub1, ID ~ indice)[,-1],use='pairwise.complete.obs')
HITcorsub1[HITcorsub1<0.99 & HITcorsub1>-0.99] <- NA
write.csv(HITcorsub1, file.path(outdir, 'HITcorsub.csv'), row.names=T)

#Subset 2
sub2 <- c('dh1','dh2','dh4','dh6','dh9','dh10','dh11','dh12','dh13','dl1','dl2','dl4','dl5','dl6','dl7','dl8','dl9',
          'dl11','dl12','dl20','ma2','ma7','ma10','ma36','ma37','ma38','ma40','mh15','mh22','mh25','mh26','mh27',
          'ml15','ml16','ml19','ml21','fh8')
HITo15ysub2 <- droplevels(HITo15y[!(HITo15y$indice %in% sub2),]) 
HITcorsub2 <- cor(dcast(HITo15ysub2, ID ~ indice)[,-1],use='pairwise.complete.obs')
HITcorsub2[HITcorsub2<0.90 & HITcorsub2>-0.90] <- NA
write.csv(HITcorsub2, file.path(outdir, 'HITcorsub.csv'), row.names=T)

#subset 3
sub3 <- c('dh1','dh2','dh3','dh4','dh5','dh6','dh8','dh9','dh10','dh11','dh12','dh13','dh20','dl1','dl2','dl3','dl4','dl5',
          'dl6','dl7','dl8','dl9','dl10','dl11','dl12','dl13','dl14','dl20','ma2','ma7','ma10','ma11','ma36','ma37','ma38',
          'ma40','ma42','ma45','mh13','mh14','mh15','mh17','mh22','mh25','mh26','mh27','ml13','ml15','ml16','ml19','ml21','fh8')
length(sub3)
HITo15ysub3 <- droplevels(HITo15y[!(HITo15y$indice %in% sub3),]) 

#Ordinate metrics
HITpca<- prcomp(HITdf_cast[,-1], scale=T)
summary(HITpca)
ordiplot(HITpca, choices=c(1,2), type='text', display='sites')
arrows(0,0,HITpca$rotation[,1]*100,HITpca$rotation[,2]*100, col='grey')
text(HITpca$rotation[,1]*100,HITpca$rotation[,2]*100, row.names(HITpca$rotation), col='red')
text(HITpca$rotation[,1][row.names(HITpca$rotation) %in% sub3]*100,HITpca$rotation[,2][row.names(HITpca$rotation) %in% sub3]*100, row.names(HITpca$rotation)[row.names(HITpca$rotation) %in% sub3], col='orange')

ordiplot(HITpca, choices=c(3,4), type='text', display='sites')
arrows(0,0,HITpca$rotation[,3]*100,HITpca$rotation[,4]*100, col='grey')
text(HITpca$rotation[,3]*100,HITpca$rotation[,4]*100, row.names(HITpca$rotation), col='red')
text(HITpca$rotation[,3][row.names(HITpca$rotation) %in% sub3]*100,HITpca$rotation[,4][row.names(HITpca$rotation) %in% sub3]*100, row.names(HITpca$rotation)[row.names(HITpca$rotation) %in% sub3], col='orange')

############################################### Format hydrologic metrics to use in classification and compute Gower's distance############
HITdist <- function(HITdf, logmetrics) { 
  if (logmetrics==TRUE){
    HITdf$Value <- log(HITdf$Value+1) #log-transform metric
  }
  HITdf_format <- dcast(HITdf, ID ~ indice)
  HITdf_format <- merge(HITdf_format, gagesenvrec[,c('RGS_No','WsArea')], by.x='ID', by.y='RGS_No')
  dimindices <- c('ma1','ma2',paste('ma',seq(12,23),sep=''),paste('ml',seq(1,12),sep=''),paste('mh',seq(1,12),sep=''), 
                  paste('dl',seq(1,5),sep=''),paste('dh',seq(1,5),sep=''),'ra1','ra3','ra6','ra7') #List of dimensional indices 
  dimindices <- dimindices[dimindices %in% colnames(HITdf_format)] #Make sure they are all in the dataset
  HITdf_format <- as.data.frame(setDT(HITdf_format)[,(dimindices) := lapply(.SD, function(x) round(x/WsArea, digits=10)), .SDcols=dimindices]) #Standardize dimensional indices by drainage area
  row.names(HITdf_format) <- HITdf_format$ID 
  HITdf_format <- HITdf_format[,-which(colnames(HITdf_format) %in% c('ID','WsArea'))] #Get rid of non-indices columns
  HITdf_stand <- data.stand(HITdf_format[,2:(ncol(HITdf_format))],method='standardize',margin='column',plot=F) #z-standardize columnwise 
  gauge_gow<- gowdis(HITdf_stand, w=rep(1,ncol(HITdf_stand)), asym.bin = NULL) #Compute Gower's distance so that missing values will not be taken in account
  return(gauge_gow)
}

######################################### CLASSIFICATION BASED ON ENTIRE PERIOD > 15 YEARS OF DATA ################################
################################################ Classify based on all indices and diagnostic ############################################
gaugegow_o15y <- HITdist(HITo15y, logmetrics=TRUE) #Format hydro metrics and compute Gower's distance matrix
gaugecla_ward <-hclust(gaugegow_o15y, method='ward.D') #Classify using hierarchical agglomerative using Ward's minimum variance method
gaugecla_ward2 <-hclust(gaugegow_o15y, method='ward.D2') #Classify using hierarchical agglomerative using Ward's minimum variance method
gaugecla_UPGMA <-hclust(gaugegow_o15y, method='average') #Classify using  UPGMA

#Classification diagnostics
cluster_diagnostic <- function(clusterres, clusname, gowdis) {
  #hclus.table(clusterres)
  coef.hclust(clusterres) #Compute agglomerative coefficient
  #cor(gowdis, cophenetic(clusterres)) #Compute cophenetic coefficient
  #Plot cophenetic relationship 
  png(file.path(outdir, paste('o15y_r6r_cophe',clusname,'.png',sep="")), width=8, height=8, units='in',res=300)
  hclus.cophenetic(gowdis, clusterres) 
  dev.off()
  #Scree plot
  png(file.path(outdir, paste('o15y_r6r_scree',clusname,'.png',sep="")), width=8, height=8, units='in',res=300)
  hclus.scree(clusterres) 
  dev.off()
  #Plot dendogram
  png(file.path(outdir, paste('o15y_r6r_dendogram',clusname,'.png',sep="")), width=8, height=8, units='in',res=300)
  plot(clusterres, main=paste(clusname, "gauge cluster dendogram",sep=" "), xlab='Gauge ID', ylab="Gower's distance", hang=-1)   
  rect.hclust(clusterres, k=4) #Draw rectangle around k classes
  rect.hclust(clusterres, k=5) 
  rect.hclust(clusterres, k=6) 
  rect.hclust(clusterres, k=7) 
  rect.hclust(clusterres, k=8) 
  dev.off()
}
cluster_diagnostic(gaugecla_ward, "Ward's D", gaugegow_o15y)
cluster_diagnostic(gaugecla_ward2, "Ward's D2", gaugegow_o15y)
cluster_diagnostic(gaugecla_UPGMA, "UPGMA", gaugegow_o15y)
#UPGMA leads to too much chaining, and Ward's D has higher cophenetic correlation and reaches an elbow after 7 (rather than 8 classes for D2)

#Test significance of classes
clus.stab <- pvclust(t(HITdf_cast[,-1]), method.hclust='ward.D', method.dist='cor',use.cor="pairwise.complete.obs", nboot=4999)
plot(clus.stab)
pvrect(clus.stab, alpha=0.90)

#TO DO:
# - Compare classifications based on adjusted Rand Index (ARI)

#Get gauge classes
classr7 <-cutree(gaugecla_ward, k=7, order_clusters_as_data = FALSE)
classr7_df <- data.frame(ID=names(classr7), gclass=classr7) 
outdirclass <- file.path(outdir,'classo15y_ward_raw')
if (dir.exists(outdirclass )) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdirclass))
  dir.create(outdirclass )
}
write.csv(classr7_df, file.path(outdirclass,'class_rawgow_ward_7.csv'), row.names=F)
classcol<- c("#176c93","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#7a5614") #7 classes with darker color (base blue-green from Colorbrewer2 not distinguishable on printed report and ppt)

#Make good looking dendogram
prettydend <- function(gaugecla, outdir, imgname, colors=classcol, kclass=7) {
  gaugecla_ward_name <- gaugecla
  gaugecla_ward_name$labels <- with(gagesenvrec[gagesenvrec$RGS_No %in% gaugecla_ward_name$labels,], 
                                    paste(RGS_No,"-",RGS_Loc," River at ", RGS_Name,sep=""))
  dendname <- as.dendrogram(gaugecla_ward_name)
  png(file.path(outdirclass,imgname),width = 8, height=8,units='in',res=300)
  par(mar=c(3,3,0,17)) #bottom left top right
  dendname %>% set("branches_lwd", 2.5) %>% 
    color_branches(k=kclass, col=colors, groupLabels=T) %>% 
    #color_branches(clusters=as.numeric(temp_col), col=levels(temp_col), groupLabels=as.character(as.numeric(temp_col))) %>% 
    color_labels(k=kclass, col=colors) %>%
    plot(horiz=TRUE,xlab="Gower's distance", ylab="Gauge ID - River at Location",mgp=c(1.5,0.5,0))
  dev.off()
}
prettydend(gaugecla_ward, outdir=outdirclass,imgname='7class_dendrogram.png')

################################################ Classify based on subsetted indices and diagnostic ############################
#Subset 1
gaugegow_o15ysub1 <- HITdist(HITo15ysub1, logmetrics=TRUE) #Format hydro metrics and compute Gower's distance matrix
gaugecla_wardsub1 <-hclust(gaugegow_o15ysub1, method='ward.D') #Classify using hierarchical agglomerative using Ward's minimum variance method
gaugecla_ward2sub1 <-hclust(gaugegow_o15ysub1, method='ward.D2') #Classify using hierarchical agglomerative using Ward's minimum variance method
gaugecla_UPGMAsub1 <-hclust(gaugegow_o15ysub1, method='average') #Classify using UPGMA
cluster_diagnostic(gaugecla_wardsub1, "Ward's D sub1", gaugegow_o15ysub1)
cluster_diagnostic(gaugecla_ward2sub1, "Ward's D2 sub1", gaugegow_o15ysub1)
cluster_diagnostic(gaugecla_UPGMAsub1, "UPGMA sub1", gaugegow_o15ysub1)

#Subset 2
gaugegow_o15ysub2 <- HITdist(HITo15ysub2, logmetrics=TRUE) 
gaugecla_wardsub2 <-hclust(gaugegow_o15ysub2, method='ward.D') 
cluster_diagnostic(gaugecla_wardsub2, "Ward's D sub2", gaugegow_o15ysub2)

#Subset 3
gaugegow_o15ysub3 <- HITdist(HITo15ysub3, logmetrics=TRUE) 
gaugecla_wardsub3 <-hclust(gaugegow_o15ysub3, method='ward.D') 
gaugecla_ward2sub3 <-hclust(gaugegow_o15ysub3, method='ward.D2')
gaugecla_UPGMAsub3 <-hclust(gaugegow_o15ysub3, method='average')
cluster_diagnostic(gaugecla_wardsub3, "Ward's D sub3", gaugegow_o15ysub3)
cluster_diagnostic(gaugecla_ward2sub3, "Ward's D2 sub3", gaugegow_o15ysub3)
cluster_diagnostic(gaugecla_UPGMAsub3, "UPGMA sub3", gaugegow_o15ysub3)
#The only class difference is that 1KA11 changes group to be more spatially contiguous with headwater G. Ruaha. Between Ward's D and D2,
#Difference in relateness of major groups.

#Test significance of classes (#Doesn't really work I think due to small sample size for several classes)
HITo15ysub3_filled <- replace.missing(HITo15ysub3, method='mean')
HITdfsub3_cast <- as.data.frame(dcast(HITo15ysub3_filled, ID ~ indice))
rownames(HITdfsub3_cast) <- as.character(HITdfsub3_cast$ID)
clus.stab <- pvclust(t(HITdfsub3_cast[,-1]), method.hclust='ward.D', method.dist='euclidean',use.cor="pairwise.complete.obs", nboot=4999)
plot(clus.stab)
pvrect(clus.stab, alpha=0.90)
clus.stab <- pvclust(t(HITdfsub3_cast[,-1]), method.hclust='ward.D2', method.dist='cor',use.cor="pairwise.complete.obs", nboot=4999)
plot(clus.stab)
pvrect(clus.stab, alpha=0.90)

#Output and make dendogram
classr7sub3 <-cutree(gaugecla_ward2sub3, k=7, order_clusters_as_data = FALSE)
classr7sub3_df <- data.frame(ID=names(classr7sub3), gclass=classr7sub3) 
outdirclass <- file.path(outdir,'classo15y_ward2_rawsub3')
if (dir.exists(outdirclass )) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdirclass))
  dir.create(outdirclass )
}
write.csv(classr7sub3_df, file.path(outdirclass,'class_rawgow_ward2_7_sub3.csv'), row.names=F)
classcol<- c("#176c93","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#7a5614") #7 classes with darker color (base blue-green from Colorbrewer2 not distinguishable on printed report and ppt)
prettydend(gaugecla_ward2sub3, outdir=outdirclass,imgname='7class_dendrogram.png')

######################################################## Class hydrograph plots ################
rufidat_select_classr7sub3 <- merge(rufidat_select_o15y, classr7sub3_df, by="ID")
write.csv(rufidat_select_classr7sub3, file.path(outdir, 'classo15y_ward_raw/rufidat_select_classr7sub3.csv'), row.names=F)
setDT(rufidat_select_classr7sub3)[,yrmean:=mean(Flow),.(ID,hyear)] #Compute average daily flow for each station and year
#Compute statistics on long-term daily flow (average, min, max, Q10, Q25, Q75, Q90) across all stations and years for each class
classflowstats <- setDT(rufidat_select_classr7sub3)[,list(classmeanfull=mean(Flow, na.rm=T), classmean= mean(Flow/yrmean,na.rm=T),classQ75= quantile(Flow/yrmean, .75,na.rm=T),
                                                      classQ25=quantile(Flow/yrmean, .25,na.rm=T),classQ90=quantile(Flow/yrmean, .90,na.rm=T),
                                                      classQ10=quantile(Flow/yrmean, .10,na.rm=T),classmax=max(Flow/yrmean,na.rm=T),
                                                      classmin=min(Flow/yrmean,na.rm=T),classsd=sd(Flow/yrmean,na.rm=T), 
                                                      cal_hdoy=format(as.Date(hdoy, origin='2015-10-01'), "%Y-%m-%d")),
                                                .(gclass,hdoy)] 

#Superimposed non-standardized average yearly hydrograph for each class
classhydro_allfull <- ggplot(as.data.frame(classflowstats), aes(x=as.Date(cal_hdoy), y=classmeanfull, color=factor(gclass))) + 
  geom_line(size=1, alpha=0.8) + 
  scale_color_manual(name='Hydrologic class',values=classcol) +
  scale_y_continuous(name='Daily mean discharge',expand=c(0,0),limits=c(0,NA)) + 
  scale_x_date(name='Date',date_breaks = "1 month", date_labels = "%b", expand=c(0,0)) + 
  theme_classic() + 
  theme(legend.position=c(0.8,0.8),
        text=element_text(size=18))
png(file.path(outdirclass,'7class_hydrographfull.png'),width = 16, height=9,units='in',res=300)
print(classhydro_allfull)
dev.off()

#Superimposed standardized average yearly hydrograph for each class
classhydro_all <- ggplot(as.data.frame(classflowstats), aes(x=as.Date(cal_hdoy), y=classmean, color=factor(gclass))) + 
  geom_line(size=1, alpha=0.8) + 
  scale_color_manual(name='Hydrologic class',values=classcol) +
  scale_fill_manual(name='Hydrologic class',values=classcol) +
  scale_y_continuous(name='Daily mean discharge/Mean daily discharge',expand=c(0,0),limits=c(0,NA)) + 
  scale_x_date(name='Date',date_breaks = "1 month", date_labels = "%b", expand=c(0,0)) + 
  theme_classic() + 
  theme(legend.position='none',
        text=element_text(size=18))

#Facetted standardized average yearly hydrograph for each class + Q90-Q10 ribbon
classhydro_facet <-ggplot(as.data.frame(classflowstats), aes(x=as.Date(cal_hdoy))) + 
  #geom_ribbon(aes(ymin=ifelse(classmean-2*classsd>=0,classmean-2*classsd,0), ymax=classmean+2*classsd,
  #                fill=factor(gclass)),alpha=0.3) +
  geom_ribbon(aes(ymin=classQ10, ymax=classQ90,
                  fill=factor(gclass)),alpha=0.3) +
  geom_line(aes(y=classmean, color=factor(gclass)),size=1.2) + 
  facet_grid(gclass~.,scale='free_y') +
  scale_color_manual(name='Hydrologic class',values=classcol) +
  scale_fill_manual(name='Hydrologic class',values=classcol) +
  scale_y_continuous(name='Daily mean discharge/Mean daily discharge',expand=c(0,0),limits=c(0,NA)) + 
  scale_x_date(name='Date',date_breaks = "1 month", date_labels = "%b", expand=c(0,0)) + 
  annotate("segment", x=as.Date('2015-10-01'), xend=as.Date('2016-09-30'), y=0, yend=0)+ 
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        text=element_text(size=18))

p1 <- ggplot_gtable(ggplot_build(classhydro_all))
p2 <- ggplot_gtable(ggplot_build(classhydro_facet))
lay= t(c(1,1,2,2))
png(file.path(outdirclass,'7class_hydrograph.png'),width = 16, height=9,units='in',res=300)
print(grid.arrange(p1,p2, ncol=2, layout_matrix = lay))
dev.off()

######################################################## Class boxplots ################
classHIT <- merge(HITo15y, classr7sub3_df, by="ID")
#Plot subset of metrics by name for selection of illustrative ones
# classHITplot_sub <-ggplot(setDT(classHIT)[(classHIT$indice %like% "ml"),], aes(x=as.factor(gclass), y=value, color=as.factor(gclass))) + 
#   scale_y_log10(name='Metric value') +
#   geom_boxplot() +
#   scale_x_discrete(name='Metric number (Appendix 1)')+
#   theme_classic() +
#   facet_grid(~indice, scale='free_y') +
#   theme(axis.title = element_text(size=9),
#         axis.text.y = element_text(size=9),
#         strip.text = element_text(size = 9),
#         legend.position='none')
# classHITplot_sub

HITselplot <- c('ma41','ma8','ml14','mh16','fl1','fh1','dl12','dh18','tl1','th1','ra1','ra3') #Selection and ordering of hydrologic metrics
HITselplotname <- c('Annual runoff', 'Q25/Q75','Min.flow/median flow','Q10/Q50',
  'Low flood pulse count', 'High flood pulse count','Annual min. 7-day flow','# of zero flow days', 
  'Date of annual min.','Date of annual max.','Rise rate','Fall rate') #Name for selected hydrologic metrics
classHITsel <- classHIT[classHIT$indice %in% HITselplot,] #Subset metrics
classHITsel$indice <- factor(classHITsel$indice, levels = HITselplot) #Order metrics
HIT_labels<-setNames(paste(HITselplot,HITselplotname,sep=": "),HITselplot) #Set metrics labels

#Boxplot of metrics for each class
classHITplot <-ggplot(classHITsel, aes(x=as.factor(gclass), y=value+0.01, color=as.factor(gclass))) + 
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(name='Metric value',expand=c(0.05,0)) +
  scale_x_discrete(name='Hydrologic class')+
  scale_colour_manual(values=classcol) + 
  theme_classic() +
  theme(axis.title = element_text(size=14),
        axis.text.y = element_text(size=12, angle=90),
        strip.text = element_text(size=10),
        legend.position='none') +
  facet_wrap(~indice, scales='free',ncol=4,labeller=as_labeller(HIT_labels)) 
png(file.path(outdirclass,'7class_boxplot.png'),width = 8.5, height=11.5,units='in',res=300)
print(classHITplot)
dev.off()

head(classHIT)
classHIT_summary<- setDT(classHIT)[,classmean:=mean(value, na.rm=T), .(indice, gclass)] #Get mean value of hydrologic metric for each class
classHIT_summary <- dcast(classHIT_summary, indice~gclass, value.var='classmean', fun=mean)
write.csv(classHIT_summary, file.path(outdirclass, 'classHIT_summary.csv'), row.names=F)

#To do:
#- Table reporting mean (SD) for each metric for all classes
#- Kruskal Wallis statistics for each hydrologic metrics to assess importance on classification
#- NMDS (Ward’s distance) of gauges, superimposing the dominant metrics

############################################### Format environmental data to be used in predictions##################################
####Make subset of data/remove uneeded columns
colnames(rufienv)
outcols <- c(1:5,7,8,10,46:71,94:114, which(colnames(rufienv) %in% c('CatFlowAcc','CatElvMin','CatDen','CatDamDen','CatFlowAcc','CatLCMaj',
                                                                     'WsPAPer','WsDamDen','WsGeolMaj','WsLCMaj','ReaElvMin',
                                                                     'ReaElvMax','SUM_LENGTH_GEO','Shape_Leng') |
                                              !is.na(str_match(colnames(rufienv),'DirSum*'))))
rufienvsub <- rufienv[,-outcols]
rufienvsub$ReaDirMaj <- as.factor(rufienvsub$ReaDirMaj)
#colnames(rufienvsub)

####Data transformation
#str(rufienvsub)
colnames(rufienvsub)

factcol <- c(1,2,54,55,56,57,156,157,160) #Columns that should be considered as factor
#Make factor colums factors
rufienvsub[,factcol] <- sapply(rufienvsub[,factcol], as.factor) #Factorize 
#hist.plots(rufienvsub[,-factcol]) #Inspect data
logcols <- c('CatPopDen','ReaSloAvg','WsArea','WsPopDen') #Columns to be log-transformed
rufienvsub[,logcols] <- data.trans(data.frame(rufienvsub[,logcols]), method = 'log', plot = F)
sqrtcols <- c('CatAIAvg', 'CatBio14Av','CatBio17Av','CatBio19Av','CatElvMax', 'CatElvAvg','CatSloAvg','CatSloStd','CatLen_1','CatPAPer',
              'CatRoadDen','CatWatcha','CatMineDen','CatWatOcc','ReaPAPer','ReaElvAvg','WsBio14Av','WsBio17Av','WsBio19Av','WsElvMax',
              'WsElvAvg','WsEroAvg','WsSloAvg','WsSloStd','WsDen','WsRoadDen','WsWatcha','WsMineDen','WsWatOcc','WsWatSea') #Columns to be sqrt transform
rufienvsub[,sqrtcols][rufienvsub[,sqrtcols]<0] <- 0 #A few precipitation values are negative, correct back to 0
rufienvsub[,sqrtcols] <- data.trans(rufienvsub[,sqrtcols], method = 'power',exp=.5, plot = F)

#Do not transform proportional data
# logitcols <- c('CatFLosSum_1', paste('LCSum',c(1,2,3,4,5,6,7,8,10,12,23,34,45,56,67,78,89,'10_11'),sep='_'),'CatWatExt','CatResInd','CatLakInd','WsFLosSum_1',
#               'WsWatExt','WsResInd','WsLakInd', 'WsVegPer', 'WsAgriPer') #Columns to be logit transformed
# rufienvsub[,logitcols] <- data.trans(rufienvsub[,logitcols], method = 'logit', plot = F)


###Standardization to mean of 0 and unit variance by variable
rufienvsub_std <- rufienvsub[,-factcol]
rufienvsub_std <- cbind(data.stand(rufienvsub_std, method = "standardize", margin = "column", plot = F),
                        rufienvsub[,factcol])
###Join standardized columns to gages
gagesenv_format <- gagesenv[,c('RGS_No','GridID')]
gagesenv_format <- merge(gagesenv_format,  rufienvsub_std, by='GridID')

################################################ Predict based on raw-hydro metrics classification and raw environmental predictors############################
#Set #1
pred_envar1 <-c('WsArea','CatSloAvg','CatWatExt','CatWatOcc','CatWatSea','CatDRocAvg','CatPopDen','ReaElvAvg','ReaSloAvg','WsLakInd',
                'WsBio01Av','WsBio07Av','WsBio12Av','LCSum_45')

#Set #2 - take out CatGeolMaj, WsDirMaj as not fully representing dataset
pred_envar2 <-c('WsArea','WsSloAvg','CatWatExt','CatWatOcc','CatWatSea','WsDRocAvg','WsPopDen','ReaElvAvg','WsLakInd','WsBio01Av','WsBio04Av',
               'WsBio07Av','WsBio08Av','WsBio09Av','WsBio12Av','WsBio13Av','WsBio14Av','WsBio15Av','WsBio16Av','WsBio17Av','WsBio18Av','WsBio19Av',
               'WsAIAvg','WsPermAvg','WsPoroAvg','LCSum_12','LCSum_23','LCSum_34','LCSum_45','LCSum_67','LCSum_78','LCSum_89','WsFLosSum_1',
               colnames(gagesenv)[which(!is.na(str_match(colnames(gagesenv),'GeolSum*')))][30:52])

#Set #3 - take out CatGeolMaj, WsDirMaj was not fully representing dataset and try without geology
pred_envar3 <-c('WsArea','WsSloAvg','CatWatExt','CatWatOcc','CatWatSea','WsDRocAvg','WsPopDen','ReaElvAvg','WsLakInd','WsBio01Av','WsBio04Av',
               'WsBio07Av','WsBio08Av','WsBio09Av','WsBio12Av','WsBio13Av','WsBio14Av','WsBio15Av','WsBio16Av','WsBio17Av','WsBio18Av','WsBio19Av',
               'WsAIAvg','WsPermAvg','WsPoroAvg','LCSum_12','LCSum_23','LCSum_34','LCSum_45','LCSum_67','LCSum_78','LCSum_89','WsFLosSum_1')
pred_envarname3 <- c("area", ' average slope', 'catchment water extent', 'catchment water occurrence', "catchment water seasonality",
                     " average depth to bedrock",' average pop. density', 'reach elevation',' lake index',' mean annual temp.',
                     " seasonality (temp. sd)", " temp. annual range", " mean temp. of wettest quarter", " mean temp. of driest quarter",
                     " annual precipitation", " precip. of wettest month", " precip. of driest month", " precip. seasonality"," precip. of wettest quarter",
                     " precip. of driest quarter", " precip. of warmest quarter", " precip. of coldest quarter", " average aridity index",
                     " average subsoil permeability", " average subsoil porosity", " percentage tree cover", " percentage shrub cover", " percentage grassland",
                     " percentage cropland", " percentage sparse vegetation", " percentage bare areas", " percentage urban areas", 
                     " percentage forest loss 2000-2016")

#Set #4 (Julian's selection) 
pred_envar <- c('ReaElvAvg','WsArea','WsDen','WsElvAvg','WsSloAvg','WsWatOcc','WsWatSea','WsBio10Av','WsBio11Av','WsBio12Av',
               'WsBio15Av','WsBio16Av','WsBio17Av','WsPETAvg','WsDRocAvg','WsPermAvg','WsPoroAvg','WsVegPer','WsAgriPer', 'LCSum_89','WsLakInd')
#Get labels for variable importance plot
pred_envarname <- c('reach elevation', "area", "drainage density","watershed average elevation","watershed average slope", 
                     'watershed water occurrence', "watershed water seasonality"," Mean temperature of warmest quarter", 
                     "Mean temperature of coldest quarter", "Annual precipitation", " precip. seasonality"," precip. of wettest quarter",
                     " precip. of driest quarter", " potential evapotranspiration"," average depth to bedrock"," average subsoil permeability", 
                     " average subsoil porosity", " percentage vegetation cover", "percentage agricultural land cover", " percentage urban cover",
                     ' lotic index')
pred_envlabel <- data.frame(var=pred_envar,label=pred_envarname) #Prepare labels

#Format data for prediction
gagesenvsel <- gagesenv_format[gagesenv_format$RGS_No %in% unique(rufidat_select_o15y$ID),] #Subset gauges environmental data
gagesenv_r7 <- merge(gagesenvsel,classr7sub3_df, by.x='RGS_No', by.y='ID') #Merge with class assignment df
rownames(gagesenv_r7) <- gagesenv_r7$RGS_No
gagesenv_r7 <- gagesenv_r7[,-which(colnames(gagesenv_r7) %in% c('RGS_No','GridID'))]  
gagesenv_r7$gclass <- as.factor(gagesenv_r7$gclass) #Factorize gclass
rownames(rufienvsub_std) <- rufienvsub_std$GridID
rufienvsub_std <- rufienvsub_std[,-which(colnames(rufienvsub_std) %in% 'GridID')] 

#Single tree
cat.r7r <- rpart(gclass~., data=gagesenv_r7[,c('gclass',pred_envar)], method='class',control=rpart.control(minsplit=2, minbucket=2, cp=0.025))
summary(cat.r7r)
prp(cat.r7r, col=classcol)
#rpart.plot(cat.r7r, cex=0.8, type=3, extra=1,box.palette = classcol[c(1,1,1,1,1,1,1)]) #To troubleshoot

#Boosted tree
adaboost.r7r <- boosting(gclass~., data=gagesenv_r7[,c('gclass',pred_envar)], boos=TRUE, mfinal=2000,  control=rpart.control(minsplit=2, minbucket=2, cp=0.05))
varimp <- data.frame(imp=adaboost.r7r$imp[order(adaboost.r7r$imp, decreasing = TRUE)])
varimp$var <- rownames(varimp)
varimp <- merge(varimp, pred_envlabel, by='var')
png(file.path(outdirclass,'7class_predict_envar3_imp.png'),width = 6, height=4,units='in',res=300)
print(
  ggplot(varimp[varimp$imp>0,],aes(x=reorder(label, -imp),y=imp)) + geom_bar(stat='identity') +
    theme_classic()+
    theme(axis.text.x=element_text(angle=45, hjust=1, size=10)) + 
    scale_y_continuous(name='Variable relative importance (%)', expand=c(0,0)) +
    scale_x_discrete(name="Variable name")
)
dev.off()
#CV boosted tree
adaboostcv.r7r <- boosting.cv(gclass~., data=gagesenv_r7[,c('gclass',pred_envar)], boos=TRUE, mfinal=1000,  
                              control=rpart.control(minsplit=2, minbucket=2, cp=0.05), v=14)
adaboostcv.r7r

#Predict and output 
rufi_r7r<- predict.boosting(adaboost.r7r, newdata=rufienvsub_std[,pred_envar], newmfinal=length(adaboost.r7r$trees))
#rufi_r7rmaxprob <- adply(rufi_r7r$prob, 1, max)
#qplot(rufi_r7rmaxprob$V1)
rufi_r7r_pred <- data.frame(GridID=as.integer(rownames(rufienvsub_std)),gclass=rufi_r7r$class)
rufi_r7r_pred_env <- merge(rufi_r7r_pred, rufienv, by='GridID')

#Identify those areas of the network where environmental variables are outside of gauges' range
#Keep all parts of the network within a standard deviation of the minimum and maximum values of all variables
#used in the predictions
rangesubset <- function(subset_df, range_df, classcol, subset_cols) {
  subset_df$gclasssub <- as.numeric(subset_df[,classcol])
  subset_df$subvar <- NA
  for (i in subset_cols) {
    if (i <= ncol(subset_df)) {
      varname <- colnames(subset_df)[i]
      print(varname)
      if (varname %in% colnames(range_df)) {
        if (is.numeric(subset_df[,i]) & is.numeric(range_df[,varname])) {
          subset_df[(subset_df[,i] <= (min(range_df[,varname], na.rm=T) - sd(range_df[,varname], na.rm=T))  | 
                      subset_df[,i] >= (max(range_df[,varname], na.rm=T)+ sd(range_df[,varname], na.rm=T))) & 
                      !is.na(subset_df[,i]), 'gclasssub'] <- 0
          subset_df[(subset_df[,i] <= (min(range_df[,varname], na.rm=T) - sd(range_df[,varname], na.rm=T))  | 
                       subset_df[,i] >= (max(range_df[,varname], na.rm=T)+ sd(range_df[,varname], na.rm=T))) & 
                      !is.na(subset_df[,i]), 'subvar'] <- varname
        } else {warning('Selected column is not numeric')}
      } else {warning('Column names do not match between data frames')}
    } else {warning('Column index out of subset_df range')}
  }
  return(subset_df)
}
colnames(rufi_r7r_pred_env[,c('GridID','gclass',pred_envar)])
rufi_r7r_predsub <- rangesubset(rufi_r7r_pred_env[,c('GridID','gclass',pred_envar)], 
                                gagesenvrec[gagesenvrec$RGS_No %in% unique(rufidat_select_o15y$ID),pred_envar], 
                                'gclass', 3:22)
length(which(rufi_r7r_predsub$gclasssub>0))
write.dbf(rufi_r7r_predsub, file.path(outdirclass, "predict_r7r_env4sub.dbf"))