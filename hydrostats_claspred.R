#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/30/2018
#Date last updated: 04/02/2018

#Purpose: compute hydrological metrics and classify stream gauges

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
datadir = file.path(getwd(),'rufiji_hydrodatafilter') #UPDATE
origdatadir = file.path(rootdir,"data") 
outdir=file.path(getwd(),'rufiji_hydrodatastats')
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}

rufidat_clean <- read.csv(file.path('rufiji_hydrodatainspect','rufidat_clean.csv'), colClasses=c('factor','Date','numeric','character','character'))
rufidat_impute <- read.csv(file.path('rufiji_hydrodataimpute', 'rufidat_interp.csv'), colClasses=c('Date',rep(c('numeric','numeric','character'),39)))
colnames(rufidat_impute)[seq(2,ncol(rufidat_impute),3)] <- substr(colnames(rufidat_impute),2,10)[seq(2,ncol(rufidat_impute),3)]
rufidat_gapsummary <- read.csv(file.path(datadir, 'rufidat_gapsummary.csv'))
rufidat_post1991<-read.csv(file.path(datadir, 'gageselect_post1991comp90.csv'))
rufidat_o15y<-read.csv(file.path(datadir, 'gageselect_o15comp90.csv'))

gagesenv <- read.csv(file.path(getwd(),'gages_netjoinclean.csv'))
gagesenvrec <- merge(gagesenv, unique(rufidat_clean[,c('ID','SYM')]), by.x='RGS_No', by.y='ID', all.x=F)
rufienv <- read.csv(file.path(getwd(),'streamnet118_rufiji_finaltabclean.csv'))

predsmelt <-melt(setDT(as.data.frame(rufidat_impute)[,c(1,which(colnames(rufidat_impute) %like% "^1K*"))]),id.vars = 'Date',value.name='Flow',variable.name='ID')
predsmelt <- predsmelt[,c(2,1,3)]
predsmelt$year <- as.numeric(format(predsmelt$Date, "%Y"))
predsmelt$month <- as.numeric(format(predsmelt$Date, "%m"))
predsmelt<-hyear.internal(predsmelt,hyrstart=10) #hdoy doesn't work?
predsmelt$doy <- as.numeric(format(predsmelt$Date,"%j"))
predsmelt[,hdoy:=ifelse(month>=10,
                        doy-as.numeric(format(as.Date(paste(year,'-10-01',sep="")),"%j")),
                        doy+as.numeric(format(as.Date(paste(year-1,'-12-31',sep="")),"%j"))-as.numeric(format(as.Date(paste(year-1,'-10-01',sep="")),"%j")))] #Compute hydrologic day
predsmelt <- merge(predsmelt, rufidat_gapsummary, by=c('ID','hyear'),all.x=T)
predsmelt <- merge(predsmelt, rufidat_post1991, by='ID',all.x=T)
predsmelt <- merge(predsmelt, rufidat_o15y, by='ID',all.x=T)

############################################### Compute hydrologic metrics##########################################
allHITcomp <- function(dfhydro, dfenv, gageID, templateID='1KA9',hstats="all", floodquantile=0.95) {
  #Get template
  dailyQClean <- validate_data(dfhydro[dfhydro$ID==templateID,c("Date", "Flow")], yearType="water")
  #Calculate all hit stats
  HITall_template <- calc_allHIT(dailyQClean, yearType="water", stats=hstats, digits=10, pref="mean",
                                drainArea=dfenv[dfenv$RGS_No==templateID,'WsArea'], floodThreshold = quantile(dailyQClean$discharge, floodquantile))
  colnames(HITall_template)[2] <- templateID
  HITall <- data.frame(indice=HITall_template$indice) 
  #Compute metrics for all gages
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

rufidat_select_o15y <- predsmelt[predsmelt$gap_per<=0.1 & predsmelt$max_gap < 37 & predsmelt$hyear<2017 & 
                                   (predsmelt$ycount_o15>=15 | predsmelt$ID=='1KB32') | 
                                   (predsmelt$ID=='1KA41' & predsmelt$max_gap < 273 & predsmelt$gap_per<0.75 & predsmelt$hyear<1996) |
                                   (predsmelt$ID=='1KA42A' & predsmelt$max_gap < 273 & predsmelt$gap_per<0.75 & predsmelt$hyear>1958 & predsmelt$hyear<2017),]
HITo15y <- allHITcomp(as.data.frame(rufidat_select_o15y), gagesenv, 'ID')
#write.csv(HITo15y, file.path(outdir, 'HITo15y.csv'),row.names=F)

############################################### Compute redundancy in hydrologic metrics#############################

############################################### Box plot of metrics##############################
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

########################################Format environmental data to be used in predictions##################################
#Make subset of data/remove uneeded columns
colnames(rufienv)
outcols <- c(1:5,7,8,10,46:71,94:114, which(colnames(rufienv) %in% c('CatFlowAcc','CatElvMin','CatDen','CatDamDen','CatFlowAcc','CatLCMaj',
                                                                    'WsPAPer','WsDamDen','WsGeolMaj','WsLCMaj','ReaElvMin',
                                                                    'ReaElvMax','SUM_LENGTH_GEO','Shape_Leng') |
                                             !is.na(str_match(colnames(rufienv),'DirSum*'))))
rufienvsub <- rufienv[,-outcols]
rufienvsub$ReaDirMaj <- as.factor(rufienvsub$ReaDirMaj)
#colnames(rufienvsub)

#Data transformation
#Transform catchment variables
str(rufienvsub)
colnames(rufienvsub)
factcol <- c(1,2,54,55,56,57,156,157,160)
#Make factor colums factors
rufienvsub[,factcol] <- sapply(rufienvsub[,factcol], as.factor)
#hist.plots(rufienvsub[,-factcol]) #Inspect data
logcols <- c('CatPopDen','ReaSloAvg','WsArea','WsPopDen')
rufienvsub[,logcols] <- data.trans(data.frame(rufienvsub[,logcols]), method = 'log', plot = F)
asincols <- c('CatFLosSum_1', paste('LCSum',c(1,2,3,4,5,6,7,8,10,12,23,34,45,56,67,78,89,'10_11'),sep='_'),'CatWatExt','CatResInd','CatLakInd','WsFLosSum_1',
              'WsWatExt','WsResInd','WsLakInd')
rufienvsub[,asincols] <- data.trans(rufienvsub[,asincols], method = 'asin', plot = F)
sqrtcols <- c('CatAIAvg', 'CatBio14Av','CatBio17Av','CatBio19Av','CatElvMax', 'CatElvAvg','CatSloAvg','CatSloStd','CatLen_1','CatPAPer',
              'CatRoadDen','CatWatcha','CatMineDen','CatWatOcc','ReaPAPer','ReaElvAvg','WsBio14Av','WsBio17Av','WsBio19Av','WsElvMax',
              'WsElvAvg','WsEroAvg','WsSloAvg','WsSloStd','WsDen','WsRoadDen','WsWatcha','WsMineDen','WsWatOcc','WsWatSea')
rufienvsub[,sqrtcols] <- data.trans(rufienvsub[,sqrtcols], method = 'power',exp=.5, plot = F)
#Transform watershed variables
#hist.plots(envsubws) #Inspect data

#Then standardize to mean of 0 and unit variance
rufienvsub_std <- rufienvsub[,-factcol]
rufienvsub_std <- cbind(data.stand(rufienvsub_std, method = "standardize", margin = "column", plot = F),
                       rufienvsub[,factcol])
#Join standardized columns to gages
gagesenv_format <- gagesenv[,c('RGS_No','GridID')]
gagesenv_format <- merge(gagesenv_format,  rufienvsub_std, by='GridID')

########################################Format hydrologic metrics to use in classification and compute Gower's distance############
HITdist <- function(HITdf, logmetrics) {
  if (logmetrics==TRUE){
    HITdf$Value <- log(HITdf$Value+1)
  }
  HITdf_format <- dcast(HITdf, ID ~ indice)
  HITdf_format <- merge(HITdf_format, gagesenvrec[,c('RGS_No','WsArea')], by.x='ID', by.y='RGS_No')
  dimindices <- c('ma1','ma2',paste('ma',seq(12,23),sep=''),paste('ml',seq(1,12),sep=''),paste('mh',seq(1,12),sep=''), 
                  paste('dl',seq(1,5),sep=''),paste('dh',seq(1,5),sep=''),'ra1','ra3','ra6','ra7') #List of dimensional indices from Kennen et al. 2007
  HITdf_format <- as.data.frame(HITdf_format[,(dimindices) := lapply(.SD, function(x) round(x/WsArea, digits=10)), .SDcols=dimindices]) #Standardize dimensional indices by drainage area
  row.names(HITdf_format) <- HITdf_format$ID
  HITdf_format <- HITdf_format[,-which(colnames(HITdf_format) %in% c('ID','WsArea'))] #Get rid of non-indices columns
  HITdf_stand <- data.stand(HITdf_format[,2:(ncol(HITdf_format))],method='standardize',margin='column',plot=F) #z-standardize columnwise 
  gauge_gow<- gowdis(HITdf_stand, w=rep(1,ncol(HITdf_stand)), asym.bin = NULL) #Compute Gower's distance so that missing values will not be taken in account
  return(gauge_gow)
}

########################################CLASSIFICATION BASED ON POST-1991 > 10 YEARS OF DATA ################################
########################################Classify based on raw indices############################################
gaugegow1991 <- HITdist(HITpost1991)
gaugecla_ward <-hclust(gaugegow1991, method='ward.D') #Classify using hierarchical agglomerative using Ward's minimum variance method

#Diagnostics
hclus.table(gaugecla_ward)
plot(gaugecla_ward, main="Ward's distance gauge cluster dendogram", xlab='Gauge ID', ylab="Gower's distance", hang=-1)   #par(mfrow=c(1,1))
coef.hclust(gaugecla_ward) #Compute agglomerative coefficient
cor(gauge_gow, cophenetic(gaugecla_ward)) #Compute cophenetic coefficient
hclus.cophenetic(gaugegow1991, gaugecla_ward) #Plot cophenetic relationship 
hclus.scree(gaugecla_ward) #Check out scree plot
plot(gaugecla_ward, main="Ward's distance gauge cluster dendogram", xlab='Gauge ID', ylab="Gower's distance", hang=-1)   #par(mfrow=c(1,1))
rect.hclust(gaugecla_ward, k=5) #Draw rectangle around k classes
#Test significance of classes
clus.stab <- pvclust(t(HITall_stand), method.hclust='ward.D', method.dist='cor',use.cor="pairwise.complete.obs", nboot=4999)
plot(clus.stab)
pvrect(clus.stab, alpha=0.90)
#Get gauge classes
classr5 <-cutree(gaugecla_ward, k=5)
classr5_df <- data.frame(ID=names(classr5), gclass=classr5) 
outdirclass <- file.path(outdir,'class1991_ward_raw')
if (dir.exists(outdirclass )) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdirclass))
  dir.create(outdirclass )
}
write.csv(classr5, file.path(outdirclass,'class_rawgow_ward_5.csv'), row.names=F)
#Boxplots
HITall_new <- cbind(classr5,HITall_format)
#box.plots(HITall_new, by='classr5')

#rufidat_select_classr5 <- merge(rufidat_select, classr5_df, by="ID")
#write.csv(rufidat_select_classr5, file.path(outdir, 'class_ward_raw/rufidat_select_classr5.csv'), row.names=F)

########################################Predict based on raw-hydro metrics classification and raw environmental predictors############################
pred_envar1 <-c('WsArea','CatSloAvg','CatWatExt','CatWatOcc','CatWatSea','CatDRocAvg','CatPopDen','ReaElvAvg','ReaSloAvg','WsLakInd',
                'WsBio01Av','WsBio07Av','WsBio12Av','LCSum_45')
gagesenv_r5 <- merge(gagesenvrec[gagesenvrec$RGS_No %in% rufidat_select$ID, c('RGS_No',pred_envar1)],classr5_df, by.x='RGS_No', by.y='ID')
rownames(gagesenv_r5) <- gagesenv_r5$RGS_No
gagesenv_r5 <- gagesenv_r5[,-which('RGS_No' %in% colnames(gagesenv_r5))]
gagesenv_r5$gclass <- as.factor(gagesenv_r5$gclass)
rufienvsel_r5r <- rufienv[,c('GridID',pred_envar)]
rownames(rufienvsel_r5r) <- rufienvsel_r5r$GridID
rufienvsel_r5r <- rufienvsel_r5r[,-which('GridID' %in% colnames(rufienvsel_r5r))]
#Single tree
cat.r5r <- rpart(gclass~., data=gagesenv_r5, method='class',control=rpart.control(minsplit=1, minbucket=1, cp=0.05))
summary(cat.r5r)
rpart.plot(cat.r5r)
#Boosted tree
adaboost.r5r <- boosting(gclass~., data=gagesenv_r5, boos=TRUE, mfinal=100,  control=rpart.control(minsplit=1, minbucket=1, cp=0.1))
barplot(adaboost.r5r$imp[order(adaboost.r5r$imp, decreasing = TRUE)],
        ylim = c(0, 100), main = "Variables Relative Importance",
        col = "lightblue")
adaboostcv.r5r <- boosting.cv(gclass~., data=gagesenv_r5, boos=TRUE, mfinal=10,  control=rpart.control(minsplit=1, minbucket=1, cp=0.1), v=19)
adaboostcv.r5r
#Predict
rufi_r5r<- predict.boosting(adaboost.r5r, newdata=rufienvsel_r5r, newmfinal=length(adaboost.r5r$trees))
#rufi_r5rmaxprob <- adply(rufi_r5r$prob, 1, max)
#qplot(rufi_r5rmaxprob$V1)
rufi_r5r_pred <- cbind(rufienvsel_r5r,gclass=rufi_r5r$class)
rufi_r5r_pred$GridID <- as.integer(rownames(rufi_r5r_pred)) 
write.dbf(rufi_r5r_pred[,c('GridID','gclass')], file.path(outdir, "class_ward_raw/predict_r5r.dbf"))


######################################### CLASSIFICATION BASED ON ENTIRE PERIOD > 15 YEARS OF DATA ################################
################################################Classify based on raw indices############################################
gaugegow_o15y <- HITdist(HITo15y, logmetrics=TRUE)
gaugecla_ward <-hclust(gaugegow_o15y, method='ward.D') #Classify using hierarchical agglomerative using Ward's minimum variance method
gaugecla_ward2 <-hclust(gaugegow_o15y, method='ward.D2') #Classify using hierarchical agglomerative using Ward's minimum variance method
gaugecla_UPGMA <-hclust(gaugegow_o15y, method='average') #Classify using hierarchical agglomerative using Ward's minimum variance method

#Diagnostics
cluster_diagnostic <- function(clusterres, clusname, gowdis) {
  #hclus.table(clusterres)
  coef.hclust(clusterres) #Compute agglomerative coefficient
  #cor(gowdis, cophenetic(clusterres)) #Compute cophenetic coefficient
  png(file.path(outdir, paste('o15y_r6r_cophe',clusname,'.png',sep="")), width=8, height=8, units='in',res=300)
  hclus.cophenetic(gowdis, clusterres) #Plot cophenetic relationship 
  dev.off()
  png(file.path(outdir, paste('o15y_r6r_scree',clusname,'.png',sep="")), width=8, height=8, units='in',res=300)
  hclus.scree(clusterres) #Check out scree plot
  dev.off()
  png(file.path(outdir, paste('o15y_r6r_dendogram',clusname,'.png',sep="")), width=8, height=8, units='in',res=300)
  plot(clusterres, main=paste(clusname, "gauge cluster dendogram",sep=" "), xlab='Gauge ID', ylab="Gower's distance", hang=-1)   #par(mfrow=c(1,1))
  rect.hclust(clusterres, k=4) #Draw rectangle around k classes
  rect.hclust(clusterres, k=5) #Draw rectangle around k classes
  rect.hclust(clusterres, k=6) #Draw rectangle around k classes
  rect.hclust(clusterres, k=7) #Draw rectangle around k classes
  rect.hclust(clusterres, k=8) #Draw rectangle around k classes
  dev.off()
}
cluster_diagnostic(gaugecla_ward, "Ward's D", gaugegow_o15y)
cluster_diagnostic(gaugecla_ward2, "Ward's D2", gaugegow_o15y)
cluster_diagnostic(gaugecla_UPGMA, "UPGMA", gaugegow_o15y)
#UPGMA leads to too much chaining, and Ward's D has higher cophenetic correlation and reaches an elbow after 7 (rather than 8 classes for D2)

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
classcol<- c('#1b9e77',"#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02")

classcol<- c("#176c93","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#7a5614") #7 classes with darker color

######################################################## Dendogram plot#######################################
#Get the dend
dend <- as.dendrogram(gaugecla_ward)
png(file.path(outdirclass,'7class_dendrogram.png'),width = 8, height=8,units='in',res=300)
par(mar=c(3,3,0,3)) #bottom left top right
dend %>% set("branches_lwd", 2.5) %>% 
  color_branches(k=7, col=classcol, groupLabels=T) %>% 
  #color_branches(clusters=as.numeric(temp_col), col=levels(temp_col), groupLabels=as.character(as.numeric(temp_col))) %>% 
  color_labels(k=7, col=classcol) %>%
  plot(horiz=TRUE,xlab="Gower's distance", ylab="Gauge ID",mgp=c(1.5,0.5,0))
dev.off()

gaugecla_ward_name <- gaugecla_ward
gaugecla_ward_name$labels <- with(gagesenvrec[gagesenvrec$RGS_No %in% gaugecla_ward_name$labels,], 
                                  paste(RGS_No,"-",RGS_Loc," River at ", RGS_Name,sep=""))
dendname <- as.dendrogram(gaugecla_ward_name)
png(file.path(outdirclass,'7class_dendrogram_names.png'),width = 8, height=8,units='in',res=300)
par(mar=c(3,3,0,17)) #bottom left top right
dendname %>% set("branches_lwd", 2.5) %>% 
  color_branches(k=7, col=classcol, groupLabels=T) %>% 
  #color_branches(clusters=as.numeric(temp_col), col=levels(temp_col), groupLabels=as.character(as.numeric(temp_col))) %>% 
  color_labels(k=7, col=classcol) %>%
  plot(horiz=TRUE,xlab="Gower's distance", ylab="Gauge ID",mgp=c(1.5,0.5,0))
dev.off()

######################################################## Class hydrograph plots ################
rufidat_select_classr7 <- merge(rufidat_select_o15y, classr7_df, by="ID")
#write.csv(rufidat_select_classr7, file.path(outdir, 'class_ward_raw/rufidat_select_classr7.csv'), row.names=F)
setDT(rufidat_select_classr7)[,yrmean:=mean(Flow),.(ID,hyear)]
classflowstats <- setDT(rufidat_select_classr7)[,list(classmeanfull=mean(Flow, na.rm=T), classmean= mean(Flow/yrmean,na.rm=T),classQ75= quantile(Flow/yrmean, .75,na.rm=T),
                                                      classQ25=quantile(Flow/yrmean, .25,na.rm=T),classQ90=quantile(Flow/yrmean, .90,na.rm=T),
                                                      classQ10=quantile(Flow/yrmean, .10,na.rm=T),classmax=max(Flow/yrmean,na.rm=T),
                                                      classmin=min(Flow/yrmean,na.rm=T),classsd=sd(Flow/yrmean,na.rm=T), 
                                                      cal_hdoy=format(as.Date(hdoy, origin='2015-10-01'), "%Y-%m-%d")),
                                                .(gclass,hdoy)]

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



classhydro_all <- ggplot(as.data.frame(classflowstats), aes(x=as.Date(cal_hdoy), y=classmean, color=factor(gclass))) + 
  geom_line(size=1, alpha=0.8) + 
  scale_color_manual(name='Hydrologic class',values=classcol) +
  scale_fill_manual(name='Hydrologic class',values=classcol) +
  scale_y_continuous(name='Daily mean discharge/Mean daily discharge',expand=c(0,0),limits=c(0,NA)) + 
  scale_x_date(name='Date',date_breaks = "1 month", date_labels = "%b", expand=c(0,0)) + 
  theme_classic() + 
  theme(legend.position='none',
        text=element_text(size=18))

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
classHIT <- merge(HITo15y, classr7_df, by="ID")
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

HITselplot <- c('ma41','ma8','ml14','mh16','fl1','fh1','dl12','dh18','tl1','th1','ra1','ra3')
HITselplotname <- c('Annual runoff', 'Q25/Q75','Min.flow/median flow','Q10/Q50',
  'Low flood pulse count', 'High flood pulse count','Annual min. 7-day flow','# of zero flow days', 
  'Date of annual min.','Date of annual max.','Rise rate','Fall rate')
classHITsel <- classHIT[classHIT$indice %in% HITselplot,]
classHITsel$indice <- factor(classHITsel$indice, levels = HITselplot)
HIT_labels<-setNames(paste(HITselplot,HITselplotname,sep=": "),HITselplot)
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
classHIT_summary<- setDT(classHIT)[,classmean:=mean(value, na.rm=T), .(indice, gclass)]
classHIT_summary <- dcast(classHIT_summary, indice~gclass, value.var='classmean', fun=mean)
write.csv(classHIT_summary, file.path(outdirclass, 'classHIT_summary.csv'), row.names=F)

################################################Predict based on raw-hydro metrics classification and raw environmental predictors############################
#Set #1
pred_envar1 <-c('WsArea','CatSloAvg','CatWatExt','CatWatOcc','CatWatSea','CatDRocAvg','CatPopDen','ReaElvAvg','ReaSloAvg','WsLakInd',
                'WsBio01Av','WsBio07Av','WsBio12Av','LCSum_45')

#Set #2 - take out CatGeolMaj, WsDirMaj as not fully representing dataset
pred_envar2 <-c('WsArea','WsSloAvg','CatWatExt','CatWatOcc','CatWatSea','WsDRocAvg','WsPopDen','ReaElvAvg','WsLakInd','WsBio01Av','WsBio04Av',
               'WsBio07Av','WsBio08Av','WsBio09Av','WsBio12Av','WsBio13Av','WsBio14Av','WsBio15Av','WsBio16Av','WsBio17Av','WsBio18Av','WsBio19Av',
               'WsAIAvg','WsPermAvg','WsPoroAvg','LCSum_12','LCSum_23','LCSum_34','LCSum_45','LCSum_67','LCSum_78','LCSum_89','WsFLosSum_1',
               colnames(gagesenv)[which(!is.na(str_match(colnames(gagesenv),'GeolSum*')))][30:52])

#Set #3 - take out CatGeolMaj, WsDirMaj was not fully representing dataset and try without geology
pred_envar <-c('WsArea','WsSloAvg','CatWatExt','CatWatOcc','CatWatSea','WsDRocAvg','WsPopDen','ReaElvAvg','WsLakInd','WsBio01Av','WsBio04Av',
               'WsBio07Av','WsBio08Av','WsBio09Av','WsBio12Av','WsBio13Av','WsBio14Av','WsBio15Av','WsBio16Av','WsBio17Av','WsBio18Av','WsBio19Av',
               'WsAIAvg','WsPermAvg','WsPoroAvg','LCSum_12','LCSum_23','LCSum_34','LCSum_45','LCSum_67','LCSum_78','LCSum_89','WsFLosSum_1')
pred_envarname <- c("area", ' average slope', 'catchment water extent', 'catchment water occurrence', "catchment water seasonality",
                    " average depth to bedrock",' average pop. density', 'reach elevation',' lake index',' mean annual temp.',
                    " seasonality (temp. sd)", " temp. annual range", " mean temp. of wettest quarter", " mean temp. of driest quarter",
                    " annual precipitation", " precip. of wettest month", " precip. of driest month", " precip. seasonality"," precip. of wettest quarter",
                    " precip. of driest quarter", " precip. of warmest quarter", " precip. of coldest quarter", " average aridity index",
                    " average subsoil permeability", " average subsoil porosity", " percentage tree cover", " percentage shrub cover", " percentage grassland",
                    " percentage cropland", " percentage sparse vegetation", " percentage bare areas", " percentage urban areas", 
                    " percentage forest loss 2000-2016")
pred_envlabel <- data.frame(var=pred_envar,label=pred_envarname)
length(pred_envar)

gagesenvsel <- gagesenv_format[gagesenv_format$RGS_No %in% unique(rufidat_select_o15y$ID),]
gagesenv_r7 <- merge(gagesenvsel,classr7_df, by.x='RGS_No', by.y='ID')
rownames(gagesenv_r7) <- gagesenv_r7$RGS_No
gagesenv_r7 <- gagesenv_r7[,-which(colnames(gagesenv_r7) %in% c('RGS_No','GridID'))]
gagesenv_r7$gclass <- as.factor(gagesenv_r7$gclass)
rownames(rufienvsub_std) <- rufienvsub_std$GridID
rufienvsub_std <- rufienvsub_std[,-which(colnames(rufienvsub_std) %in% 'GridID')]

#Single tree
cat.r7r <- rpart(gclass~., data=gagesenv_r7[,c('gclass',pred_envar)], method='class',control=rpart.control(minsplit=4, minbucket=1, cp=0.05))
summary(cat.r7r)
prp(cat.r7r, col=classcol)
rpart.plot(cat.r7r, cex=0.8, type=3, extra=1,box.palette = classcol[c(1,1,1,1,1,1,1)])

#Boosted tree
adaboost.r7r <- boosting(gclass~., data=gagesenv_r7[,c('gclass',pred_envar)], boos=TRUE, mfinal=2000,  control=rpart.control(minsplit=1, minbucket=1, cp=0.1))
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
adaboostcv.r7r <- boosting.cv(gclass~., data=gagesenv_r7[,c('gclass',pred_envar)], boos=TRUE, mfinal=100,  
                              control=rpart.control(minsplit=1, minbucket=1, cp=0.1), v=23)
adaboostcv.r7r

#Predict
rufi_r7r<- predict.boosting(adaboost.r7r, newdata=rufienvsub_std[,pred_envar], newmfinal=length(adaboost.r7r$trees))
#rufi_r7rmaxprob <- adply(rufi_r7r$prob, 1, max)
#qplot(rufi_r7rmaxprob$V1)
rufi_r7r_pred <- data.frame(GridID=as.integer(rownames(rufienvsub_std)),gclass=rufi_r7r$class)
write.dbf(rufi_r7r_pred, file.path(outdirclass, "predict_r7r_env3.dbf"))




################################## DUMP #####################################
# gaugecla_warddat <- dendro_data(as.dendrogram(gaugecla_ward), type='rectangle')
# ggplot(segment(gaugecla_warddat)) +
#   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_text(data = label(gaugecla_warddat),
#             aes(x = x, y = y, label = label), vjust = -0.5, size = 3) +
#   coord_flip() +
#   scale_y_reverse(name="Gower's distance", expand = c(0.2, 0)) +
#   scale_x_continuous(name='Gauge ID') +
#   theme_classic() +
#   theme(axis.text.y = element_blank(),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank())
