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

# library(zoo)
# library(stargazer)

library(stringr)
library(fastcluster)
library(cluster)

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
tzenv[,'WsVegPer'] <-  with(tzenv, WsTreePer+WsShruPer+WsGrasPer) #Sum % trees, scrubs, grassland into a vegetation variable
tzenv[,'WsAgriPer'] <-  with(tzenv, WsCropPer+WsBarePer) #Sum % cropland and bare ground into an agriculture variable
tzenv[tzenv$WsVegPer>1,'WsVegPer']<- 1 #Some rounding must have led to insignificant exceedance of 1
tzenv[tzenv$WsAgriPerr>1,'WsAgriPer'] <- 1

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
factcol <- c('GridID','ReaOrd','CatDirMaj','CatGeolMaj','CatSoilMaj','WsDirMaj','WsSoilMaj','ReaDirMaj') #Columns that should be considered as factor
#Make factor colums factors
tzenvsub[,factcol] <- sapply(tzenvsub[,factcol], as.factor) #Factorize 
#hist.plots(tzenvsub[,-factcol]) #Inspect data
logcols <- c('CatPopDen','ReaSloAvg','WsArea','WsPopDen','WsAIAvg') #Columns to be log-transformed
tzenvsub[,logcols] <- data.trans(data.frame(tzenvsub[,logcols]), method = 'log', plot = F)
sqrtcols <- c('CatAIAvg', 'CatBio13Av','CatBio14Av','CatBio16Av','CatBio17Av','CatBio18Av','CatBio19Av','CatElvMax','CatDRocAvg', 'CatElvAvg','CatSloAvg','CatSloStd',
              'CatEroAvg','CatRoadDen','CatWatcha','CatMineDen','CatWatOcc','ReaPAPer','ReaElvAvg','WsBio13Av','WsBio14Av','WsBio17Av','WsBio19Av','WsElvMax',
              'WsElvAvg','WsEroAvg','WsSloAvg','WsSloStd','WsDen','WsRoadDen','WsWatcha','WsMineDen','WsWatOcc','WsWatSea','WsDRocAvg','ReaElvAvg') #Columns to be sqrt transform
tzenvsub[,sqrtcols][tzenvsub[,sqrtcols]<0] <- 0 #A few precipitation values are negative, correct back to 0
tzenvsub[,sqrtcols] <- data.trans(tzenvsub[,sqrtcols], method = 'power',exp=.5, plot = F)

###Standardization to mean of 0 and unit variance by variable
tzenvsub_std <- tzenvsub[,which(!(colnames(tzenvsub) %in% factcol))]
tzenvsub_std <- cbind(data.stand(tzenvsub_std, method = "standardize", margin = "column", plot = F),
                        tzenvsub[,factcol])
###Join standardized columns to gages
gagesenv_format <- gagesenv[,c('RGS_No','GridID')]
gagesenv_format <- merge(gagesenv_format,  tzenvsub_std, by='GridID')

################## Subset and format environmental variables ##########################
#Set #4 (Julian's selection without drainage density) 
pred_envar <- c('ReaElvAvg','WsArea','WsElvAvg','WsSloAvg','WsWatOcc','WsWatSea','WsBio10Av','WsBio11Av','WsBio12Av',
                'WsBio15Av','WsBio16Av','WsBio17Av','WsPETAvg','WsDRocAvg','WsPermAvg','WsPoroAvg','WsVegPer','WsAgriPer', 'WsUrbPer')
#Get labels for variable importance plot
pred_envarname <- c('Reach elevation', "Catchment area","Average elevation","Average slope", 
                    'Water occurrence', "Water seasonality"," Temp. warmest quarter", 
                    "Temp. coldest quarter", "Annual rainfall", " Precip. seasonality"," Precip. wettest quarter",
                    "Rainfall driest quarter", " Potential ET"," Depth to bedrock"," Subsoil permeability", 
                    "Subsoil porosity", "Vegetation % cover", "Agricultural % cover", "Urban % cover")
tzenvsub_std[,pred_envar] <- replace.missing(tzenvsub_std[,pred_envar], method='mean')
summary(tzenvsub_std[,pred_envar])
rm("tzenv","tzenvsub")
gc()
################## Classify #########################################
tzward <- fastcluster::hclust.vector(tzenvsub_std[,pred_envar], method='ward', metric="euclidean")

cluster_diagnostic <- function(clusterres, clusname) {
  #hclus.table(clusterres)
  print(coef.hclust(clusterres)) #Compute agglomerative coefficient
  #Scree plot
  png(file.path(outdir, paste('deductive_scree',clusname,'.png',sep="")), width=8, height=8, units='in',res=300)
  hclus.scree(clusterres, xlim=c(0,30)) 
  dev.off()
  #Plot dendogram
  png(file.path(outdir, paste('deductive_ward_dendogram',clusname,'.png',sep="")), width=40, height=40, units='in',res=300)
  plot(clusterres, main=paste(clusname, "gauge cluster dendogram",sep=" "), xlab='Gauge ID', ylab="Gower's distance", hang=-1)  
  rect.hclust(clusterres, k=5) 
  rect.hclust(clusterres, k=7) 
  rect.hclust(clusterres, k=8) 
  rect.hclust(clusterres, k=10)
  dev.off()
}
cluster_diagnostic(tzward, "Ward's D")

kclass=10
classr <-cutree(tzward, k=10)
classr_df <- data.frame(GridID=tzenvsub_std$GridID, gclass=classr) 
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}
write.csv(classr_df, file.path(outdir,paste0(kclass,'classtab.csv')), row.names=F)



gauge_gow<- gowdis(tzenvsub_std, w=rep(1,ncol(tzenvsub_std)), asym.bin = NULL) #Compute Gower's distance so that missing values will not be taken in account


gaugecla_ward <-hclust(gaugegow_o15y, method='ward.D') #Classify using hierarchical agglomerative using Ward's minimum variance method
gaugecla_ward2 <-hclust(gaugegow_o15y, method='ward.D2') #Classify using hierarchical agglomerative using Ward's minimum variance method
gaugecla_UPGMA <-hclust(gaugegow_o15y, method='average') #Classify using  UPGMA

#################### Make table ###############################
classenv <- merge(tzenv[,c('GridID',pred_envar)], classr_df, by="GridID")
classenv_melt <- melt(classenv, id.vars=c('GridID','gclass'),value.name='Value')
#pred_envar[which(!(pred_envar %in% colnames(tzenv)))]

df <- classenv_melt

classtableformat <- function(df, tabname) {
  classenv_stats<- setDT(df)[,`:=`(classmean=mean(value, na.rm=T),classsd=sd(value,na.rm=T)), .(variable, gclass)] #Get mean and SD of hydrologic metric for each class
  classenv_stats <- classenv_stats[!duplicated(classenv_stats[,c('variable','gclass')]),]
  classenv_statsmean <- as.data.frame(dcast(classenv_stats, gclass~variable, value.var='classmean'))
  classenv_statssd <- as.data.frame(dcast(classenv_stats, gclass~variable, value.var='classsd'))
  #Format digits for mean and sd
  classenv_stats_meanform <- melt(setDT(digitform(classenv_statsmean,
                                                  cols=2:(ncol(classenv_statsmean)), extradigit=1)),id.vars='gclass',variable.name='Variable', value.name='classmean')
  classenv_stats_sdform <- melt(setDT(digitform(classenv_statssd,
                                                cols=2:(ncol(classenv_statssd)), extradigit=1)),id.vars='gclass',variable.name='Variable', value.name='classsd')
  classenv_format <- merge(classenv_stats_meanform, classenv_stats_sdform, by=c('gclass','Variable'))
  classenv_format$tabcol <- with(classenv_format, paste0(classmean,' (',classsd,')'))
  
  classenv_format[is.na(classenv_format$classsd),'tabcol']  <- '–'
  table <- as.data.frame(dcast(classenv_format, Variable~gclass, value.var='tabcol'))
  table <- table[with(Variables, order(type, num)),] #Order rows by Variables number
  
  kable(table) %>%
    kable_styling(bootstrap_options = "striped", font_size = 10) %>%
    row_spec(which(table$Metric %in% envo15ysub3$indice), bold=T) %>%
    column_spec(2,color=classcol[1]) %>%
    column_spec(3,color=classcol[2]) %>%
    column_spec(4,color=classcol[3]) %>%
    column_spec(5,color=classcol[4]) %>%
    column_spec(6,color=classcol[5]) %>%
    column_spec(7,color=classcol[6]) %>%
    column_spec(8,color=classcol[7]) %>%
    save_kable(tabname, self_contained=T)
}
classtableformat(classenv, KWtab=metricKW, tabname='deductive_envgclass_cast.doc')

classcol=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
           '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')

env_labels<-setNames(pred_envarname,pred_envar) #Set metrics labels
#Boxplot of metrics for each class
classenvplot <-ggplot(classenv_melt[!(classenv_melt$variable %in% c('WsArea','ReaElvAvg','WsWatSea')),], aes(x=as.factor(gclass), y=value+0.01, color=as.factor(gclass))) + 
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(name='Metric value',expand=c(0.05,0)) +
  scale_x_discrete(name='Hydrologic class')+
  scale_colour_manual(values=classcol) + 
  theme_classic() +
  theme(axis.title = element_text(size=14),
        axis.text.y = element_text(size=12, angle=90),
        strip.text = element_text(size=12),
        legend.position='none') +
  facet_wrap(~variable, scales='free',ncol=4,labeller=as_labeller(env_labels))
#dir='classo15y_ward_rawsub3'
#outdirclass <- file.path(outdir,dir)
png(file.path(outdir,'10class_boxplot.png'),width = 9, height=12.5,units='in',res=300)
print(classenvplot)
dev.off()
