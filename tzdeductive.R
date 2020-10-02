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
library(data.table)
library(fastcluster)
library(cluster)
library(ggplot2)

rootdir="D:/Tanzania/Tanzania" #####UPDATE THIS TO MATCH YOUR ROOT PROJECT FOLDER #######

source(file.path(rootdir,"bin/outside_src/Biostats.R"))
#source(file.path(rootdir,"bin/outside_src/mjcgraphics.R"))
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
tzenv <- fread('streamnet118_final.csv') %>%
  as.data.frame

gagesenv <- read.csv(file.path(getwd(),'gages_netjoinclean.csv')) #Import gauges' environmental data


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
              'WsElvAvg','WsEroAvg','WsSloAvg','WsSloStd','WsDen','WsRoadDen','WsWatcha','WsMineDen','WsWatOcc','WsWatSea','WsDRocAvg') #Columns to be sqrt transform
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
                'WsBio15Av','WsBio16Av','WsBio17Av','WsPETAvg','WsDRocAvg','WsPermAvg','WsPoroAvg','WsVegPer','WsAgriPer', 'WsUrbPer','WsLakInd') #Add lake index
#Get labels for variable importance plot
pred_envarname <- c('Reach elevation (m)', "Catchment area (km^2)",
                    "Catchment elevation (m)","Catchment slope (°)", 
                    'Water occurrence', "Water seasonality"," Temp. warmest quarter (°C)", 
                    "Temp. coldest quarter (°C)", "Annual precipitation (mm)",
                    "Precip. seasonality"," Precip. wettest quarter (mm)",
                    "Precip. driest quarter (mm)", "Annual PET (mm)",
                    "Depth to bedrock (cm)","Subsoil permeability", 
                    "Subsoil porosity", "Vegetation cover (%)", 
                    "Agricultural cover (%)", "Built up areas (%)","Lake index")
pred_envlabel <- data.frame(var=pred_envar,label=pred_envarname) #Prepare labels
tzenvsub_std[,pred_envar] <- replace.missing(tzenvsub_std[,pred_envar], method='mean')
summary(tzenvsub_std[,pred_envar])
rm("tzenv","tzenvsub")
gc()
################## Classify #########################################
flist <- list.files()
deducimg <- grep('tz_deductive.*[.]RData', flist)
if (deducimg) {
  print('A saved classification already exists, loading...')
  load(flist[max(deducimg)])
} else {
  print('There is no existing image of the deductive classification, computing...')
  # tzward <- fastcluster::hclust.vector(tzenvsub_std[,pred_envar], method='ward', metric="euclidean")
  # save.image(file = paste0("tz_deductive", Sys.Date(),".RData"))
}

hclus.scree <- function(x,ylabel,...){
  old.par<-par(no.readonly=TRUE)
  par(ps = 14, cex = 1, cex.main = 1,mai=c(1,1,0.2,0.2))
  z1<-seq(length(x$height),1)
  z<-as.data.frame(cbind(z1,sort(x$height)))
  plot(z[,1],z[,2],type='o',lwd=1.5,pch=19,col='blue',
       main=NULL,ylab=ylabel,xlab='Number of clusters',...)
  par(old.par)
}

cluster_diagnostic <- function(clusterres, ylabel,clusname,dendo) {
  #hclus.table(clusterres)
  #print(paste0('Agglomerative coefficient: ', coef.hclust(clusterres))) #Compute agglomerative coefficient
  #print(paste0('Cophenetic correlation coefficient: ',cor(gowdis, cophenetic(clusterres)))) #Compute cophenetic coefficient
  #Scree plot
  pdf(file.path(outdir, paste('deductive_scree',clusname,'.pdf',sep="")), width=8, height=8)
  hclus.scree(clusterres, ylabel=ylabel, frame.plot=FALSE,cex=1, xlim=c(0,31), ylim=c(0,max(clusterres$height)+50), xaxs='i', yaxs='i') 
  dev.off()
}
cluster_diagnostic(tzward, ylabel="Average within-cluster dissimilarity (Euclidean distance)",clusname="Ward's D_20180706", dendo=T)

kclass=11
classr <-dendextend::cutree(tzward, k=kclass, order_clusters_as_data = T)
classr_df <- data.frame(GridID=tzenvsub_std$GridID, gclass=classr) 
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}
write.csv(classr_df, file.path(outdir,paste0(kclass,'classtab_20180701.csv')), row.names=F)
classr_df <-read.csv(file.path(outdir,paste0(kclass,'classtab_20180701.csv')))

## cut the tree into k clusters and reconstruct the upper part of the
cent <- NULL
for(k in 1:kclass){
  cent <- rbind(cent, colMeans(setDT(tzenvsub_std[,pred_envar])[classr == k, , drop = FALSE]))
}
hc1 <- hclust(dist(cent)^2, method='ward',members = table(classr))
dendname<- as.dendrogram(hc1)
colorder <- NULL
if (is.null(colorder)) colorder = 1:kclass

png(file.path(outdir, paste('deductive_ward_dendogram_k11.png',sep="")), width=8, height=8, units='in',res=300)
dendname %>% set("branches_lwd", 1) %>% 
  #color_branches(clusters=as.numeric(temp_col), col=levels(temp_col), groupLabels=as.character(as.numeric(temp_col))) %>% 
  #color_labels(k=kclass, col=classcol[colorder]) %>%
  plot(horiz=TRUE,xlab="Euclidean distance", ylab="",mgp=c(1.5,0.5,0))
dev.off()


#################### Make table ###############################
#RE-ORDER CLASSES BASED ON TREE
classes <- c('Urban', 'LkInf', 'LkInf', 'CoaWt', 'LWtSe',
             'MWtSe', 'Monta', 'AgrPl', 'DrySt', 'DrySe',
             'VDryS')
reord_k11 <- data.frame(k_orig=c(3,4,5,9,10,11,2,6,1,7,8), 
                        #k_order=c(1,2,2,4,5,6,7,8,9,10,11), 
                        k_label = factor(classes, levels=unique(classes))
                        )
#k_label=c('A','B','B','C','D','E','F','G','H','I','J'))

# classcol=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
#            '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
#            '#969696', '#252525')

classcol=c(
  '#000000',
  '#6a3d9a',
  '#a6cee3',
  '#1f78b4',
  '#003366',
  '#33a02c',
  '#b2df8a',
  '#fdbf6f',
  '#ff7f00',
  '#e31a1c')

classr_df_format <- merge(classr_df, reord_k11, by.x='gclass', by.y='k_orig')
classenv <- merge(tzenv[,c('GridID',pred_envar)], classr_df_format, by="GridID")
classenv_melt <- melt(classenv, id.vars=c('GridID','k_label'),value.name='Value')

#df <- classenv_melt
classtableformat <- function(df, tabname) {
  classenv_stats<- setDT(df)[,`:=`(classmean=mean(value, na.rm=T),classsd=sd(value,na.rm=T)), .(variable, k_label)] #Get mean and SD of hydrologic metric for each class
  classenv_stats <- classenv_stats[!duplicated(classenv_stats[,c('variable','k_label')]),]
  classenv_statsmean <- as.data.frame(dcast(classenv_stats, k_label~variable, value.var='classmean'))
  classenv_statssd <- as.data.frame(dcast(classenv_stats, k_label~variable, value.var='classsd'))
  #Format digits for mean and sd
  classenv_stats_meanform <- melt(setDT(digitform(classenv_statsmean,
                                                  cols=2:(ncol(classenv_statsmean)), extradigit=1)),id.vars='k_label',variable.name='Variable', value.name='classmean')
  classenv_stats_sdform <- melt(setDT(digitform(classenv_statssd,
                                                cols=2:(ncol(classenv_statssd)), extradigit=1)),id.vars='k_label',variable.name='Variable', value.name='classsd')
  classenv_format <- merge(classenv_stats_meanform, classenv_stats_sdform, by=c('k_label','Variable'))
  classenv_format$tabcol <- with(classenv_format, paste0(classmean,' (',classsd,')'))
  
  classenv_format[is.na(classenv_format$classsd),'tabcol']  <- '–'
  classenv_format <- merge(as.data.frame(classenv_format[!(classenv_format$Variable=='gclass'),]), pred_envlabel, by.x='Variable', by='var')
  table <- as.data.frame(dcast(setDT(classenv_format), Variable+label~k_label, value.var='tabcol'))
  #table <- as.data.frame(dcast(setDT(classenv_format), k_label~Variable+label, value.var='classmean'))

  kable(table) %>%
    kable_styling(bootstrap_options = "striped", font_size = 10) %>%
    column_spec(3,color=classcol[1]) %>%
    column_spec(4,color=classcol[2]) %>%
    column_spec(5,color=classcol[3]) %>%
    column_spec(6,color=classcol[4]) %>%
    column_spec(7,color=classcol[5]) %>%
    column_spec(8,color=classcol[6]) %>%
    column_spec(9,color=classcol[7]) %>%
    column_spec(10,color=classcol[8]) %>%
    column_spec(11,color=classcol[9]) %>%
    column_spec(12,color=classcol[10]) %>%
    save_kable(tabname, self_contained=T)
}
classtableformat(classenv_melt,tabname=file.path(outdir,'deductive_envgclass_cast_20180704.doc'))

env_labels<-setNames(pred_envarname,pred_envar) #Set metrics labels

setDT(classenv_melt) 
classenv_melt[,.N,.(k_label)]

#Boxplot of metrics for each class
#Remove reach elevation, catchment area, water occurrence and seasonality, temp. coldest quarter, 
excludvar <- c('ReaElvAvg','WsArea','WsWatOcc','WsWatSea','WsBio11Av','gclass')
tz_classenvplot <-ggplot(classenv_melt[!(variable %in% excludvar) & 
                                         !is.na(Value),],
                      aes(x=k_label, y=Value, color=k_label)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~as.factor(variable), scales='free',ncol=4,labeller=as_labeller(env_labels))+
  scale_y_continuous(name='Metric value',expand=c(0.05,0))+
  scale_x_discrete(name='River class')+
  scale_colour_manual(values=classcol) + 
  theme_classic() +
  theme(axis.title = element_text(size=14),
        axis.text.y = element_text(size=11, angle=90),
        axis.text.x = element_text(angle=90, vjust=0.1,  hjust=1.1, color=classcol),
        strip.text = element_text(size=10),
        strip.background = element_rect(color=NA, fill='#f0f0f0'),
        legend.position='none') 
#dir='classo15y_ward_rawsub3'
#outdirclass <- file.path(outdir,dir)
png(file.path(outdir,'11class_boxplot_20201002.png'),width = 8.5, height=10,units='in',res=600)
print(tz_classenvplot)
dev.off()
