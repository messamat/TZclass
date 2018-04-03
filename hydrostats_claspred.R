#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12
#Subtask B-1 and B-2

#Authors: Mathis L. Messager and Dr. Julian D. Olden
#Contact info: messamat@uw.edu
#Date created: 03/30/2018
#Date last updated: 04/02/2018

#Purpose: compute hydrological metrics and classify stream gauges

library(plyr)
library(dplyr)
library(hydrostats)
library(data.table)
devtools::install_github("messamat/EflowStats") #Intro Vignette: https://cdn.rawgit.com/USGS-R/EflowStats/9507f714/inst/doc/intro.html
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
library(ggplot2)
library(ggdendro)

rootdir="F:/Tanzania/Tanzania" #UPDATE
setwd(file.path(rootdir,"results")) 
datadir = file.path(getwd(),paste('rufiji_hydrodatafilter','20180327',sep='_')) #UPDATE
origdatadir = file.path(rootdir,"data") 
outdir=file.path(getwd(),'rufiji_hydrodatastats_20180331')
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}

rufidat_clean <- read.csv(file.path('rufiji_hydrodatainspect_20180326','rufidat_clean.csv'), colClasses=c('factor','Date','numeric','character','character'))
rufidat_impute <- read.csv(file.path('rufiji_hydrodataimpute_20180329', 'rufidat_interp.csv'), colClasses=c('Date',rep('numeric',34)))
colnames(rufidat_impute)[2:(ncol(rufidat_impute))] <- substr(colnames(rufidat_impute),2,10)[2:(ncol(rufidat_impute))]
rufidat_gapsummary <- read.csv(file.path(datadir, 'rufidat_gapsummary.csv'))
rufidat_post1991<-read.csv(file.path(datadir, 'gageselect_post1991comp90.csv'))

gagesenv <- read.dbf(file.path(getwd(),'gages_netjoinclean.dbf'))
gagesenvrec <- merge(gagesenv, unique(rufidat_clean[,c('ID','SYM')]), by.x='RGS_No', by.y='ID', all.x=F)
rufienv <- read.dbf(file.path(getwd(),'streamnet118_rufiji_finaltabclean.dbf'))

predsmelt <-melt(setDT(rufidat_impute),id.vars = 'Date',value.name='Flow',variable.name='ID')
predsmelt <- predsmelt[,c(2,1,3)]
predsmelt$year <- as.numeric(format(predsmelt$Date, "%Y"))
predsmelt$month <- as.numeric(format(predsmelt$Date, "%m"))
predsmelt<-hyear.internal(predsmelt,hyrstart=10) #Ignore hdoy
predsmelt <- merge(predsmelt, rufidat_gapsummary, by=c('ID','hyear'),all.x=T)
predsmelt <- merge(predsmelt, rufidat_post1991, by='ID',all.x=T)

rufidat_select <- predsmelt[predsmelt$hyear>=1991 & predsmelt$hyear<2017 & predsmelt$ycount>=10 & predsmelt$gap_per<=0.1,]
str(rufidat_select)

#####################################################################################
# Compute hydrologic metrics
HITcomp <- function(dfhydro, dfenv, gageID, hstats="all", floodquantile=0.95) {
  #Check data for completeness
  dailyQClean <- validate_data(dfhydro[dfhydro$ID==gageID,c("Date", "Flow")], yearType="water")
  #Calculate all hit stats
  calc_allHITout <- calc_allHIT(dailyQClean, yearType="water", stats=hstats, digits=10, pref="mean",
                                drainArea=dfenv[dfenv$RGS_No==gageID,'WsArea'], floodThreshold = quantile(dailyQClean$discharge, floodquantile))
  return(calc_allHITout)
}
#Get template
HITall_template <- HITcomp(rufidat_select, gagesenvrec, '1KA9')
colnames(HITall_template)[2] <- '1KA9'
HITall <- data.frame(indice=HITall_template$indice) 
#Compute metrics for all gages
for (gage in unique(rufidat_select$ID)) {
  print(gage)
  try({
    calc_allHITout <- HITcomp(rufidat_select, gagesenvrec, gage)
    colnames(calc_allHITout)[2] <- gage
    HITall <- merge(HITall, calc_allHITout, by='indice')
  })
}

#Calculate mag7 stats
#magnifStatsOut <- calc_magnifSeven(dailyQClean,yearType="water",digits=3)

######################################################################
#Format data to use in classification
HITall_format <-melt(setDT(HITall), id.vars = "indice",variable.name = "ID") 
HITall_format[is.infinite(HITall_format$value),'Value'] <- NA
HITall_format[is.nan(HITall_format$value),'Value'] <- NA
HITall_format <- dcast(HITall_format, ID ~ indice)
HITall_format <- merge(HITall_format, gagesenvrec[,c('RGS_No','WsArea')], by.x='ID', by.y='RGS_No')
dimindices <- c('ma1','ma2',paste('ma',seq(12,23),sep=''),paste('ml',seq(1,12),sep=''),paste('mh',seq(1,12),sep=''), 
                paste('dl',seq(1,5),sep=''),paste('dh',seq(1,5),sep=''),'ra1','ra3','ra6','ra7') #List of dimensional indices from Kennen et al. 2007
HITall_format <- as.data.frame(HITall_format[,(dimindices) := lapply(.SD, function(x) round(x/WsArea, digits=10)), .SDcols=dimindices]) #Standardize dimensional indices by drainage area
row.names(HITall_format) <- HITall_format$ID
HITall_format <- HITall_format[,-which(colnames(HITall_format) %in% c('ID','WsArea'))] #Get rid of non-indices columns
HITall_stand <- data.stand(HITall_format[,2:(ncol(HITall_format))],method='standardize',margin='column',plot=F) #z-standardize columnwise 
statHIT<- summary(HITall_stand)
gauge_gow<- gowdis(HITall_stand, w=rep(1,ncol(HITall_stand)), asym.bin = NULL) #Compute Gower's distance so that missing values will not be taken in account

########################################Classify based on raw indices############################################
gaugecla_ward <-hclust(gauge_gow, method='ward.D') #Classify using hierarchical agglomerative using Ward's minimum variance method

#Diagnostics
hclus.table(gaugecla_ward)
plot(gaugecla_ward, main="Ward's distance gauge cluster dendogram", xlab='Gauge ID', ylab="Gower's distance", hang=-1)   #par(mfrow=c(1,1))
coef.hclust(gaugecla_ward) #Compute agglomerative coefficient
cor(gauge_gow, cophenetic(gaugecla_ward)) #Compute cophenetic coefficient
hclus.cophenetic(gauge_gow, gaugecla_ward) #Plot cophenetic relationship 
hclus.scree(gaugecla_ward) #Check out scree plot
plot(gaugecla_ward, main="Ward's distance gauge cluster dendogram", xlab='Gauge ID', ylab="Gower's distance", hang=-1)   #par(mfrow=c(1,1))
rect.hclust(gaugecla_ward, k=5) #Draw rectangle around k classes
#Test significance of classes
clus.stab <- pvclust(t(HITall_stand), method.hclust='ward.D', method.dist='cor',use.cor="pairwise.complete.obs", nboot=4999)
plot(clus.stab)
pvrect(clus.stab, alpha=0.90)
#Get gauge classes
gauge_class <-cutree(gaugecla_ward, k=5)
write.csv(gauge_class, file.path(outdir,'class_rawgow_ward_5.csv'))
#Boxplots
HITall_new <- cbind(gauge_class,HITall_format)
#box.plots(HITall_new, by='gauge_class')

########################################Classify based on PCoA synthetic indices######################################
###PCOA
#Convert gowdis output to matrix
gauge_gowmat <- as.matrix(gauge_gow)
#Run PCoA without a constant added using the stats package
gauge_pcoa <- cmdscale(gauge_gowmat, k = 5, eig = T, add = F)
#Diagnostic PCOA
eig1 <-data.frame(x=1:length(gauge_pcoa$eig),y=100*gauge_pcoa$eig/sum(gauge_pcoa$eig))
stick1 <- data.frame(x=1:length(gauge_pcoa$eig),y=100*bstick(length(gauge_pcoa$eig)))
ggplot(eig1, aes(x, y)) + geom_point(color='red') + geom_point(data=stick1)
sum(eig1[1:5,2])
#Calculate PC loadings using correlation analysis 
vec_ga <- envfit(scores(gauge_pcoa),  HITall_format, perm = 1000, na.rm = T)
vec_ga
vec <- data.frame(indice=colnames(HITall_format))
vec <- cbind(vec, vec_ga$vectors$arrows)
vec$pval <- vec_ga$vectors$pvals
#Visualize scores
#ordiplot(gauge_pcoa, choices=c(1,2), type="text",display='sites', xlab='PCo 1 (48%)', ylab='PCo 2 (24%)')
pcoa_scores <- as.data.frame(gauge_pcoa$points)
colnames(pcoa_scores) <- c("PC1", "PC2")
ggplot(pcoa_scores, aes(x = PC1, y = PC2, label = rownames(pcoa_scores))) + geom_label() +
  geom_segment(data=vec[vec$pval<=0.1,], aes(x = rep(0, 36), y = rep(0, 36), xend = Dim1 , yend = Dim2, label=indice)) + 
  geom_text(data=vec[vec$pval<=0.1,], aes(x = Dim1, y = Dim2, label = indice),  color = "red") + 
  theme_classic() +
  labs(title='PCoA gauges', x='PCo 1 (48%)', y='PCo2 (24%)')

###Classification
gauge_pceuc <- vegdist(pcoa_scores, method='euclidean')
gaugecla_pcward <-hclust(gauge_pceuc, method='ward.D') #Classify using hierarchical agglomerative using Ward's minimum variance method

#Diagnostics
hclus.table(gaugecla_pcward)
plot(gaugecla_pcward, main="Ward's distance gauge cluster dendogram", xlab='Gauge ID', ylab="Gower's distance", hang=-1)   #par(mfrow=c(1,1))
coef.hclust(gaugecla_pcward) #Compute agglomerative coefficient
cor(gauge_pceuc, cophenetic(gaugecla_pcward)) #Compute cophenetic coefficient
hclus.cophenetic(gauge_pceuc, gaugecla_pcward) #Plot cophenetic relationship 
hclus.scree(gaugecla_pcward) #Check out scree plot
plot(gaugecla_pcward, main="Ward's distance gauge cluster dendogram", xlab='Gauge ID', ylab="Gower's distance", hang=-1)   #par(mfrow=c(1,1))
rect.hclust(gaugecla_pcward, k=3) #Draw rectangle around k classes
rect.hclust(gaugecla_pcward, k=6) #Draw rectangle around k classes
#Test significance of classes
clus.stab <- pvclust(t(pcoa_scores), method.hclust='ward.D', method.dist='cor',use.cor="pairwise.complete.obs", nboot=1999)
plot(clus.stab)
pvrect(clus.stab, alpha=0.90)
#Get gauge classes
gauge_classpc5 <-cutree(gaugecla_pcward, k=5)
write.csv(gauge_classpc5, file.path(outdir,'class_pcoagow_ward_5.csv'))
gauge_classpc6 <-cutree(gaugecla_pcward, k=6)
write.csv(gauge_classpc6, file.path(outdir,'class_pcoagow_ward_6.csv'))
#Boxplots
HITpcoa_new <- cbind(gauge_class,pcoa_scores)
#box.plots(HITpcoa_new, by='gauge_class')

#Check cluster similarity
cluster_similarity(gauge_class, gauge_classpc5, similarity='rand',method='independence')
comembership_table(gauge_class, gauge_classpc5)

########################################Predict based on PCoA classification and raw environmental predictors############################
pred_envar <-c('WsArea','CatSloAvg','CatWatExt','CatWatOcc','CatWatSea','CatDRocAvg','CatPopDen','ReaElvAvg','ReaSloAvg','WsLakInd','WsBio01Av','WsBio07Av','WsBio12Av','LCSum_45')
gagesenv_pc6 <- merge(gagesenvrec[gagesenvrec$RGS_No %in% rufidat_select$ID, c('RGS_no',pred_envar)],gauge_classpc6, by.x='RGS_No', by.y='ID')
adaboost.pc6r <- boosting(gauge_class~#pcoa 6-classes with raw environmental variables
which(pred_envar <-c('WsArea','CatSloAvg','CatWatExt','CatWatOcc','CatWatSea','CatDRocAvg','CatPopDen','ReaElvAvg','ReaSloAvg','WsLakInd','WsBio01Av','WsBio07Av','WsBio12Av','LCSum_45') %in% colnames(gagesenvrec))















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
