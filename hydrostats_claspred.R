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
library(ggdendro)
source(file.path(rootdir,"bin/outside_src/Biostats.R"))
source(file.path(rootdir,"bin/outside_src/Flowscreen.hyear.internal.R"))

rootdir="F:/Tanzania/Tanzania" #UPDATE
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
rufidat_impute <- read.csv(file.path('rufiji_hydrodataimpute', 'rufidat_interp.csv'), colClasses=c('Date',rep('numeric',34)))
colnames(rufidat_impute)[2:(ncol(rufidat_impute))] <- substr(colnames(rufidat_impute),2,10)[2:(ncol(rufidat_impute))]
rufidat_gapsummary <- read.csv(file.path(datadir, 'rufidat_gapsummary.csv'))
rufidat_post1991<-read.csv(file.path(datadir, 'gageselect_post1991comp90.csv'))
rufidat_o15y<-read.csv(file.path(datadir, 'gageselect_o15comp90.csv'))

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
predsmelt <- merge(predsmelt, rufidat_o15y, by='ID',all.x=T)

rufidat_select <- predsmelt[predsmelt$gap_per<=0.1 & (predsmelt$ycount_o15>=15 | predsmelt$ycount1991>=10) & predsmelt$hyear<2017,]
str(rufidat_select)

#####################################################################################
# Compute hydrologic metrics
allHITcomp <- function(dfhydro, dfenv, gageID, templateID='1KA9',hstats="all", floodquantile=0.95) {
  #Get template
  HITall_template <- singleHITcomp(dfhydro, dfenv, templateID)
  colnames(HITall_template)[2] <- templateID
  HITall <- data.frame(indice=HITall_template$indice) 
  #Compute metrics for all gages
  for (gage in unique(dfhydro[,gageID])) {
    print(gage)
    try({
      #Check data for completeness
      dailyQClean <- validate_data(dfhydro[dfhydro$ID==gageID,c("Date", "Flow")], yearType="water")
      #Calculate all hit stats
      calc_allHITout <- calc_allHIT(dailyQClean, yearType="water", stats=hstats, digits=10, pref="mean",
                                    drainArea=dfenv[dfenv$RGS_No==gageID,'WsArea'], floodThreshold = quantile(dailyQClean$discharge, floodquantile))
      colnames(calc_allHITout)[2] <- gage
      HITall <- merge(HITall, calc_allHITout, by='indice')
    })
  }
  HITall_formatmelt <-melt(setDT(HITall), id.vars = "indice",variable.name = "ID") 
  HITall_formatmelt[is.infinite(HITall_formatmelt$value),'value'] <- NA
  HITall_formatmelt[is.nan(HITall_formatmelt$value),'value'] <- NA
  return(HITall_formatmelt)
}
HITpost1991 <- allHITcomp(rufidat_select, gagesenv, 'ID')

#Calculate mag7 stats
#magnifStatsOut <- calc_magnifSeven(dailyQClean,yearType="water",digits=3)

#######################################################################################
# Box plot of metrics
HITallbox<- HITall_formatmelt
HITallbox$group1 <- as.factor(substr(HITall_formatmelt$indice,1,1))
HITallbox$group2 <- as.factor(substr(HITall_formatmelt$indice,2,2))
HITallbox$indice_sub <- substr(HITall_formatmelt$indice,3,5)
HITallbox$indice_sub <- factor(HITallbox$indice_sub, levels = unique(HITallbox$indice_sub[order(as.numeric(as.character(HITallbox$indice_sub)))]))
HITallbox$group1 <- factor(HITallbox$group1, levels = c('m','f','d','t','r'), 
                           labels = c("Magnitude-m", "Frequency-f", "Duration-d",'Timing-t',"Rate of change-r"))
HITallbox$group2 <- factor(HITallbox$group2, levels = c('h','a','l'), labels=c("High flow-h",'Average flow-a',"Low flow-l"))
lab <- ddply(HITallbox, .(indice), summarize, 
             labels=median(value)-1.58*(quantile(value,0.75,na.rm=T)-quantile(value,0.25,na.rm=T))/sqrt(length(value)),
             labels2=min(value))
HITallbox <-  merge(HITallbox,lab, by='indice')

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
png(file.path(outdir,'HITallboxplot.png'),width=20, height=12,units='in',res=300)
HITallboxplot
dev.off()

#HITallbox[HITallbox$indice_sub %in% seq(1,25,2),'label'] <- as.character(HITallbox[HITallbox$indice_sub %in% seq(1,25,2),'indice_sub'])
#HITallbox[!(HITallbox$indice_sub %in% seq(1,25,2)),'label'] <- ''
#scale_x_discrete(name='Metric number (Appendix 1)',aes(breaks=indice_sub),labels = HITallbox$label)

#####################Format data to use in classification and compute Gower's distance############
HITdist <- function(HITdf, gagelist) {
  HITdfsel <- HITdf[HITdf$ID %in% gagelist,]
  HITdf_format <- dcast(HITdfsel, ID ~ indice)
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
gaugegow1991 <- HITdist(HITall_formatmelt, unique(rufidat_post1991[rufidat_post1991$ycount1991>=10,'ID']))
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
colnames(pcoa_scores) <- paste("PC",seq(1,ncol(pcoa_scores)),sep="")
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
outdirclass <- file.path(outdir,'class1991_ward_pcoa')
if (dir.exists(outdirclass )) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdirclass))
  dir.create(outdirclass )
}
classpc5 <-cutree(gaugecla_pcward, k=5)
write.csv(classpc5, file.path(outdirclass,'class_pcoagow_ward_5.csv'))
classpc6 <-cutree(gaugecla_pcward, k=6)
classpc6_df <- data.frame(ID=names(classpc6), gclass=classpc6) 
write.csv(gauge_classpc6, file.path(outdirclass,'class_pcoagow_ward_6.csv'))
#Boxplots
HITpcoa_new <- cbind(gauge_class,pcoa_scores)
#box.plots(HITpcoa_new, by='gauge_class')

#Check cluster similarity
cluster_similarity(classr5, classpc5, similarity='rand',method='independence')
comembership_table(classr5, classpc5)
########################################Predict based on raw-hydro metrics classification and raw environmental predictors############################
pred_envar <-c('WsArea','CatSloAvg','CatWatExt','CatWatOcc','CatWatSea','CatDRocAvg','CatPopDen','ReaElvAvg','ReaSloAvg','WsLakInd','WsBio01Av','WsBio07Av','WsBio12Av','LCSum_45')
gagesenv_r5 <- merge(gagesenvrec[gagesenvrec$RGS_No %in% rufidat_select$ID, c('RGS_No',pred_envar)],classr5_df, by.x='RGS_No', by.y='ID')
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

########################################Predict based on PCoA classification and raw environmental predictors############################
pred_envar <-c('WsArea','CatSloAvg','CatWatExt','CatWatOcc','CatWatSea','CatDRocAvg','CatPopDen','ReaElvAvg','ReaSloAvg','WsLakInd','WsBio01Av','WsBio07Av','WsBio12Av','LCSum_45')
gagesenv_pc6 <- merge(gagesenvrec[gagesenvrec$RGS_No %in% rufidat_select$ID, c('RGS_No',pred_envar)],classpc6_df, by.x='RGS_No', by.y='ID')
rownames(gagesenv_pc6) <- gagesenv_pc6$RGS_No
gagesenv_pc6 <- gagesenv_pc6[,-which('RGS_No' %in% colnames(gagesenv_pc6))]
gagesenv_pc6$gclass <- as.factor(gagesenv_pc6$gclass)
rufienvsel_pc6r <- rufienv[,c('GridID',pred_envar)]
rownames(rufienvsel_pc6r) <- rufienvsel_pc6r$GridID
rufienvsel_pc6r <- rufienvsel_pc6r[,-which('GridID' %in% colnames(rufienvsel_pc6r))]
#Single tree
cat.pc6r <- rpart(gclass~., data=gagesenv_pc6, method='class',control=rpart.control(minsplit=1, minbucket=1, cp=0.1))
summary(cat.pc6r)
rpart.plot(cat.pc6r)
#Boosted tree
adaboost.pc6r <- boosting(gclass~., data=gagesenv_pc6, boos=TRUE, mfinal=100,  control=rpart.control(minsplit=1, minbucket=1, cp=0.1))
barplot(adaboost.pc6r$imp[order(adaboost.pc6r$imp, decreasing = TRUE)],
        ylim = c(0, 100), main = "Variables Relative Importance",
        col = "lightblue")
adaboostcv.pc6r <- boosting.cv(gclass~., data=gagesenv_pc6, boos=TRUE, mfinal=10,  control=rpart.control(minsplit=1, minbucket=1, cp=0.1), v=19)
adaboostcv.pc6r
#Predict
rufi_pc6r<- predict.boosting(adaboost.pc6r, newdata=rufienvsel_pc6r, newmfinal=length(adaboost.pc6r$trees))
#rufi_pc6rmaxprob <- adply(rufi_pc6r$prob, 1, max)
#qplot(rufi_pc6rmaxprob$V1)
rufi_pc6r_pred <- cbind(rufienvsel_pc6r,gclass=rufi_pc6r$class)
rufi_pc6r_pred$GridID <- as.integer(rownames(rufi_pc6r_pred)) 
write.dbf(rufi_pc6r_pred[,c('GridID','gclass')], file.path(outdir, "class_ward_pcoa/predict_pc6r.dbf"))













######################################### CLASSIFICATION BASED ON ENTIRE PERIOD > 15 YEARS OF DATA ################################
########################################Classify based on raw indices############################################
gaugegow_o15y <- HITdist(HITall_formatmelt, unique(rufidat_o15y[rufidat_o15y$ycount_o15>=15,'ID']))
gaugecla_ward <-hclust(gaugegow_o15y, method='ward.D') #Classify using hierarchical agglomerative using Ward's minimum variance method

#Diagnostics
hclus.table(gaugecla_ward)
plot(gaugecla_ward, main="Ward's distance gauge cluster dendogram", xlab='Gauge ID', ylab="Gower's distance", hang=-1)   #par(mfrow=c(1,1))
coef.hclust(gaugecla_ward) #Compute agglomerative coefficient
cor(gaugegow_o15y, cophenetic(gaugecla_ward)) #Compute cophenetic coefficient
hclus.cophenetic(gaugegow_o15y, gaugecla_ward) #Plot cophenetic relationship 
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
