#USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#Contract No. EDH-I-00-08-00023-00 621-TO 12

#Authors: Dr. Julian D. Olden | Mathis L. Messager 
#Contact info: olden@uw.edu | messamat@uw.edu

#Purpose: Summarize literature review on river classifications

library(ggplot2)
library(vegan)
library(gridExtra)
library(topicmodels)

rootdir="F:/Tanzania/Tanzania"
setwd(file.path(rootdir,"results")) #UPDATE
datadir = file.path(rootdir,'data/LitReview') #UPDATE
outdir=file.path(getwd(),'litreview')
if (dir.exists(outdir)) {
  print('Directory already exists')
} else {
  print(paste('Create new directory:',outdir))
  dir.create(outdir)
}


hydro <- read.csv(file.path(datadir,'litreview.csv'))
hydro[hydro==""]  <- NA 

#subset the data for just hydrological classifications
hydro.sub <- subset(hydro, Type == "Hydrological" | Type == "Hydrogeomorphological") 
summary(hydro.sub)

summary(hydro)
names(hydro)

# FIGURE
# plot for the # studies over time

cumsum(hydro.sub$Year)
summ<-data.frame(table(hydro.sub$Year))
summary(summ)
summ$Var1<-as.numeric(summ$Var1)

# cumulative number of classifications over time
cumnum <- ggplot(hydro.sub, aes(x=Year,color='red',size=2)) + 
  stat_bin(aes(y=cumsum(..count..)),geom="step") + 
  ylab("Cumulative number of publications") + xlab("Year") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks=seq(1960,2018,by=4)) +
  theme (
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(colour="black",hjust=1,size=12,angle=45), 
    axis.text.y = element_text(colour="black",size=11),
    axis.title.y=element_text(colour="black",size=11), 
    axis.title.x=element_text(colour="black",size=11),
    legend.position="none"
  )
png(file.path(outdir,'cumulutative_studies.png'),width = 6.5, height=4.5,units='in',res=300)
print(cumnum)
dev.off()


# Deductive vs. Inductive
hydro.sub$Approach1 <- factor(hydro.sub$Approach1, levels = c('Inductive','Deductive','Mixed'))
a<-ggplot(hydro.sub, aes(Approach1)) + geom_bar(aes(y=(..count..)/sum(..count..)*100),fill='darkblue') + ylab("% studies") +
  scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0, 100, by=10)) + ggtitle("(a) Approach") +
  theme (
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="black",size=11), 
    axis.text.y = element_text(colour="black",size=11),
    axis.title.y=element_text(colour="black",vjust=1.5, size=12), 
    axis.title.x=element_blank(),
    plot.title=element_text(hjust=-0.12,vjust=2.5, size=14),
    legend.position="none",
    plot.margin = unit(c(0.25, 0, 0.25, 0), "cm")
  )

# Membership
#hydro.sub<-hydro.sub[!is.na(hydro.sub$Approach2),] # omitting rows with NAs
b<-ggplot(hydro.sub, aes(Membership)) + geom_bar(aes(y=(..count..)/sum(..count..)*100),fill='darkblue') + ylab("% studies") +
  scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0, 100, by=10)) + ggtitle("(b) Membership") +
  theme (
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="black",size=11), 
    axis.text.y = element_text(colour="black",size=11),
    axis.title.x=element_blank(),
    plot.title=element_text(hjust=-0.12,vjust=2.5, size=14),
    legend.position="none",
    plot.margin = unit(c(0.25, 0, 0.25, 0), "cm"),
    axis.title.y= element_blank()
  )

# HydroDataType
hydro.sub2<-hydro.sub[!is.na(hydro.sub$HydroDataType),] # omitting rows with NAs
hydro.sub$HydroDataType <- factor(hydro.sub$HydroDataType , levels = c('Observed','Modeled','Obs/Mod'), labels = c('Observed', 'Modeled', 'Mixed'))
c<-ggplot(hydro.sub2, aes(HydroDataType)) + geom_bar(aes(y=(..count..)/sum(..count..)*100),fill='darkblue') + ylab("% studies") +
  scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0, 100, by=10)) + ggtitle("(c) Hydrologic Data") +
  theme (
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="black",size=11), 
    axis.text.y = element_text(colour="black",size=11),
    axis.title.y=element_text(colour="black",vjust=1.5, size=12), 
    axis.title.x=element_blank(),
    plot.title=element_text(hjust=-0.12,vjust=2.5, size=14),
    plot.margin = unit(c(0.25, 0, 0.25, 0), "cm"),
    legend.position="none"
  )

# TemporalGrain
hydro.sub2<-hydro.sub[!is.na(hydro.sub$TemporalGrain),] # omitting rows with NAs
hydro.sub$TemporalGrain <- factor(hydro.sub$TemporalGrain , levels = c('Annual','Monthly','Daily','Hourly'))
d<-ggplot(hydro.sub2, aes(TemporalGrain)) + geom_bar(aes(y=(..count..)/sum(..count..)*100),fill='darkblue') + ylab("% studies") +
  scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0, 100, by=10)) + ggtitle("(d) Temporal Grain") +
  theme (
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="black",size=11), 
    axis.title.x=element_blank(),
    axis.title.y= element_blank(),
    plot.title=element_text(hjust=-0.12,vjust=2.5, size=14),
    plot.margin = unit(c(0.25, 0, 0.25, 0), "cm"),
    legend.position="none"
  )

# SpatialGrain
#hydro.sub2<-hydro.sub[!is.na(hydro.sub$HydroDataType),] # omitting rows with NAs
hydro.sub$SpatialGrain <-factor(hydro.sub$SpatialGrain,levels=c("Pixel","Reach","Catchment","Region"))
e<-ggplot(hydro.sub, aes(SpatialGrain)) + geom_bar(aes(y=(..count..)/sum(..count..)*100),fill='darkblue') + ylab("% studies") +
  scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0, 100, by=10)) + ggtitle("(e) Spatial Grain") +
  theme (
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="black",size=11), 
    axis.text.y = element_text(colour="black",size=11),
    axis.title.y=element_text(colour="black",vjust=1.5, size=12), 
    axis.title.x=element_blank(),
    plot.title=element_text(hjust=-0.12,vjust=2.5, size=14),
    plot.margin = unit(c(0.25, 0, 0.25, 0), "cm"),
    legend.position="none"
  )


# Prediction
hydro.sub2<-hydro.sub[!is.na(hydro.sub$Prediction),] # omitting rows with NAs
hydro.sub2$Prediction<- factor(hydro.sub2$Prediction , levels = c('None','ClassFirst','PredFirst'), labels=c('None','Classify-first', 'Predict-first'))
f<-ggplot(hydro.sub2, aes(Prediction)) + 
  geom_bar(aes(y=(..count..)/sum(..count..)*100),fill='darkblue') + ylab("% studies") +
  scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0, 100, by=10)) + ggtitle("(f) Prediction") +
  theme (
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="black",size=11), 
    axis.title.x=element_blank(),
    axis.title.y= element_blank(),
    plot.title=element_text(hjust=-0.12,vjust=2.5, size=14),
    plot.margin = unit(c(0.25, 0, 0.25, 0), "cm"),
    legend.position="none"
  )

# SpatialExtent
g<-ggplot(hydro.sub, aes(SpatifalExtent)) + geom_bar(aes(y=(..count..)/sum(..count..)*100),fill='darkblue') + ylab("% studies") +
  scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0, 100, by=10)) + ggtitle("Spatial Extent") +
  theme (
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="black",size=11), 
    axis.text.y = element_text(colour="black",size=11),
    axis.title.y=element_text(colour="black",vjust=1.5, size=12), 
    axis.title.x=element_blank(),
    plot.title=element_text(hjust=-0.12,vjust=2.5, size=14),
    plot.margin = unit(c(0.25, 0, 0.25, 0), "cm"),
    legend.position="none"
  )


# FIGURE
png(file.path(outdir,'studies_characs.png'),width = 6.5, height=7.5,units='in',res=300)
grid.draw(gtable_cbind(gtable_rbind(ggplotGrob(a), ggplotGrob(c), ggplotGrob(e)),
                       gtable_rbind(ggplotGrob(b), ggplotGrob(d), ggplotGrob(f))))
dev.off()

