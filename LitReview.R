# Hydrologic Classification Literature Review - Summary Figure Code (04/03/18)

setwd("~/Dropbox/DOCUMENTS/Projects/Contracts/Tanzania/Literature/Summary Figures")

hydro <- read.csv('litreview.csv', header=TRUE)
hydro[hydro==""]  <- NA 

#subset the data for just hydrological classifications
hydro.sub <- subset(hydro, Type == "Hydrological" | Type == "Hydrogeomorphological") 
summary(hydro.sub)

summary(hydro)
names(hydro)

library(ggplot2)
library(vegan)
library("gridExtra", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
library("topicmodels", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")

# FIGURE
# plot for the # studies over time

cumsum(hydro.sub$Year)
summ<-data.frame(table(hydro.sub$Year))
summary(summ)
summ$Var1<-as.numeric(summ$Var1)

# cumulative number of classifications over time
ggplot(hydro.sub, aes(x=Year,color='red',size=3)) + stat_bin(aes(y=cumsum(..count..)),geom="step") + ylab("# of publications") + xlab("Year") + scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks=seq(1960,2018,by=4)) +
  theme (
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(colour="black",hjust=1,size=16,angle=45), 
    axis.text.y = element_text(colour="black",size=16),
    axis.title.y=element_text(colour="black",size=16), 
    axis.title.x=element_text(colour="black",size=16),
    legend.position="none"
  )

# Deductive vs. Inductive
a<-ggplot(hydro.sub, aes(Approach1)) + geom_bar(aes(y=(..count..)/sum(..count..)*100),fill='darkblue') + ylab("% studies") +
  scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0, 100, by=10)) + ggtitle("Approach") +
  theme (
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="black",size=12), 
    axis.text.y = element_text(colour="black",size=12),
    axis.title.y=element_text(colour="black",vjust=1.5, size=12), 
    axis.title.x=element_blank(),
    plot.title=element_text(hjust=-0.12,vjust=2.5, size=18),
    legend.position="none"
  )

# Membership
#hydro.sub<-hydro.sub[!is.na(hydro.sub$Approach2),] # omitting rows with NAs
b<-ggplot(hydro.sub, aes(Membership)) + geom_bar(aes(y=(..count..)/sum(..count..)*100),fill='darkblue') + ylab("% studies") +
  scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0, 100, by=10)) + ggtitle("Membership") +
  theme (
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="black",size=12), 
    axis.text.y = element_text(colour="black",size=12),
    axis.title.y=element_text(colour="black",vjust=1.5, size=12), 
    axis.title.x=element_blank(),
    plot.title=element_text(hjust=-0.12,vjust=2.5, size=18),
    legend.position="none"
  )

# HydroDataType
hydro.sub2<-hydro.sub[!is.na(hydro.sub$HydroDataType),] # omitting rows with NAs
c<-ggplot(hydro.sub2, aes(HydroDataType)) + geom_bar(aes(y=(..count..)/sum(..count..)*100),fill='darkblue') + ylab("% studies") +
  scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0, 100, by=10)) + ggtitle("Hydrologic Data") +
  theme (
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="black",size=12), 
    axis.text.y = element_text(colour="black",size=12),
    axis.title.y=element_text(colour="black",vjust=1.5, size=12), 
    axis.title.x=element_blank(),
    plot.title=element_text(hjust=-0.12,vjust=2.5, size=18),
    legend.position="none"
  )

# TemporalGrain
hydro.sub2<-hydro.sub[!is.na(hydro.sub$TemporalGrain),] # omitting rows with NAs
d<-ggplot(hydro.sub2, aes(TemporalGrain)) + geom_bar(aes(y=(..count..)/sum(..count..)*100),fill='darkblue') + ylab("% studies") +
  scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0, 100, by=10)) + ggtitle("Temporal Grain") +
  theme (
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="black",size=12), 
    axis.text.y = element_text(colour="black",size=12),
    axis.title.y=element_text(colour="black",vjust=1.5, size=12), 
    axis.title.x=element_blank(),
    plot.title=element_text(hjust=-0.12,vjust=2.5, size=18),
    legend.position="none"
  )

# SpatialGrain
#hydro.sub2<-hydro.sub[!is.na(hydro.sub$HydroDataType),] # omitting rows with NAs
#hydro.sub<-factor(hydro.sub$SpatialGrain,c("Pixel","Reach","Catchment","Region"))
e<-ggplot(hydro.sub, aes(SpatialGrain)) + geom_bar(aes(y=(..count..)/sum(..count..)*100),fill='darkblue') + ylab("% studies") +
  scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0, 100, by=10)) + ggtitle("Spatial Grain") +
  theme (
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="black",size=12), 
    axis.text.y = element_text(colour="black",size=12),
    axis.title.y=element_text(colour="black",vjust=1.5, size=12), 
    axis.title.x=element_blank(),
    plot.title=element_text(hjust=-0.12,vjust=2.5, size=18),
    legend.position="none"
  )


# Prediction
hydro.sub2<-hydro.sub[!is.na(hydro.sub$Prediction),] # omitting rows with NAs
f<-ggplot(hydro.sub2, aes(Prediction)) + geom_bar(aes(y=(..count..)/sum(..count..)*100),fill='darkblue') + ylab("% studies") +
  scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0, 100, by=10)) + ggtitle("Prediction") +
  theme (
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="black",size=12), 
    axis.text.y = element_text(colour="black",size=12),
    axis.title.y=element_text(colour="black",vjust=1.5, size=12), 
    axis.title.x=element_blank(),
    plot.title=element_text(hjust=-0.12,vjust=2.5, size=18),
    legend.position="none"
  )

# SpatialExtent
g<-ggplot(hydro.sub, aes(SpatialExtent)) + geom_bar(aes(y=(..count..)/sum(..count..)*100),fill='darkblue') + ylab("% studies") +
  scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0, 100, by=10)) + ggtitle("Spatial Extent") +
  theme (
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="black",size=12), 
    axis.text.y = element_text(colour="black",size=12),
    axis.title.y=element_text(colour="black",vjust=1.5, size=12), 
    axis.title.x=element_blank(),
    plot.title=element_text(hjust=-0.12,vjust=2.5, size=18),
    legend.position="none"
  )


# FIGURE
grid.arrange(a,b,c,d,e,f, ncol=2,nrow=3)
