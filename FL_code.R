renv::init(repos = "https://eur01.safelinks.protection.outlook.com/?url=https%3A%2F%2Fpackagemanager.posit.co%2Fcran%2F2023-10-03&data=05%7C02%7Cxin.ding.16%40ucl.ac.uk%7C5396216764284499ccea08dc545ef112%7C1faf88fea9984c5b93c9210a11d9a5c2%7C0%7C0%7C638478013958602799%7CUnknown%7CTWFpbGZsb3d8eyJWIjoiMC4wLjAwMDAiLCJQIjoiV2luMzIiLCJBTiI6Ik1haWwiLCJXVCI6Mn0%3D%7C0%7C%7C%7C&sdata=Xvpcc0e%2BU0jse2oF%2FddXfsJ2wtiPisd1oIJbyFDwwFA%3D&reserved=0")

library(maptools)
library(spatstat)
library(colorRamps)
library(rgdal)
library(rgeos)
library(raster)
library(circular)
library(corrplot)
library(Hmisc)
library(spdep)
library(ggbiplot)
library(ggplot2)
library(ggridges)


##Read data
setwd("/Users/XinDing/r")
alldata <- read.csv(file="FL_data/FL_basic data.csv", header=TRUE, na="NA", fileEncoding="UTF-8")
T1 <- alldata[alldata$Trench == "T1", ]
T2 <- alldata[alldata$Trench == "T2", ]
T3 <- alldata[alldata$Trench == "T3", ]
TOK <- alldata[alldata$Trench == "TOK", ]

##Data preparation
T1F <- T1[T1$Gtype == "B",]
T1L <- T1[T1$Gtype == "S",]
T1XY <- readOGR(dsn="FL_outline", layer="FL_T1_outline")
T1XYmainowin <- as.owin(T1XY)
T1XZmainowin <- owin(xrange=c(808.0583418, 808.063533), yrange=c(0.934456, 0.935035))
T1YZmainowin <- owin(yrange=c(629.7346941, 629.740503), xrange=c(0.934456, 0.935035))

T2 <- T2[T2$z <= 0.934653, ]
T2F <- T2[T2$Gtype == "B",]
T2L <- T2[T2$Gtype == "S",]
T2XY <- readOGR(dsn="FL_outline", layer="FL_T2_outline")
T2XYmainowin <- as.owin(T2XY)
T2XZmainowin <- owin(xrange=c(808.0319, 808.0379), yrange=c(0.933872, 0.934653))
T2YZmainowin <- owin(yrange=c(629.7437, 629.7496), xrange=c(0.933872, 0.934653))

T3 <- T3[T3$z <= 0.9355, ]
T3F <- T3[T3$Gtype == "B",]
T3L <- T3[T3$Gtype == "S",]
T3XY <- readOGR(dsn="FL_outline", layer="FL_T3_outline")
T3XYmainowin <- as.owin(T3XY)
T3XZmainowin <- owin(xrange=c(808.0938, 808.0990), yrange=c(0.934352, 0.935475))
T3YZmainowin <- owin(yrange=c(629.7165, 629.7219), xrange=c(0.934352, 0.935475))

TOK <- TOK[TOK$z >= 0.937, ]
TOKF <- TOK[TOK$Gtype == "B",]
TOKL <- TOK[TOK$Gtype == "S",]
TOKXY <- readOGR(dsn="FL_outline", layer="FL_TOK_outline")
TOKXYmainowin <- as.owin(TOKXY)
TOKXZmainowin <- owin(xrange=c(808.1184, 808.1262), yrange=c(0.937925, 0.938961))
TOKYZmainowin <- owin(yrange=c(629.7052, 629.7140), xrange=c(0.937925, 0.938961))

##Rounding
#The whole assemblage
ggplot(alldata, aes(x= Roundness,  group=Trench)) + 
  geom_bar(aes(y = after_stat(prop), fill = factor(after_stat(x))), stat="count") +
  geom_text(aes(label = paste0("n=", after_stat(count)),
                y= ..prop..), stat= "count", vjust = 1.5, size=3.9) +
  geom_text(aes(label = scales::percent(..prop..),
                y= ..prop.. ), stat= "count", vjust = -.5, size=3.9) +
  labs(y = "Percent", x="Degree of edge abrasion") +
  facet_grid(~Trench) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks =c(0, 1, 2, 3), labels = c("Fresh", "Slightly\nabraded", "Abraded", "Heavily\nabraded")) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(margin = margin(t = 10)))

#Length group
data <- alldata[complete.cases(alldata$Roundness), ]
ggplot(data, aes(x=Lgroup, group=Roundness)) +
  geom_bar(aes(fill=factor(after_stat(group))), position = position_fill(reverse = TRUE)) +
  facet_grid(~Trench) +
  labs(fill="Degree of\nedge abrasion") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = paste0("n=", after_stat(count))), stat= "count", 
            vjust = 0.75, size=3.9, position=position_fill(0.9, reverse = TRUE)) +
  #geom_text(aes(label = scales::percent(..prop..)), stat="count", vjust = 0.75, size=3.9, position=position_fill(0.9, reverse = TRUE)) +
  labs(y = "Percent", x="Length group") +
  scale_x_discrete(limits=c("s", "m", "l"), labels=c("<20mm", "20-50mm", ">50mm")) +
  scale_fill_discrete(labels=c("Fresh", "Slightly abraded", "Abraded", "Heavily abraded")) +
  theme(axis.title.x = element_text(margin = margin(t = 10))) 

#Weight group
ggplot(data, aes(x=Wgroup, group=Roundness)) +
  geom_bar(aes(fill=factor(after_stat(group))), position = position_fill(reverse = TRUE)) +
  facet_grid(~Trench) +
  labs(fill="Degree of\nedge abrasion") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = paste0("n=", after_stat(count))), stat= "count", 
            vjust = 0.75, size=3.9, position=position_fill(0.9, reverse = TRUE)) +
  labs(y = "Percent", x="Weight group") +
  scale_x_discrete(limits=c("s", "m", "l"), labels=c("<5g", "5-50g", ">50g")) +
  scale_fill_discrete(labels=c("Fresh", "Slightly abraded", "Abraded", "Heavily abraded")) +
  theme(axis.title.x = element_text(margin = margin(t = 10))) 

#Technological group
ggplot(data, aes(x=Technological.category, group=Roundness)) +
  geom_bar(aes(fill=factor(after_stat(group))), position = position_fill(reverse = TRUE)) +
  facet_grid(~Trench) +
  labs(fill="Degree of\nedge abrasion") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = paste0("n=", after_stat(count))), stat= "count", 
            vjust = 0.5, size=3.9, position=position_fill(0.9, reverse = TRUE)) +
  labs(y = "Percent", x="Technological group") +
  scale_fill_discrete(labels=c("Fresh", "Slightly abraded", "Abraded", "Heavily abraded")) +
  theme(axis.title.x = element_text(margin = margin(t = 10))) +
  scale_x_discrete(limits=c("BP", "BPF", "Chunk", "Core", "Tool"), labels=c("Flake", "Flake\nfragment", "Chunk", "Core", "Retouched\ntool")) 

##Refitting
#Read data
T1_refits_weight <- read.csv("FL_data/FL_Refit_T1_weight.csv", header=TRUE, na="NA", fileEncoding = "UTF-8")
T1XY <- readOGR(dsn="FL_outline", layer="FL_T1_outline")
TOK_refits_weight <- read.csv("FL_data/FL_Refit_TOK_weight.csv", header=TRUE, na="NA", fileEncoding = "UTF-8")
TOKXY <- readOGR(dsn="FL_outline", layer="FL_TOK_outline")

#Create function and draw refit map
#Create function
drawrefits3viewweight.fl <- function(dataframe, refitsdataframe, XY){
  par(mar=c(2,1,1,1))
  dataframef <- dataframe[dataframe$Type == "B", ]
  dataframel <- dataframe[dataframe$Type == "S", ]
  #NS
  coordinates(dataframe) <- ~x+z
  coordinates(dataframef) <- ~x+z
  coordinates(dataframel) <- ~x+z
  par(fig=c(0.3, 0.99, 0, 0.3))
  plot(dataframe, pch=1, cex=0.3, col="white", axes=TRUE)
  points(dataframef, pch=19, cex=0.3, col=rgb(253, 231, 76, maxColorValue=255))
  points(dataframel, pch=19, cex=0.3, col=rgb(155, 197, 61, maxColorValue=255))
  #title("N-S", adj=0.5, line=-1)
  for (x in (1:nrow(refitsdataframe))){
    arrows(refitsdataframe[x,]$x1, refitsdataframe[x,]$z1, refitsdataframe[x,]$x2, refitsdataframe[x,]$z2, length=0.12, col=rgb(229, 89, 52, maxColorValue=255))
  }
  dataframe <- as.data.frame(dataframe)
  dataframef <- as.data.frame(dataframef)
  dataframel <- as.data.frame(dataframel)
  #WE
  coordinates(dataframe) <- ~z+y
  coordinates(dataframef) <- ~z+y
  coordinates(dataframel) <- ~z+y
  par(fig=c(0.1, 0.3, 0.3, 1), new=TRUE)
  plot(dataframe, pch=1, cex=0.3, col="white", axes=TRUE)
  points(dataframef, pch=19, cex=0.3, col=rgb(253, 231, 76, maxColorValue=255))
  points(dataframel, pch=19, cex=0.3, col=rgb(155, 197, 61, maxColorValue=255))
  #title("W-E", adj=0.5, line=0)
  for (x in (1:nrow(refitsdataframe))){
    arrows(refitsdataframe[x,]$z1, refitsdataframe[x,]$y1, refitsdataframe[x,]$z2, refitsdataframe[x,]$y2, length=0.12, col=rgb(229, 89, 52, maxColorValue=255))
  }
  dataframe <- as.data.frame(dataframe)
  dataframef <- as.data.frame(dataframef)
  dataframel <- as.data.frame(dataframel)
  #plan view
  coordinates(dataframe) <- ~x+y
  coordinates(dataframef) <- ~x+y
  coordinates(dataframel) <- ~x+y
  par(fig=c(0.3, 0.99, 0.3, 1), new=TRUE)
  plot(XY, border="grey25", col=FALSE, axes=TRUE)
  points(dataframef, pch=19, cex=0.5, col=rgb(253, 231, 76, maxColorValue=255))
  points(dataframel, pch=19, cex=0.5, col=rgb(155, 197, 61, maxColorValue=255))
  #title("T1 main plan", adj=0.5, line=0)
  for (x in (1:nrow(refitsdataframe))){
    arrows(refitsdataframe[x,]$x1, refitsdataframe[x,]$y1, refitsdataframe[x,]$x2, refitsdataframe[x,]$y2, length=0.12,col=rgb(229, 89, 52, maxColorValue=255))
  }
  #legend
  par(fig=c(0, 0.3, 0, 0.3), new=TRUE)
  plot(dataframef, col="white")
  dataframe <- as.data.frame(dataframe)
  dataframef <- as.data.frame(dataframef)
  dataframel <- as.data.frame(dataframel)
  legend("center", legend=c("Lithics", "Fossils", "refit lines"), col=c(rgb(155, 197, 61, maxColorValue=255), rgb(253, 231, 76, maxColorValue=255), rgb(229, 89, 52, maxColorValue=255)), bty="n", pch=c(19, 19, NA), cex=0.5)
  par(font=5)
  legend("left", legend=c(NA, NA, NA), col=c(rgb(155, 197, 61, maxColorValue=255), rgb(253, 231, 76, maxColorValue=255), rgb(229, 89, 52, maxColorValue=255)), bty="n", pch=c(NA, NA, 174), cex=0.5)
  par(font=1)
  par(mfrow=c(1,1))
}

drawrefits3viewweight.fl(T1, T1_refits_weight, T1XY)
title("T1")
drawrefits3viewweight.fl(TOK, TOK_refits_weight, TOKXY)
title("TOK")

#Draw rose diagram
T1wrefitsCO <- T1_refits_weight[complete.cases(T1_refits_weight$azimuth),]
rose(T1wrefitsCO$azimuth, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="")

TOKwrefitsCO <- TOK_refits_weight[complete.cases(TOK_refits_weight$azimuth),]
rose(TOKwrefitsCO$azimuth, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="")

#Statistic test
T1wrefitsCORC <- circular(T1wrefitsCO$azimuth, type="angles", units="degree")
rayleigh.test(T1wrefitsCORC)
rayleigh.test(T1wrefitsCO$azimuth)
kuiper.test(T1wrefitsCORC)

TOKwrefitsCORC <- circular(TOKwrefitsCO$azimuth, type="angles", units="degree")
rayleigh.test(TOKwrefitsCORC)
rayleigh.test(TOKwrefitsCO$azimuth)
kuiper.test(TOKwrefitsCORC)

###Spatial analysis
##Kernel density
#kdens for all and fossils
#create function
kdensitymap.xy.75contour4fossils <- function(dataframe, xyowin){
  par(mar=c(1,1,1,1))
  dataframef <- dataframe[dataframe$Type == "B", ]
  #plan view
  coordinates(dataframef) <- ~x+y
  dfxyppp <- as.ppp(coordinates(dataframef), as.owin(xyowin))
  dfxydensity <- density(dfxyppp, sigma=0.0005, edge=TRUE)
  dfxykdens <- dfxydensity/max(dfxydensity)
  dfconxy <- rasterToContour(raster(dfxykdens), nlevels=20)
  dfconxy <- dfconxy[dfconxy$level == 0.75, ]
  
  coordinates(dataframe) <- ~x+y
  xyppp <- as.ppp(coordinates(dataframe), as.owin(xyowin))
  xydensity <- density(xyppp, edge=TRUE, sigma=0.0005)
  plot(xydensity, legend=TRUE, useRaster=FALSE, main="", xlab=FALSE, ylab=FALSE)
  plot(dfconxy, col="grey25", add=TRUE, lwd=3.9)
  #title(dataframe, adj=0.5, line=0)
  dataframe <- as.data.frame(dataframe)
  dataframef <- as.data.frame(dataframef)
}

kdensitymap.xy.75contour4fossils(T1, T1XYmainowin)
#scale bar
polygon(c(808.059,808.060,808.060,808.059,808.059),c(629.7353,629.7353,629.73535,629.73535,629.7353))#0.001#0.00005
text(808.0595, 629.7351, "1 m", cex=0.6, col="black")#-0.0002
#north arrow
lines(c(808.0595,808.0595),c(629.7355,629.736))#0.0002#0.0005
polygon(c(808.0595,808.05945,808.0595,808.05955,808.0595), c(629.736,629.7358,629.736,629.7358,629.736), col="black", border="black")
#0.00005#0.0002
text(808.0595,629.7362, "N", cex=0.6, col="black")#0.0002

kdensitymap.xy.75contour4fossils(T2, T2XYmainowin)
#scale bar
polygon(c(808.032,808.033,808.033,808.032,808.032),c(629.7443,629.7443,629.74435,629.74435,629.7443))#0.001#0.00005
text(808.0325, 629.7441, "1 m", cex=0.6, col="black")#-0.0002
#north arrow
lines(c(808.0325,808.0325),c(629.7445,629.745))#0.0002#0.0005
polygon(c(808.0325,808.03245,808.0325,808.03255,808.0325), c(629.745,629.7448,629.745,629.7448,629.745), col="black", border="black")
#+—0.00005#0.0002
text(808.0325,629.7452, "N", cex=0.6, col="black")#0.0002

kdensitymap.xy.75contour4fossils(T3, T3XYmainowin)
#scale bar
polygon(c(808.094,808.095,808.095,808.094,808.094),c(629.7165,629.7165,629.71655,629.71655,629.7165))#0.001#0.00005
text(808.0945, 629.7163, "1 m", cex=0.6, col="black")#-0.0002
#north arrow
lines(c(808.0945,808.0945),c(629.7167,629.7172))#0.0002#0.0005
polygon(c(808.0945,808.09445,808.0945,808.09455,808.0945), c(629.7172,629.717,629.7172,629.717,629.7172), col="black", border="black")
#+—0.00005#0.0002
text(808.0945,629.7174, "N", cex=0.6, col="black")#0.0002

kdensitymap.xy.75contour4fossils(TOK, TOKXYmainowin)
#scale bar
polygon(c(808.1188,808.1198,808.1198,808.1188,808.1188),c(629.7057,629.7057,629.70575,629.70575,629.7057))#0.001#0.00005
text(808.1193, 629.70545, "1 m", cex=0.55, col="black")#-0.0002
#north arrow
lines(c(808.1193,808.1193),c(629.7059,629.7064))#0.0002#0.0005
polygon(c(808.1193,808.11925,808.1193,808.11935,808.1193), c(629.7064,629.7062,629.7064,629.7062,629.7064), col="black", border="black")
#+—0.00005#0.0002
text(808.1193,629.7067, "N", cex=0.55, col="black")#0.0002

#kdens for lithics
kdensitymap.xy <- function(dataframe, xyowin){
  par(mar=c(1,1,1,1))
  #plan view
  coordinates(dataframe) <- ~x+y
  xyppp <- as.ppp(coordinates(dataframe), as.owin(xyowin))
  xydensity <- density(xyppp, edge=TRUE, sigma=0.0005)
  plot(xydensity, legend=TRUE, useRaster=FALSE, main="", xlab=FALSE, ylab=FALSE)
  #title(dataframe, adj=0.5, line=0)
  dataframe <- as.data.frame(dataframe)
}#create function

kdensitymap.xy(T1L, T1XYmainowin)
#scale bar
polygon(c(808.059,808.060,808.060,808.059,808.059),c(629.7353,629.7353,629.73535,629.73535,629.7353))#0.001#0.00005
text(808.0595, 629.7351, "1 m", cex=0.6, col="black")#-0.0002
#north arrow
lines(c(808.0595,808.0595),c(629.7355,629.736))#0.0002#0.0005
polygon(c(808.0595,808.05945,808.0595,808.05955,808.0595), c(629.736,629.7358,629.736,629.7358,629.736), col="black", border="black")
#0.00005#0.0002
text(808.0595,629.7362, "N", cex=0.6, col="black")#0.0002
#title("T1 kernel density surface of lithics")

kdensitymap.xy(T2L, T2XYmainowin)
#scale bar
polygon(c(808.032,808.033,808.033,808.032,808.032),c(629.7443,629.7443,629.74435,629.74435,629.7443))#0.001#0.00005
text(808.0325, 629.7441, "1 m", cex=0.6, col="black")#-0.0002
#north arrow
lines(c(808.0325,808.0325),c(629.7445,629.745))#0.0002#0.0005
polygon(c(808.0325,808.03245,808.0325,808.03255,808.0325), c(629.745,629.7448,629.745,629.7448,629.745), col="black", border="black")
#+—0.00005#0.0002
text(808.0325,629.7452, "N", cex=0.6, col="black")#0.0002
#title("T2 kernel density surface of lithics")

kdensitymap.xy(T3L, T3XYmainowin)
#scale bar
polygon(c(808.094,808.095,808.095,808.094,808.094),c(629.7165,629.7165,629.71655,629.71655,629.7165))#0.001#0.00005
text(808.0945, 629.7163, "1 m", cex=0.6, col="black")#-0.0002
#north arrow
lines(c(808.0945,808.0945),c(629.7167,629.7172))#0.0002#0.0005
polygon(c(808.0945,808.09445,808.0945,808.09455,808.0945), c(629.7172,629.717,629.7172,629.717,629.7172), col="black", border="black")
#+—0.00005#0.0002
text(808.0945,629.7174, "N", cex=0.6, col="black")#0.0002
#title("T3 kernel density surface of lithics")

kdensitymap.xy(TOKL, TOKXYmainowin)
#scale bar
polygon(c(808.1188,808.1198,808.1198,808.1188,808.1188),c(629.7057,629.7057,629.70575,629.70575,629.7057))#0.001#0.00005
text(808.1193, 629.70545, "1 m", cex=0.55, col="black")#-0.0002
#north arrow
lines(c(808.1193,808.1193),c(629.7059,629.7064))#0.0002#0.0005
polygon(c(808.1193,808.11925,808.1193,808.11935,808.1193), c(629.7064,629.7062,629.7064,629.7062,629.7064), col="black", border="black")
#+—0.00005#0.0002
text(808.1193,629.7067, "N", cex=0.55, col="black")#0.0002
#title("TOK kernel density surface of lithics")

##Local G
#create function
localGplot <- function(dataframe, vector, XY){
  par(mar=c(1,1,1,1))
  #plan view_xy
  #generate a matrix of inverse distance weights 
  dataframe.xydis <- as.matrix(dist(cbind(dataframe$x, dataframe$y)))
  dataframe.xydis.inv <- 1/dataframe.xydis
  diag(dataframe.xydis.inv) <- 0
  dataframe.xydis.inv[is.infinite(dataframe.xydis.inv)] <- 0
  dataframe.xydis.inv.listw <- mat2listw(dataframe.xydis.inv)
  #localG
  g.dataframe <- localG(x=vector, listw=dataframe.xydis.inv.listw, zero.policy=FALSE)
  g.dataframe <- as.numeric(g.dataframe)
  g.dataframe <- as.data.frame(g.dataframe)
  g.dataframe <- as.data.frame(cbind(dataframe, g.dataframe))
  coordinates(g.dataframe) <- ~x+y
  breaks = c(-100, -2.58, -1.96, 1.96, 2.58, 100)
  palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
  col = palette[cut(g.dataframe$g.dataframe, breaks)]
  
  plot(XY, border="black", col=FALSE, axes=FALSE)
  plot(g.dataframe, pch=1, cex=1.5, col="black", add=TRUE)
  points(g.dataframe, pch=19, cex=1.4, col=col)
  #legend
  #par(fig=c(0.05, 0.25, 0, 0.3), new=TRUE)
  #plot(1,1, col="white", axes=FALSE)
  #legend("left", legend=c("Hot spot-99% confidence", "Hot spot-95% confidence", "Not significant", "Cold spot-95% confidence", "Cold spot-99% confidence"), col=c("#FF000080","#FF808080","grey25","#8080FF80","#0000FF80"), pch=c(19, 19, 1, 19, 19), cex=0.7, y.intersp=0.3, x.intersp=0.2, bty="n")
  par(mfrow=c(1,1))
}

#Rounding
T1wr <- T1[complete.cases(T1$Roundness), ]
T1wr <- T1wr[, c("x", "y", "z", "Roundness")]
localGplot(T1wr, T1wr$Roundness, T1XY)
#scale bar
polygon(c(808.059,808.060,808.060,808.059,808.059),c(629.7353,629.7353,629.73535,629.73535,629.7353))
text(808.0595, 629.7351, "1 m", cex=0.6, col="black")
#north arrow
lines(c(808.0595,808.0595),c(629.7355,629.736))
polygon(c(808.0595,808.05945,808.0595,808.05955,808.0595), c(629.736,629.7358,629.736,629.7358,629.736), col="black", border="black")
text(808.0595,629.7362, "N", cex=0.6, col="black")
#title(main="FL T1 Roundness localG")#FL_T1_localG_roundness

T2wr <- T2[complete.cases(T2$Roundness), ]
T2wr <- T2wr[, c("x", "y", "z", "Roundness")]
localGplot(T2wr, T2wr$Roundness, T2XY)
#scale bar
polygon(c(808.032,808.033,808.033,808.032,808.032),c(629.7443,629.7443,629.74435,629.74435,629.7443))#0.001#0.00005
text(808.0325, 629.7441, "1 m", cex=0.6, col="black")#-0.0002
#north arrow
lines(c(808.0325,808.0325),c(629.7445,629.745))#0.0002#0.0005
polygon(c(808.0325,808.03245,808.0325,808.03255,808.0325), c(629.745,629.7448,629.745,629.7448,629.745), col="black", border="black")
#+—0.00005#0.0002
text(808.0325,629.7452, "N", cex=0.6, col="black")#0.0002
#title(main="FL T2 Roundness localG")

T3wr <- T3[complete.cases(T3$Roundness), ]
T3wr <- T3wr[, c("x", "y", "z", "Roundness")]
localGplot(T3wr, T3wr$Roundness, T3XY)
#scale bar
polygon(c(808.094,808.095,808.095,808.094,808.094),c(629.7165,629.7165,629.71655,629.71655,629.7165))#0.001#0.00005
text(808.0945, 629.7163, "1 m", cex=0.6, col="black")#-0.0002
#north arrow
lines(c(808.0945,808.0945),c(629.7167,629.7172))#0.0002#0.0005
polygon(c(808.0945,808.09445,808.0945,808.09455,808.0945), c(629.7172,629.717,629.7172,629.717,629.7172), col="black", border="black")
#+—0.00005#0.0002
text(808.0945,629.7174, "N", cex=0.6, col="black")#0.0002
#title(main="FL T3 Roundness localG")

TOKwr <- TOK[complete.cases(TOK$Roundness), ]
TOKwr <- TOKwr[, c("x", "y", "z", "Roundness")]
localGplot(TOKwr, TOKwr$Roundness, TOKXY)
#scale bar
polygon(c(808.1188,808.1198,808.1198,808.1188,808.1188),c(629.7057,629.7057,629.70575,629.70575,629.7057))#0.001#0.00005
text(808.1193, 629.70545, "1 m", cex=0.55, col="black")#-0.0002
#north arrow
lines(c(808.1193,808.1193),c(629.7059,629.7064))#0.0002#0.0005
polygon(c(808.1193,808.11925,808.1193,808.11935,808.1193), c(629.7064,629.7062,629.7064,629.7062,629.7064), col="black", border="black")
#+—0.00005#0.0002
text(808.1193,629.7067, "N", cex=0.55, col="black")#0.0002
#title(main="FL TOK Roundness localG")#239X200

#max length
T1wl <- T1[complete.cases(T1$max), ]
T1wl <- T1wl[, c("x", "y", "z", "max")]
localGplot(T1wl, T1wl$max, T1XY)
#scale bar
polygon(c(808.059,808.060,808.060,808.059,808.059),c(629.7353,629.7353,629.73535,629.73535,629.7353))
text(808.0595, 629.7351, "1 m", cex=0.6, col="black")
#north arrow
lines(c(808.0595,808.0595),c(629.7355,629.736))
polygon(c(808.0595,808.05945,808.0595,808.05955,808.0595), c(629.736,629.7358,629.736,629.7358,629.736), col="black", border="black")
text(808.0595,629.7362, "N", cex=0.6, col="black")
#title(main="FL T1 Max Length localG")

T2wl <- T2[complete.cases(T2$max), ]
T2wl <- T2wl[, c("x", "y", "z", "max")]
localGplot(T2wl, T2wl$max, T2XY)
#scale bar
polygon(c(808.032,808.033,808.033,808.032,808.032),c(629.7443,629.7443,629.74435,629.74435,629.7443))#0.001#0.00005
text(808.0325, 629.7441, "1 m", cex=0.6, col="black")#-0.0002
#north arrow
lines(c(808.0325,808.0325),c(629.7445,629.745))#0.0002#0.0005
polygon(c(808.0325,808.03245,808.0325,808.03255,808.0325), c(629.745,629.7448,629.745,629.7448,629.745), col="black", border="black")
#+—0.00005#0.0002
text(808.0325,629.7452, "N", cex=0.6, col="black")#0.0002
#title(main="FL T2 Max Length localG")

T3wl <- T3[complete.cases(T3$max), ]
T3wl <- T3wl[, c("x", "y", "z", "max")]
localGplot(T3wl, T3wl$max, T3XY)
#scale bar
polygon(c(808.094,808.095,808.095,808.094,808.094),c(629.7165,629.7165,629.71655,629.71655,629.7165))#0.001#0.00005
text(808.0945, 629.7163, "1 m", cex=0.6, col="black")#-0.0002
#north arrow
lines(c(808.0945,808.0945),c(629.7167,629.7172))#0.0002#0.0005
polygon(c(808.0945,808.09445,808.0945,808.09455,808.0945), c(629.7172,629.717,629.7172,629.717,629.7172), col="black", border="black")
#+—0.00005#0.0002
text(808.0945,629.7174, "N", cex=0.6, col="black")#0.0002
#title(main="FL T3 Max Length localG")

TOKwl <- TOK[complete.cases(TOK$max), ]
TOKwl <- TOKwl[, c("x", "y", "z", "max")]
localGplot(TOKwl, TOKwl$max, TOKXY)
#scale bar
polygon(c(808.1188,808.1198,808.1198,808.1188,808.1188),c(629.7057,629.7057,629.70575,629.70575,629.7057))#0.001#0.00005
text(808.1193, 629.70545, "1 m", cex=0.55, col="black")#-0.0002
#north arrow
lines(c(808.1193,808.1193),c(629.7059,629.7064))#0.0002#0.0005
polygon(c(808.1193,808.11925,808.1193,808.11935,808.1193), c(629.7064,629.7062,629.7064,629.7062,629.7064), col="black", border="black")
#+—0.00005#0.0002
text(808.1193,629.7067, "N", cex=0.55, col="black")#0.0002
#title(main="FL TOK Max length localG")

#Weight
T1ww <- T1[complete.cases(T1$Weight), ]
T1ww <- T1ww[, c("x", "y", "z", "Weight")]
localGplot(T1ww, T1ww$Weight, T1XY)
#scale bar
polygon(c(808.059,808.060,808.060,808.059,808.059),c(629.7353,629.7353,629.73535,629.73535,629.7353))
text(808.0595, 629.7351, "1 m", cex=0.6, col="black")
#north arrow
lines(c(808.0595,808.0595),c(629.7355,629.736))
polygon(c(808.0595,808.05945,808.0595,808.05955,808.0595), c(629.736,629.7358,629.736,629.7358,629.736), col="black", border="black")
text(808.0595,629.7362, "N", cex=0.6, col="black")
#title(main="FL T1 Weight localG")

T2ww <- T2[complete.cases(T2$Weight), ]
T2ww <- T2ww[, c("x", "y", "z", "Weight")]
localGplot(T2ww, T2ww$Weight, T2XY)
#scale bar
polygon(c(808.032,808.033,808.033,808.032,808.032),c(629.7443,629.7443,629.74435,629.74435,629.7443))#0.001#0.00005
text(808.0325, 629.7441, "1 m", cex=0.6, col="black")#-0.0002
#north arrow
lines(c(808.0325,808.0325),c(629.7445,629.745))#0.0002#0.0005
polygon(c(808.0325,808.03245,808.0325,808.03255,808.0325), c(629.745,629.7448,629.745,629.7448,629.745), col="black", border="black")
#+—0.00005#0.0002
text(808.0325,629.7452, "N", cex=0.6, col="black")#0.0002
#title(main="FL T2 Weight localG")

T3ww <- T3[complete.cases(T3$Weight), ]
T3ww <- T3ww[, c("x", "y", "z", "Weight")]
localGplot(T3ww, T3ww$Weight, T3XY)
#scale bar
polygon(c(808.094,808.095,808.095,808.094,808.094),c(629.7165,629.7165,629.71655,629.71655,629.7165))#0.001#0.00005
text(808.0945, 629.7163, "1 m", cex=0.6, col="black")#-0.0002
#north arrow
lines(c(808.0945,808.0945),c(629.7167,629.7172))#0.0002#0.0005
polygon(c(808.0945,808.09445,808.0945,808.09455,808.0945), c(629.7172,629.717,629.7172,629.717,629.7172), col="black", border="black")
#+—0.00005#0.0002
text(808.0945,629.7174, "N", cex=0.6, col="black")#0.0002
#title(main="FL T3 Weight localG")

TOKww <- TOK[complete.cases(TOK$Weight), ]
TOKww <- TOKww[, c("x", "y", "z", "Weight")]
localGplot(TOKww, TOKww$Weight, TOKXY)
#scale bar
polygon(c(808.1188,808.1198,808.1198,808.1188,808.1188),c(629.7057,629.7057,629.70575,629.70575,629.7057))#0.001#0.00005
text(808.1193, 629.70545, "1 m", cex=0.55, col="black")#-0.0002
#north arrow
lines(c(808.1193,808.1193),c(629.7059,629.7064))#0.0002#0.0005
polygon(c(808.1193,808.11925,808.1193,808.11935,808.1193), c(629.7064,629.7062,629.7064,629.7062,629.7064), col="black", border="black")
#+—0.00005#0.0002
text(808.1193,629.7067, "N", cex=0.55, col="black")#0.0002
#title(main="FL TOK Weight localG")

##Chi-square test
#Create function
chisquaretest.2half.fl <- function(dataframe, XY){
  dataframef <- dataframe[dataframe$Gtype == "B", ]
  dataframel <- dataframe[dataframe$Gtype == "S", ]
  dataframe <- rbind(dataframef, dataframel)
  
  a <- (XY@bbox[1,1]+XY@bbox[1,2])/2
  b <- (XY@bbox[2,1]+XY@bbox[2,2])/2
  c <- (max(dataframe$z)-min(dataframe$z))/2+min(dataframe$z)
  
  upper <- dataframe[dataframe$z >= c, ]
  lower <- dataframe[dataframe$z < c, ]
  uppertable <- as.data.frame(table(upper$Gtype))
  lowertable <- as.data.frame(table(lower$Gtype))
  upplowtable <- data.frame(uppertable$Freq, lowertable$Freq, row.names=c("B", "S")) 
  upplow <- chisq.test(upplowtable)
  
  west <- dataframe[dataframe$x < a, ]
  east <- dataframe[dataframe$x >= a, ]
  westtable <- as.data.frame(table(west$Gtype))
  easttable <- as.data.frame(table(east$Gtype))
  westeasttable <- data.frame(westtable$Freq, easttable$Freq, row.names=c("B", "S")) 
  westeast <- chisq.test(westeasttable)
  
  north <- dataframe[dataframe$y >= b, ]
  south <- dataframe[dataframe$y < b, ]
  northtable <- as.data.frame(table(north$Gtype))
  southtable <- as.data.frame(table(south$Gtype))
  norsoutable <- data.frame(northtable$Freq, southtable$Freq, row.names=c("B", "S")) 
  norsou <- chisq.test(norsoutable)
  
  result <- c(upplow, westeast, norsou, use.names=TRUE)
  return(result)
}

chisquaretest.2half.fl(T1, T1XY)
chisquaretest.2half.fl(T2, T2XY)
chisquaretest.2half.fl(T3, T3XY)
chisquaretest.2half.fl(TOK, TOKXY)

##Correlation analysis
col2 = colorRampPalette(c('#67001F', '#B2182B', '#D6604D', '#F4A582',
                          '#FDDBC7', '#FFFFFF', '#D1E5F0', '#92C5DE',
                          '#4393C3', '#2166AC', '#053061'))

T1R <- T1[complete.cases(T1$Roundness),]
T1R <- T1R[, c(12,13,14,21,22,23,24,34)]
T1cor <- rcorr(as.matrix(T1R))
corrplot(T1cor$r, p.mat=T1cor$P, sig.level=0.05, insig="p-value", type="upper", col=rev(col2(200)), tl.pos="td", diag=FALSE, number.font = TRUE)
title("T1 (n=106)")

T2R <- T2[complete.cases(T2$Roundness),]
T2R <- T2R[complete.cases(T2R$Orientation),]
T2R <- T2R[complete.cases(T2R$Degree), ]
T2R <- T2R[, c(12,13,14,19,21,22,23,24,34)]
T2cor <- rcorr(as.matrix(T2R))
corrplot(T2cor$r, p.mat=T2cor$P, sig.level=0.05, insig="p-value", type="upper", col=rev(col2(200)), tl.pos="td", diag=FALSE, number.font = TRUE)
title("T2 (n=57)")

T3R <- T3[complete.cases(T3$Roundness),]
T3R <- T3R[, c(12,13,14,21,22,23,24,34)]
T3cor <- rcorr(as.matrix(T3R))
corrplot(T3cor$r, p.mat=T3cor$P, sig.level=0.05, insig="p-value", type="upper", col=rev(col2(200)), tl.pos="td", diag=FALSE, number.font = TRUE)
title("T3 (n=78)")

TOKR <- TOK[complete.cases(TOK$Roundness),]
TOKR <- TOKR[complete.cases(TOKR$Degree), ]
TOKR <- TOKR[, c(12,13,14,19,21,22,23,24,34)]
TOKcor <- rcorr(as.matrix(TOKR))
corrplot(TOKcor$r, p.mat=TOKcor$P, sig.level=0.05, insig="p-value", type="upper", col=rev(col2(200)), tl.pos="td", diag=FALSE, number.font = TRUE)
title("TOK (n=421)")

###Fabric analysis
##2 dimension
#Prepare data
T1CO <- T1[complete.cases(T1$Orientation),]#得全是数字
T1CO$Orientation2 <- (180+T1CO$Orientation)
T1CO$radian <- (T1CO$Orientation*pi)/180
T1CO$cradian <- T1CO$radian*2

T1COF <- T1CO[T1CO$Type == "B", ]
T1COL <- T1CO[T1CO$Type == "S", ]

T1COLM <- T1CO[which(T1CO$max > 20), ]
T1COLM <- T1COLM[which(T1COLM$lwratio >= 1.6), ]

T2CO <- T2[complete.cases(T2$Orientation),]#得全是数字
T2CO$Orientation2 <- (180+T2CO$Orientation)
T2CO$radian <- (T2CO$Orientation*pi)/180
T2CO$cradian <- T2CO$radian*2

T2COF <- T2CO[T2CO$Type == "B", ]
T2COL <- T2CO[T2CO$Type == "S", ]

T2COLM <- T2CO[which(T2CO$max > 20), ]
T2COLM <- T2COLM[which(T2COLM$lwratio >= 1.6), ]

T3CO <- T3[complete.cases(T3$Orientation),]#得全是数字
T3CO$Orientation2 <- (180+T3CO$Orientation)
T3CO$radian <- (T3CO$Orientation*pi)/180
T3CO$cradian <- T3CO$radian*2

T3COF <- T3CO[T3CO$Type == "B", ]
T3COL <- T3CO[T3CO$Type == "S", ]

T3COLM <- T3CO[which(T3CO$max > 20), ]
T3COLM <- T3COLM[which(T3COLM$lwratio >= 1.6), ]

TOKCO <- TOK[complete.cases(TOK$Orientation),]#得全是数字
TOKCO$Orientation2 <- (180+TOKCO$Orientation)
TOKCO$radian <- (TOKCO$Orientation*pi)/180
TOKCO$cradian <- TOKCO$radian*2

TOKCOF <- TOKCO[TOKCO$Type == "B", ]
TOKCOL <- TOKCO[TOKCO$Type == "S", ]

TOKCOLM <- TOKCO[which(TOKCO$max > 20), ]
TOKCOLM <- TOKCOLM[which(TOKCOLM$lwratio >= 1.6), ]

#Plot rose diagram
par(mfrow=c(4,4))
rose(T1CO$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL T1 all orientation (n=138)")
rose(T1CO$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)

rose(T1COF$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL T1 fossil orientation (n=116)")
rose(T1COF$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)

rose(T1COL$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL T1 lithic orientation (n=21)")
rose(T1COL$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)

rose(T1COLM$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL T1 qualified lithic orientation (n=6)")
rose(T1COLM$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)

rose(T2CO$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL T2 all orientation (n=544)")
rose(T2CO$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)

rose(T2COF$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL T2 fossil orientation (n=480)")
rose(T2COF$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)

rose(T2COL$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL T2 lithic orientation (n=63)")
rose(T2COL$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)

rose(T2COLM$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL T2 qualified lithic orientation (n=13)")
rose(T2COLM$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)

rose(T3CO$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL T3 all orientation (n=438)")
rose(T3CO$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)

rose(T3COF$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL T3 fossil orientation (n=379)")
rose(T3COF$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)

rose(T3COL$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL T3 lithic orientation (n=55)")
rose(T3COL$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)

rose(T3COLM$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL T3 qualified lithic orientation (n=15)")
rose(T3COLM$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)

rose(TOKCO$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL TOK all orientation (n=1239)")
rose(TOKCO$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)

rose(TOKCOF$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL TOK fossil orientation (n=738)")
rose(TOKCOF$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)

rose(TOKCOL$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL TOK lithic orientation (n=467)")
rose(TOKCOL$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)

rose(TOKCOLM$Orientation, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), main="FL TOK qualified lithic orientation (n=103)")
rose(TOKCOLM$Orientation2, col="grey", unit="degree", start="N", clockwise=TRUE, breaks=seq(0, 360, 10), add=TRUE)


#Statistic test
#Create function for Curray' L
CurrayL <- function(dataframe){
  dataframe <- dataframe[complete.cases(dataframe$Orientation), ]
  dataframe$sin2a <- sin((dataframe$Orientation*2))
  dataframe$cos2a <- cos((dataframe$Orientation*2))
  N <- length(dataframe$Orientation)
  R <- ((sum(dataframe$sin2a))^2+(sum(dataframe$cos2a))^2)^0.5
  L1 <- 100*R*(N^-1)
  P <- exp(1)^((-L1^2)*N*(10^-4))
  CurrayL <- data.frame(c(N, L1, P), row.names = c("n", "L", "p"))
  return(CurrayL)
}

#test
T1CORC <- circular(T1CO$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(T1CORC)
kuiper.test(T1CORC)
CurrayL(T1CO)

T1COFRC <- circular(T1COF$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(T1COFRC)
kuiper.test(T1COFRC)
CurrayL(T1COF)

T1COLRC <- circular(T1COL$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(T1COLRC)
kuiper.test(T1COLRC)
CurrayL(T1COL)

T1COLMRC <- circular(T1COLM$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(T1COLMRC)
kuiper.test(T1COLMRC)
CurrayL(T1COLM)

T2CORC <- circular(T2CO$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(T2CORC)
kuiper.test(T2CORC)
CurrayL(T2CO)

T2COFRC <- circular(T2COF$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(T2COFRC)
kuiper.test(T2COFRC)
CurrayL(T2COF)

T2COLRC <- circular(T2COL$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(T2COLRC)
kuiper.test(T2COLRC)
CurrayL(T2COL)

T2COLMRC <- circular(T2COLM$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(T2COLMRC)
kuiper.test(T2COLMRC)
CurrayL(T2COLM)

T3CORC <- circular(T3CO$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(T3CORC)
kuiper.test(T3CORC)
CurrayL(T3CO)

T3COFRC <- circular(T3COF$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(T3COFRC)
kuiper.test(T3COFRC)
CurrayL(T3COF)

T3COLRC <- circular(T3COL$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(T3COLRC)
kuiper.test(T3COLRC)
CurrayL(T3COL)

T3COLMRC <- circular(T3COLM$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(T3COLMRC)
kuiper.test(T3COLMRC)
CurrayL(T3COLM)

TOKCORC <- circular(TOKCO$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(TOKCORC)
kuiper.test(TOKCORC)
CurrayL(TOKCO)

TOKCOFRC <- circular(TOKCOF$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(TOKCOFRC)
kuiper.test(TOKCOFRC)
CurrayL(TOKCOF)

TOKCOLRC <- circular(TOKCOL$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(TOKCOLRC)
kuiper.test(TOKCOLRC)
CurrayL(TOKCOL)

TOKCOLMRC <- circular(TOKCOLM$cradian, type="angles", units="radians", modulo="2pi")
rayleigh.test(TOKCOLMRC)
kuiper.test(TOKCOLMRC)
CurrayL(TOKCOLM)

##3 dimensional
#Read data
setwd("/Users/XinDing/r")
fabric <- read.csv(file="FL_data/FL_fabric index.csv", header=TRUE, na="NA", fileEncoding="UTF-8")

#Plot figure
ggplot(fabric, aes(x=Type)) +
  geom_bar(aes(y=C, fill = factor(after_stat(x))), stat="identity", width=0.5) +
  facet_grid(~Trench) + #按trench分组
  labs(y = "Fabric strength", x="") +
  theme(legend.position="none") +
  geom_text(aes(label = paste0("n=", n),
                y=0), stat= "identity", vjust = 1.2, size=3.9) +
  geom_text(aes(label = round(C, digits=3),
                y= C ), stat= "identity", vjust = -.5, size=3.9) +
  scale_x_discrete(limits=c("All materials", "Fossils", "Lithics", "Qualified lithics"), labels=c("All\nmaterials", "Animal\nfossils", "Lithic\nartifacts", "Qualified\nlithics"))

ggplot(fabric, aes(x=Type)) +
  geom_bar(aes(y=F, fill = factor(after_stat(x))), stat="identity", width=0.5) +
  facet_grid(~Trench) + #按trench分组
  labs(y = "Flatness index", x="") +
  theme(legend.position="none") +
  geom_text(aes(label = paste0("n=", n),
                y=0), stat= "identity", vjust = 1.2, size=3.9) +
  geom_text(aes(label = round(F, digits=3),
                y= F ), stat= "identity", vjust = -.5, size=3.9) +
  scale_x_discrete(limits=c("All materials", "Fossils", "Lithics", "Qualified lithics"), labels=c("All\nmaterials", "Animal\nfossils", "Lithic\nartifacts", "Qualified\nlithics"))

ggplot(fabric, aes(x=Type)) +
  geom_bar(aes(y=CGI, fill = factor(after_stat(x))), stat="identity", width=0.5) +
  facet_grid(~Trench) + #按trench分组
  labs(y = "Cluster-girdle index", x="") +
  theme(legend.position="none") +
  geom_text(aes(label = paste0("n=", n),
                y=0), stat= "identity", vjust = 1.2, size=3.9) +
  geom_text(aes(label = round(CGI, digits=3),
                y=CGI), stat= "identity", vjust = -.5, size=3.9) +
  scale_x_discrete(limits=c("All materials", "Fossils", "Lithics", "Qualified lithics"), labels=c("All\nmaterials", "Animal\nfossils", "Lithic\nartifacts", "Qualified\nlithics"))

ggplot(fabric, aes(x=Type)) +
  geom_bar(aes(y=K, fill = factor(after_stat(x))), stat="identity", width=0.5) +
  facet_grid(~Trench) + #按trench分组
  labs(y = "K index", x="") +
  theme(legend.position="none") +
  geom_text(aes(label = paste0("n=", n),
                y=0), stat= "identity", vjust = 1.2, size=3.9) +
  geom_text(aes(label = round(K, digits=3),
                y= K ), stat= "identity", vjust = -.5, size=3.9) +
  scale_x_discrete(limits=c("All materials", "Fossils", "Lithics", "Qualified lithics"), labels=c("All\nmaterials", "Animal\nfossils", "Lithic\nartifacts", "Qualified\nlithics"))

