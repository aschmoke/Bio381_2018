library(plyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(reshape2)


d <- read.table(file="SampleData.csv",header=TRUE,sep=",")

head(d)
str(d)
table(d$SampleID)
length(unique(d$GeneSymbol))
length(unique(d$Reference))    
       
d$Ratio <- with(d,2^Log2Ratio)
head(d)

d$Row <- with(d, 1:nrow(d))

dSN <- subset(d,!Row %in% Row[SumSN < 20])
table(dSN$SampleID)

dMean <- ddply(dSN,.(Reference,SampleID,GeneSymbol),
               summarize,
               Count=length(Reference),
               MeanFC=mean(Ratio),
               StDev=sd(Ratio))
head(dMean)
length(unique(dSN$Reference))
length(unique(dMean$Reference))

dFilt <- subset(dMean, !Reference %in% Reference[Count<3])
dim(dFilt)
length(unique(dFilt$Reference))

dFiltWide <- dcast(dFilt,Reference ~ SampleID,value.var="Count")
dFiltWide
dFiltAB <- dFiltWide[which(complete.cases(dFiltWide)),]
dFiltAB
dFiltAB$mag <- dFiltAB$B/dFiltAB$A
dFiltAB
dFiltLargeB <- subset(dFiltAB,dFiltAB$mag>=5)
dFiltLargeB
dFiltAorB <- dFiltWide[which(!complete.cases(dFiltWide)),]
dFiltAorB
dFiltBs <- subset(dFiltAorB,dFiltAorB$B !="NA")
dFiltBs
namestosubset <- data.frame(Names=c(as.character(dFiltLargeB$Reference), as.character(dFiltBs$Reference)))
namestosubset
dFilt$Reference <- as.character(dFilt$Reference)
dFilt
dFiltBs <- subset(dFilt,dFilt$Sample=="B")
dFiltBs
dFiltNew <- subset(dFiltBs, dFiltBs$Reference %in% namestosubset$Names)
head(dFiltNew)
dim(dFiltNew)

norm <- dFiltNew$MeanFC[dFiltNew$GeneSymbol=="Prkaca"]
norm

dFiltNew$FCnorm <- with(dFiltNew,MeanFC/norm)
mNorm <- dFiltNew[dFiltNew$GeneSymbol=="Prkaca","FCnorm"]
mNorm                  
sNorm <- dFiltNew[dFiltNew$GeneSymbol=="Prkaca","StDev"]
sNorm
nNorm <- dFiltNew[dFiltNew$GeneSymbol=="Prkaca","Count"]
nNorm

dFiltNew$pVal <- with(dFiltNew,
                      2*pt(
                        -abs(
                          ((FCnorm - mNorm) /
                             sqrt((((Count-1)*StDev^2 + (nNorm-1)*sNorm^2)/
                                     (Count + nNorm - 2))*(1/Count+1/nNorm)))),
                        (Count + nNorm - 2)),
                                  lower=FALSE)
dFiltNew$FDR <- p.adjust(dFiltNew$pVal,method="BH")                           
dFiltNew$Log2FCnorm <- with(dFiltNew,log(FCnorm,2))
dFiltNew$LogFDR <- with(dFiltNew,-log(FDR,10))
dim(dFiltNew)

dFiltNew$threshold <- as.factor(abs(dFiltNew$Log2FCnorm)>2 &
                               dFiltNew$FDR < 0.05)

p1 <- ggplot(data=dFiltNew,aes(x=Log2FCnorm,y=LogFDR,color=threshold)) +
  geom_point() +
  scale_color_manual(values=c("darkgrey","red")) +
  theme_bw(base_size=12) +
  theme(legend.position="none") +
  geom_text_repel(data=subset(dFiltNew,FDR<0.05 &
                                abs(Log2FCnorm) >2),
                  aes(label=GeneSymbol),
                  size=2.5,
                  color="black",
                  box.padding=unit(0.35,"lines"),
                  point.padding=unit(0.3,"lines"))
p1                       

dFiltNew <- dFiltNew[order(dFiltNew$Log2FCnorm),]
dFiltNew$Rank <-  with(dFiltNew,1:nrow(dFiltNew))

p2 <- ggplot(data=dFiltNew,mapping=aes(x=Rank,y=Log2FCnorm,color=threshold)) +
  geom_point() +
  scale_color_manual(values=c("darkgrey","red")) +
  theme_bw(base_size=12) +
  theme(legend.position="none") +
  geom_text_repel(data=subset(dFiltNew,FDR<0.05 &
                                abs(Log2FCnorm) >2),
                  aes(label=GeneSymbol),
                  size=2.5,
                  color="black",
                  box.padding=unit(0.35,"lines"),
                  point.padding=unit(0.3,"lines"))
p2

p3 <- ggplot(data=dFiltNew,mapping=aes(Log2FCnorm,fill=I("dodgerblue"),color=I("black"))) +
  geom_histogram() +
  theme_bw(base_size=12)

p3 +{
  p2 +
    p1 +
    plot_layout(nrow=1)
} +
  plot_layout(ncol=1)
