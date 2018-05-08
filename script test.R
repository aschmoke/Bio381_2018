# Data filtering and analysis in R
# Anna Schmoker

library(plyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(reshape2)

d <- read.table(file="SampleData.csv",header=TRUE,sep=",")
head(d)

# Filters:
    # S/N
    # 3 or more peptides
    # remove peptides present in control (or 5X enriched in experimental)

# Normalize H/L to Pka-c
    # get mean for each protein
    # get SD

# p-value (increased/decreased binding)

# visual summary

table(d$SampleID)
length(unique(d$GeneSymbol))

d$Ratio <- with(d, 2^Log2Ratio)

d$Row <- with(d, 1:nrow(d))

dSN <- subset(d,!Row %in% Row[SumSN <20])

dMean <- ddply(dSN,
               .(Reference,SampleID,GeneSymbol),
               summarize,
               Count=length(Reference),
               MeanFC=mean(Ratio),
               StDev=sd(Ratio))
head(dMean)

dFilt <- subset(dMean, !Reference %in% Reference[Count<3])

# Lauren

dFiltW <- dcast(dFilt,Reference~SampleID,value.var="Count")
head(dFiltW)

dFiltAB <- dFiltW[which(complete.cases(dFiltW)),]

dFiltAB$div <- dFiltAB$B/dFiltAB$A

dFiltLB <- subset(dFiltAB,dFiltAB$div>=5)

dFiltAorB <- dFiltW[which(!complete.cases(dFiltW)),]

dFiltBs <- subset(dFiltAorB,dFiltAorB$B !="NA")

namestosubset <- data.frame(Names=c(as.character(dFiltLB$Reference), as.character(dFiltBs$Reference)))

dFiltBs <- subset(dFilt,dFilt$Sample=="B")
dFiltNew <- subset(dFiltBs,dFiltBs$Reference %in% namestosubset$Names)

head(dFiltNew)

norm <- dFiltNew$MeanFC[dFiltNew$GeneSymbol=="Prkaca"]
norm

dFiltNew$FCnorm <- with(dFiltNew, MeanFC/norm)

# p-value

mNorm <- dFiltNew[dFiltNew$GeneSymbol=="Prkaca","FCnorm"]
sNorm <- dFiltNew[dFiltNew$GeneSymbol=="Prkaca","StDev"]
nNorm <- dFiltNew[dFiltNew$GeneSymbol=="Prkaca","Count"]

dFiltNew$pVal <- with(dFiltNew,
                      2*pt(
                        -abs(
                          ((FCnorm - nNorm)/
                             sqrt((((Count-1)*StDev^2 +
                                      (nNorm-1)*sNorm^2)/
                                     (Count+nNorm-2))*(1/Count+1/nNorm)))),
                          (Count+nNorm-2)),lower=FALSE)

dFiltNew$FDR <- p.adjust(dFiltNew$pVal,method="BH")                      

dFiltNew$Log2FCnorm <- with(dFiltNew, log(FCnorm,2))
dFiltNew$LogFDR <- with(dFiltNew, -log(FDR,10))

dFiltNew$threshold <- as.factor(abs(dFiltNew$Log2FCnorm)>2 &
                                  dFiltNew$FDR < 0.05)

p1 <- ggplot(data=dFiltNew,mapping=aes(x=Log2FCnorm,y=LogFDR,color=threshold)) +
  geom_point() +
  scale_color_manual(values=c("darkgrey","red")) +
  theme_bw()
p1
