# SILAC test on Amanda's mouse data
# 28 March 2018
# AMS

# preliminaries
library(plyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(reshape2)

# ADD A SAMPLE ID LANE (A = mock, B = experimental)
# Find-and-replace "#" with "+"
# filter UniqueDeltaCorr
# filter false positives

# read table (in .csv format)
d <- read.table(file="SampleData.csv",header=TRUE,sep=",")
head(d)
tail(d)
str(d)
table(d$GeneSymbol) # Number rows with blank gene ID = 1st entry
table(d$SampleID)
length(unique(d$GeneSymbol))
length(unique(d$Reference)) # Number of unique proteins

d$Ratio <- with(d, 2^Log2Ratio)
head(d)
d$Ratio <- round(d$Ratio,4)
head(d)

d$Row <- with(d, 1:nrow(d))

dSN <- subset(d,!Row %in% Row[SumSN < 20]) # need to choose a column with unique value! or will remove duplicates (unless you want that)
head(dSN)
table(dSN$GeneSymbol)
length(unique(dSN$GeneSymbol))

# write.csv(dSN,file="dSN.csv")

dMean <- ddply(dSN,.(Reference,SampleID,GeneSymbol),
               summarize,
               Count=length(Reference),
               MeanFC=mean(Ratio),
               StDev=sd(Ratio)) # summarize data: mean Heavy-to-light ratio, other stats on that, and number of peptides identified per protien
dim(dMean)
head(dMean)
length(unique(dMean$GeneSymbol))

dFilt <- subset(dMean,!Reference %in% Reference[Count<3]) # subset excluding proteins identified by less than 3 peptides
str(dFilt)
table(dFilt$GeneSymbol)

norm <- dFilt$MeanFC[dFilt$GeneSymbol=="Prkaca"]
norm


dFilt$FCnorm <- with(dFilt, MeanFC/norm) # Normalize to loading control
head(dFilt)


mNorm <- dFilt[dFilt$GeneSymbol=="Prkaca","FCnorm"]
mNorm
sNorm <- dFilt[dFilt$GeneSymbol=="Prkaca","StDev"]
sNorm
nNorm <- dFilt[dFilt$GeneSymbol=="Prkaca","Count"]
nNorm

# Test t value

dFilt$pVal <- with(dFilt,
                   2*pt(
                     -abs(
                       ((FCnorm - mNorm)/
                          sqrt((((Count-1)*StDev^2 + (nNorm-1)*sNorm^2)/
                                  (Count+nNorm-2))*(1/Count+1/nNorm)))),
                     (Count + nNorm - 2)),lower=FALSE)

# Adjust p-value for multiple calculations

dFilt$FDR <- p.adjust(dFilt$pVal,method="BH")
head(dFilt)


dFilt$Log2FCnorm <- with(dFilt, log(FCnorm,2)) # Log2-transform the normalized fold change
head(dFilt)
tail(dFilt)

dFilt$LogFDR <- with(dFilt, -log(FDR,10))
head(dFilt)

#############################################################
# subset so only those with both A and B 
# divide B by A
# subset those that have > 5 in column
# rbind just distinct B's and B's greater than 5x their As
# Subset data with As removed based on large Bs


# Reshape into 'wide' format for easier manipulation
dFiltWide<-dcast(dFilt, Reference ~ SampleID, value.var = "Count")

# Extract those that have both A and B in Sample column
dFiltAB<-dFiltWide[which(complete.cases(dFiltWide)),]

# Create column that divides Sample B by A
dFiltAB$mag<-dFiltAB$B/dFiltAB$A

# Subset only those where B is at least 5 times greater than A
dFiltLargeB<-subset(dFiltAB, dFiltAB$mag>=5)

# Get names of those with either A or B
dFiltAorB<-dFiltWide[which(!complete.cases(dFiltWide)),]

# Subset only B samples
dFiltBs<-subset(dFiltAorB, dFiltAorB$B != "NA")

# Create a data frame that have References with either only B or where B > 5x A
namesToSubset<-data.frame(Names=c(as.character(dFiltLargeB$Reference), as.character(dFiltBs$Reference)))

dFilt$Reference<-as.character(dFilt$Reference)

# Subset original data to only include References that have B samples
dFiltBs<-subset(dFilt, dFilt$Sample=="B")
dFiltNew<-subset(dFiltBs, dFiltBs$Reference %in% namesToSubset$Names)

dFiltNew
############################################
# ggplot volcano plot
# Highlight genes that have an absolute fold change > 4 and a p-value < 0.05
dFiltNew$threshold <- as.factor(abs(dFiltNew$Log2FCnorm) > 1 &
                                  dFiltNew$FDR < 0.05)

# genes$Significant <- ifelse(genes$padj < 0.05, "FDR < 0.05", "Not Sig")

# Construct the plot object
g <- ggplot(data=dFiltNew, aes(x=Log2FCnorm, y=LogFDR, color=threshold)) +
  geom_point(size=1.75) +
  scale_color_manual(values = c("darkgrey", "red")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  xlim(c(-5, 5)) + ylim(c(0, 21)) +
  xlab("log2 fold change") + ylab("-log10 FDR")
g
g + geom_text(aes(label = dFiltNew$GeneSymbol),
              size = 3.5, color="black")
g + geom_text_repel(data=subset(dFiltNew,FDR < 0.05 &
                                  abs(Log2FCnorm) > 1),
                    aes(label = GeneSymbol),
                    size = 3.5,
                    color="black",
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines"))


dFiltNew <- dFiltNew[order(dFiltNew$Log2FCnorm),] # order rows by Log2-transformed mean fold change
head(dFiltNew)
tail(dFiltNew)

dFiltNew$Rank <- with(dFiltNew, 1:nrow(dFiltNew)) # add rank order


p1 <- ggplot(data=dFiltNew,mapping=aes(x=Rank,y=Log2FCnorm,color=threshold)) +
  geom_point(size=1.75) +
  scale_color_manual(values = c("darkgrey", "red")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  xlab("Rank") + ylab("log2 fold change") +
  geom_text_repel(data=subset(dFiltNew,FDR < 0.05 &
                                abs(Log2FCnorm) > 1),
                  aes(label = GeneSymbol),
                  size = 2.5,
                  color="black",
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"))

p2 <- ggplot(data=dFiltNew,mapping=aes(Log2FCnorm,fill=I("dodgerblue"), color=I("black"))) +
  geom_histogram() +
  theme_bw(base_size = 12)
p2

p3 <- g + geom_text_repel(data=subset(dFiltNew,FDR < 0.05 &
                                        abs(Log2FCnorm) > 1),
                          aes(label = GeneSymbol),
                          size = 2.5,
                          color="black",
                          box.padding = unit(0.35, "lines"),
                          point.padding = unit(0.3, "lines"))

p2 + {
  p1 +
    p3 +
    plot_layout(nrow=1)
} +
  plot_layout(ncol=1)


#################
# Upload onto Slack or Website
################
