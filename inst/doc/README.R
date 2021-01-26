## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(MMD)
library(e1071)
library(bigmemory)
library(plyr)

## -----------------------------------------------------------------------------
datafile <- system.file("extdata", "Campylobacter_10SNP_HlW.csv", package = "MMD")
popfile <- system.file("extdata", "Campylobacter_10SNP_HlW.pop", package = "MMD")


## -----------------------------------------------------------------------------
datafile
popfile

## -----------------------------------------------------------------------------
ToAttribute <- "Human"

## -----------------------------------------------------------------------------
sourcenames <- c("Cattle","Chicken","Pig","Sheep","WB")

## -----------------------------------------------------------------------------
NL <- 100

## -----------------------------------------------------------------------------
attribution <- MMD_attr(datafile,popfile,NL,sourcenames,ToAttribute,verbose=FALSE)

## ---- eval=F------------------------------------------------------------------
#  attribution

## -----------------------------------------------------------------------------
barplot(height = attribution[[3]]$mean,names.arg = attribution[[3]]$Sources,las=2,cex.names = 0.8,main=NULL, ylab = expression("Attribution probability, p"[s]))

## ---- eval=FALSE--------------------------------------------------------------
#  datadir <- "~/Work/Source_Attribution/R_projects/project"
#  dataname <- "Campylobacter_10SNP_HlW"
#  NL <- 100
#  sourcenames <- c("Cattle","Chicken","Pig","Sheep","WB")

## -----------------------------------------------------------------------------
SelfAttribution <- MMD_attr(datafile,popfile,NL,sourcenames,ToAttribute = "Pig",SelfA = "yes",optq = "yes", pqmin = 0, pqmax = 0.5, np=1000, fSelfA = 0.5)

## ---- eval=FALSE--------------------------------------------------------------
#  SelfAttribution

## -----------------------------------------------------------------------------
barplot(height = SelfAttribution[[3]]$mean,names.arg = SelfAttribution[[3]]$Sources,las=2,cex.names = 0.8,main=NULL, ylab = expression("Attribution probability, p"[s]))

## -----------------------------------------------------------------------------
EntropyLoci <- MMD_Entropy(datafile,popfile,NL,sourcenames)

## ---- eval=FALSE--------------------------------------------------------------
#  EntropyLoci

## -----------------------------------------------------------------------------
histHlT <- hist(EntropyLoci[[6]]$HlT,col="#1C86EE",main=NULL,xlab=expression(paste("H"^"T")),ylab="Frequency")
histHlB <- hist(EntropyLoci[[6]]$HlB,col="#1C86EE",main=NULL,xlab=expression(paste("H"^"B")),ylab="Frequency")
histHlW <- hist(EntropyLoci[[6]]$HlW,col="#1C86EE",main=NULL,xlab=expression(paste("H"^"W")),ylab="Frequency")

## -----------------------------------------------------------------------------
plot(cbind(EntropyLoci[[6]]$HlW,EntropyLoci[[6]]$HlB),col="#1C86EE",main=NULL, xlab=expression(paste("H"^"W")),ylab=expression(paste("H"^"B")))

## -----------------------------------------------------------------------------
RedundancyLoci <- MMD_Rn(datafile,popfile,NL,sourcenames)

## ---- eval=F------------------------------------------------------------------
#  RedundancyLoci

## -----------------------------------------------------------------------------
plot(RedundancyLoci[[6]],col="#EE4000",main=NULL,log="x",type="l",xlab="n-th selected locus",ylab=expression("R"[n]))

## -----------------------------------------------------------------------------
plot(RedundancyLoci[[8]],col="#EE4000",main=NULL,log="x",type="l",xlab="n-th selected locus",ylab=expression("R"[n]))

