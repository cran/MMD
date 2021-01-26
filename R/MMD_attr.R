#--- Documentation  (roxygen2 needed to compile)
#--
#' Attribution of individuals to sources using the MMD method
#'
#' @param datafile character; Name of the file *.csv (with full path in the file system) containing the genotypes (features) of individuals.
#' @param popfile character; Name of the file *.pop (with full path in the file system) containing the genotypes (features) of individuals.
#' @param NL integer; number of loci. If larger than the number of available loci in the data set, NL is reduced to the maximum available number of loci.
#' @param sourcenames a character vector listing the names of the sources.
#' @param ToAttribute character giving the name of the individuals of aknown origin (i.e. those that will be attributed to source).
#' @param SelfA character; if "no" attribution of individuals to sources is made; if "yes", self-attribution of selected individuals from sources is made. (Default "no")
#' @param fSelfA real number in the interval (0,1). When SelfA="yes", fSelfA specifies the fraction of individuals from the source specified by ToAttribute that will be assumed to be of unknown origin. (Default 0.1)
#' @param randomSelfA character only relevant if SelfA="yes". If "yes", individuals to be considered as unknown are randomly selected from the source specified by ToAttribute; if "no" a list of names for individuals is read from filepoplist. (Default "yes")
#' @param quantile real number with values in (0,1) giving the q-quantile for the MMD method. Only used if the quantile is not obatined through optimisation of the probability of correct self-attribution. (Default 0.01)
#' @param optq character; if "no", the specified quantile value is used; if "yes", the q-quantile is optimised (only meaningful for self-attribution so optq="no" automatically if SelfA="no"). (Default "no")
#' @param pqmin real number with values in (0,1); minimum value of q-quantile when optq="yes". (Default 0)
#' @param pqmax real number with values in (0,1); maximum value of q-quantile when optq="yes". (Default 0.5)
#' @param np integer giving the number of values of q-quantile in the interval (pqmin,pqmax) when optq="yes". (Default 20)
#' @param Nbootstrap integer giving the number of samples used for bootstrapping to estimate the uncertainty of the attribution probability $p_s$ bootstrap. (Default 10000)
#' @param verbose boolean (TRUE/FALSE) for the display of a progress bar  (Default FALSE)
#'
#' @return
#' If optq="yes", the output is a list with seven elements:
#' \enumerate{
#' \item Number of individuals from unknown origin.
#' \item Number of sources.
#' \item Statistics of the attribution probability to sources, $p_s$.
#' \item Probability of attribution of each unknown individual to each source $p_{u,s}$
#' \item Runtime of the calculation.
#' \item Number of loci.
#' \item Parameter q used to calculate the q-quantile of the Hamming distance in the MMD method.
#' \item Data frame giving the probability of correct attribution vs. q-quantile.
#' }
#'
#' If optq="no", the output list contains all the items in the list above except the last one.
#'
#' @examples
#' ## This example uses a small dataset stored in the MMD package
#' datafile <- system.file("extdata", "Campylobacter_10SNP_HlW.csv", package = "MMD")
#' popfile <- system.file("extdata", "Campylobacter_10SNP_HlW.pop", package = "MMD")
#'
#' NL <- 100
#' sourcenames <- c("Cattle","Chicken","Pig","Sheep","WB")
#'
#' ##----- Source attribution
#' ToAttribute <- "Human"
#' SelfA="no"
#' attribution <- MMD_attr(datafile,popfile,NL,sourcenames,ToAttribute)
#'
#' ## See more detailed examples in the vignette.
#'
#' @author Francisco J. Perez-Reche (Univeristy of Aberdeen)


## --- End of documentation ------------
## -------------------------------------

MMD_attr <- function(datafile,popfile,NL,sourcenames,ToAttribute,SelfA="no",fSelfA=0.5,randomSelfA="yes",quantile=0.01,optq="no",pqmin=0,pqmax=0.5,np=20,Nbootstrap=10000,verbose=FALSE){

  start_timeMSD <- Sys.time()

  # - Progress bar
  if(verbose) pb <- txtProgressBar(min = 0, max = (length(sourcenames)+2), initial = 0,style=3)

  rawdata <- read.big.matrix(datafile,sep=",",type="integer",header = F)
  NL <- min(NL,ncol(rawdata))
  rawdata <- rawdata[,c(1:NL)]
  popdata <- read.table(popfile)

  if(verbose) setTxtProgressBar(pb,1)

  indexSource0 <- which(popdata == ToAttribute)

  if(SelfA=="no"){
    optq<-"no"
    poplistSA <- popdata
    index_u <- which(poplistSA==ToAttribute)
    Nu <- length(index_u)
    popdatavec <- as.matrix(popdata)
    popdataNoAttr <- as.data.frame(popdatavec)
  }

  if(SelfA=="yes"){
    if(randomSelfA=="no"){
      poplistSA <- popdata
      index_u <- which(poplistSA==ToAttribute)
      Nu <- length(index_u)
    }

    if(randomSelfA=="yes"){
      Nu <- floor(fSelfA*length(indexSource0))  # number of individuals selected as unknown
      if(Nu==0){
        Nu <- 2
        rawdata <- as.matrix(rawdata)
        }
      index_u <- sample(indexSource0,Nu,replace = F)
      temp<-rep("unknown",Nu)
      popdatavec <- as.matrix(popdata)
      #poplistSA <- replace(popdatavec[,1],index_u,temp)
      popdatavec[index_u,1] <- "unknown"
      popdataNoAttr <- as.data.frame(popdatavec)
    }
  }

  # -- Random selection of loci - Not implemented here

  NumberLoci <- ncol(rawdata)
  Nsamples <- nrow(popdata)
  ListNames <- popdata
  #reservoirs <- unique(ListNames)$V1
  reservoirs <- as.matrix(unique(ListNames))
  Sources <- intersect(reservoirs,sourcenames)
  NR <- length(Sources)

  if(Nu>1) AttributionSamples <- as.big.matrix(rawdata[index_u,])
  if(Nu==1) AttributionSamples <- rawdata[index_u,]

  nunknown <- nrow(AttributionSamples)

  AttributeSourceIndex <- which.min(pmatch(Sources,ToAttribute))


  # --- Calculating the Hamming distance between individuals from unknown origin and sources
  hdlist <- vector("list", NR)
  for(k in 1:NR){
    #aux1 <-subset(rawdata,popdata==Sources[k])[,c(2:(NumberLoci+1))]
    if(Nu>1){
    aux1 <- as.big.matrix(rawdata[which(popdataNoAttr==Sources[k]),])
    hd1 <- big.matrix(nunknown,nrow(aux1))
    for(i in 1:nunknown){
      #print(i)
      for(j in 1:nrow(aux1)){
        hd1[i,j] <- hamming.distance(AttributionSamples[i,],aux1[j,])
      }
    }
    }
    if(Nu==1){
      aux1 <- rawdata[which(popdataNoAttr==Sources[k]),]
      hd1 <- vector(mode = "integer",nrow(aux1))
      for(j in 1:nrow(aux1)){
        #print(i
          hd1[j] <- hamming.distance(AttributionSamples,aux1[j,])
        }
    }


    hdlist[[k]]<-hd1
    pbcount<-k+1
    if(verbose) setTxtProgressBar(pb,pbcount)
  }

  ## --- Optimising q for source attribution
  if(optq == "yes"){
    dp = (pqmax - pqmin)/np
    pcorrect <- matrix(0,np+1,1)
    pquantile <- matrix(0,np+1,1)

    for(ip in 1:(np+1)){
      quantileTemp = pqmin+(ip-1)*dp
      pquantile[ip]=quantileTemp
      qlist <- matrix(-1,NR,nunknown)
      for(k in 1:NR){
        for(u in 1:nunknown){
          qlist[k,u] <- quantile(hdlist[[k]][u,],probs = quantileTemp)
        }
      }
      minq <- apply(qlist,2,min)
      sigmaus <- matrix(-1,NR,nunknown)
      for(k in 1:NR){
        for(u in 1:nunknown){
          cdfq <- ecdf(hdlist[[k]][u,])
          sigmaus[k,u] <- cdfq(minq[u])
        }
      }

      sumsigmaus <- apply(sigmaus,2,sum)
      pus <- matrix(-1,NR,nunknown)
      for(k in 1:NR){
        for(u in 1:nunknown){
          pus[k,u] <- sigmaus[k,u]/sumsigmaus[u]
        }
      }

      ps <- apply(pus,1,mean)
      pcorrect[ip] <-ps[AttributeSourceIndex]
    }
    argmaxposition <-which.max(pcorrect)
    quantile <- pqmin+(argmaxposition-1)*dp
    maxAttr <- pcorrect[argmaxposition]
  }

  # --- After optimising q (or when q is not optimised)
  qlist <- matrix(-1,NR,nunknown)
  for(k in 1:NR){
    for(u in 1:nunknown){
      qlist[k,u] <- quantile(hdlist[[k]][u,],probs = quantile,na.rm=T)
    }
  }

  minq <- apply(qlist,2,min)
  sigmaus <- matrix(-1,NR,nunknown)
  for(k in 1:NR){
    for(u in 1:nunknown){
      cdfq <- ecdf(hdlist[[k]][u,])
      sigmaus[k,u] <- cdfq(minq[u])
    }
  }

  sumsigmaus <- apply(sigmaus,2,sum)
  pus <- matrix(-1,NR,nunknown)

  for(k in 1:NR){
    for(u in 1:nunknown){
      pus[k,u] <- sigmaus[k,u]/sumsigmaus[u]
    }
  }

  ps <- apply(pus,1,mean)
  psSummary <- data.frame(Sources,ps)

  psStats <- matrix(-1,NR,4)
  for(k in 1:NR){
    aux3 <- replicate(Nbootstrap,mean(sample(pus[k,],nunknown,replace = T)))
    psStats[k,1]=mean(aux3)
    psStats[k,2]=sqrt(var(aux3))
    psStats[k,3]=quantile(aux3,probs=2.5/100)
    psStats[k,4]=quantile(aux3,probs=97.5/100)
  }
  colnames(psStats) <- c("mean","standard deviation","2.5 percentile","97.5 percentile")
  psStatsSummarytemp <- data.frame(Sources,psStats)
  psStatsSummary <- psStatsSummarytemp[with(psStatsSummarytemp, order(Sources)), ]

  rownames(pus)<-Sources

  aux1 <-  cbind(index_u,t(pus))
  #pus_list <- aux1[order(aux1[,1]),]
  pus_listDF <- data.frame(aux1[order(aux1[,1]),])
  colnames(pus_listDF) <- append("individual",Sources)

  end_timeMSD <- Sys.time()
  timeMessage <- end_timeMSD-start_timeMSD
  nunknouwnOut <- nunknown
  NROut <- NR

  if(optq=="yes"){
    pqvspcorrect <- cbind(pquantile,pcorrect)
    colnames(pqvspcorrect) <- c("q-quantile","prob. correct attribution")
    pqvspcorrectSummary <- data.frame(pqvspcorrect)
  }

  if(optq=="yes"){
    # output <- list(nunknouwnOut,NROut,psStatsSummary,pqvspcorrectSummary,timeMessage,NL,quantile,pqvspcorrectPlot)
    output <- list(nunknouwnOut,NROut,psStatsSummary,pus_listDF,timeMessage,NL,quantile,pqvspcorrectSummary)
  }
  if(optq=="no"){
    output <- list(nunknouwnOut,NROut,psStatsSummary,pus_listDF,timeMessage,NL,quantile)
  }

  pbcount <- pbcount+1
  if(verbose) setTxtProgressBar(pb,pbcount)

  return(output)

}


