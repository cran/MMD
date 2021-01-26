#--- Documentation  (roxygen2 needed to compile)
#--
#' Loci redundancy in sequences
#'
#' @param datafile character; Name of the file *.csv (with full path in the file system) containing the genotypes (features) of individuals.
#' @param popfile character; Name of the file *.pop (with full path in the file system) containing the genotypes (features) of individuals.
#' @param NL integer; number of loci. If larger than the number of available loci in the data set, NL is reduced to the maximum available number of loci.
#' @param sourcenames a character vector listing the names of the sources.
#' @param verbose boolean (TRUE/FALSE) for the display of a progress bar  (Default FALSE)
#'
#' @return
#NL0,Nsamples,cbind(Sources,qs),nallele,alleles,hist_s_nl,plot_Rn_original,cbind(c(2:NL),Rn_original),plot_Rn_sorted,cbind(c(2:(NL-1)),Rn_sorted[c(2:(NL-1))]),newfiles,timeMessage
#' A list with
#' \enumerate{
#' \item Number of loci.
#' \item Number of individuals in sources.
#' \item Data frame with proportional weight of each population, qs.
#' \item Number of alleles in the dataset.
#' \item Value of the alleles in the dataset.
#' \item Dataframe with two columns: (1) Index of locus. (2) Rn for loci in the original dataset.
#' \item numerical; index of loci with increasing Rn.
#' \item Dataframe with two columns: (1) Index of loci sorted by increasing Rn. (2) Value of Rn in increasing order.
#' \item Runtime.
#' }
#'
#' @examples
#' ## This example uses a small dataset stored in the MMD package
#' datafile <- system.file("extdata", "Campylobacter_10SNP_HlW.csv", package = "MMD")
#' popfile <- system.file("extdata", "Campylobacter_10SNP_HlW.pop", package = "MMD")
#'
#' NL <- 100
#' sourcenames <- c("Cattle","Chicken","Pig","Sheep","WB")
#'
#' RedundancyLoci <- MMD_Rn(datafile,popfile,NL,sourcenames)
#'
#' ## See more detailed examples in the vignette.
#'
#' @author Francisco J. Perez-Reche (Univeristy of Aberdeen)

MMD_Rn <- function(datafile,popfile,NL,sourcenames,verbose=FALSE){

  start_timeMSD <- Sys.time()


  # - Progress bar
  if(verbose) pb <- txtProgressBar(min = 0, max = (length(sourcenames)+5), initial = 0,style=3)

  rawdata <- read.big.matrix(datafile,sep=",",type="integer",header = F)
  popdata <- as.matrix(read.table(popfile))
  NL <- min(NL,ncol(rawdata))

  if(verbose) setTxtProgressBar(pb,1)
  pbcount <- 1

  # --- Select only genotypes from sources
  indexsources <- which(popdata %in% sourcenames)
  rawdata <- rawdata[indexsources,c(1:NL)]

  NumberLoci <- ncol(rawdata)
  Nsamples <- nrow(rawdata)
  ListNames <- popdata[indexsources]
  reservoirs <- as.matrix(unique(ListNames))
  Sources <- intersect(reservoirs,sourcenames)
  NR <- length(Sources)

  aux1 <-ListNames
  NumberIr<-vector("integer",NR)
  for (k in 1:NR){
    NumberIr[k]<-sum(aux1==Sources[k])
  }

  # Proportional weight of each source
  qs=NumberIr/sum(NumberIr)

  #aux2<-data.matrix(subset(rawdata,rawdata[,1]==Sources[k])[,c(2:(NumberLoci+1))])
  #aux2<-data.matrix(rawdata[,c(2:(NumberLoci+1))])

  aux2<-rawdata

  alleles<-unique(subset(as.vector(aux2),aux2>0))  # -- Vector of alleles, neglecting negative values

  nallele<-length(alleles)


  #--------------------------------
  #-- Mutual information         --
  #--------------------------------
  Ill <- array(0,c(NL,NL))
  for(l1 in 1:NL){
    for(l2 in l1:NL){
      factorIll <- array(0,c(max(alleles),max(alleles)))
      a1<-cbind(rawdata[,l1],rawdata[,l2])
      aux1 <- count(subset(a1,a1[,1]>0 & a1[,2]>0))
      #pi_alal <- data.matrix(aux1[,3]/sum(aux1[,3]))
      aux1[,3] <- aux1[,3]/sum(aux1[,3])
      pi_alal <- aux1[,3]
      for(ta in 1:nrow(aux1)){
        #for(ta2 in 1:length(nrow(aux1))){
        pi_al1 <- sum(subset(aux1[,3],aux1[,1]==aux1[,1][ta]))
        pi_al2 <- sum(subset(aux1[,3],aux1[,2]==aux1[,2][ta]))
        factorIll[aux1[,1][ta],aux1[,2][ta]]=pi_alal[ta]*log2(pi_alal[ta]/(pi_al1*pi_al2))
        #}
      }

      Ill[l1,l2] <- sum(as.vector(factorIll))
      Ill[l2,l1] <- sum(as.vector(factorIll))
    }
    pbcount <- pbcount+1
    if(verbose) setTxtProgressBar(pb,pbcount) #pb-b
  }

  ns <- which(diag(Ill)>0)
  maxns <- max(ns)

  NL0 <- NL
  NL <- length(ns)
  ### Define s_{n|l}
  s_nl <- array(-1,c(NL,NL))
  for(l1 in 1:NL){
    for(l2 in l1:NL){
      s_nl[l1,l2] <- Ill[ns[l1],ns[l2]]/Ill[ns[l1],ns[l1]]
      s_nl[l2,l1] <- Ill[ns[l2],ns[l1]]/Ill[ns[l2],ns[l2]]
    }
    pbcount <- pbcount+1
    if(verbose) setTxtProgressBar(pb,pbcount)  #pb-c
  }

  #---------------------------------------------------
  # Redundancy of ordering in the original file, Rn --
  #---------------------------------------------------
  Rn_original <- array(1,NL-1)
  for(i in 1:(NL-1)){
    Rn_original[i]=max(s_nl[i+1,c(1:i)])
  }


  #------------------------------------------------
  #-- Greedy pair redundancy sorting             --
  #-- Based on R_n = max{s_{s|l}, l=1,2,...,n-1} --
  #------------------------------------------------
  all_l <- c(1:NL)
  sortedl <-array(0,NL)
  l2 <- 1
  Rn_sorted <- array(0,NL)
  sortedl[1]=1

  for(l in 1:(NL-1)){
    remain_l<-as.vector(setdiff(all_l,sortedl))
    a <- s_nl[sortedl[c(1:l)],remain_l]
    if(l>1 && l<NL-1) Rn_remain <- apply(a,2,max)
    if(l==1) Rn_remain <- a
    if(l==(NL-1)){
      # Rn_remain <- a
      # print(a)
      # print(remain_l)
      Rn_sorted[l+1]=max(a)
      sortedl[l+1]=remain_l
      }
    if(l<(NL-1)){
      Rn_sorted[l+1]=min(Rn_remain)
    #print(cbind(l,min(Rn_remain)))
    argmin_Rn_remain <- which(Rn_remain==min(Rn_remain))
    l2 <- remain_l[argmin_Rn_remain[1]]
    sortedl[l+1] <- l2
    }
    pbcount <- pbcount+1
    if(verbose) setTxtProgressBar(pb,pbcount)  #pb-d
    #print(l2)
  }

    end_timeMSD <- Sys.time()
    timeMessage <- end_timeMSD-start_timeMSD

  qsDF <- as.data.frame(cbind(Sources,qs))
  colnames(qsDF) <- c("Sources","qs")
  Rn_originalDF <- as.data.frame(cbind(c(2:NL),Rn_original))
  colnames(Rn_originalDF) <- c("locus","Rn_original")
  Rn_sortedDF <- as.data.frame(cbind(c(2:NL),Rn_sorted[c(2:NL)]))
  colnames(Rn_sortedDF) <- c("locus","Rn_sorted")

  output <- list(NL0,Nsamples,qsDF,nallele,alleles,Rn_originalDF,sortedl,Rn_sortedDF,timeMessage)

  pbcount <- pbcount+1
  if(verbose) setTxtProgressBar(pb,pbcount)  #pb-e

  return(output)
  }
