#--- Documentation  (roxygen2 needed to compile)
#--
#' Loci entropies to measure allele diversity
#'
#' @param datafile character; Name of the file *.csv (with full path in the file system) containing the genotypes (features) of individuals.
#' @param popfile character; Name of the file *.pop (with full path in the file system) containing the genotypes (features) of individuals.
#' @param NL integer; number of loci. If larger than the number of available loci in the data set, NL is reduced to the maximum available number of loci.
#' @param sourcenames a character vector listing the names of the sources.
#' @param verbose boolean (TRUE/FALSE) for the display of a progress bar  (Default FALSE)
#'
#' @return
#' A list with
#' \enumerate{
#' \item Number of loci.
#' \item Number of individuals in sources.
#' \item Table with proportional weight of each population, qs.
#' \item Number of alleles in the dataset.
#' \item Value of the alleles in the dataset.
#' \item Data frame with three columns for entropies: HlT, HlW, HlB
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
#' EntropyLoci <- MMD_Entropy(datafile,popfile,NL,sourcenames)
#'
#' ## See more detailed examples in the vignette.
#'
#' @author Francisco J. Perez-Reche (Univeristy of Aberdeen)


MMD_Entropy <- function(datafile,popfile,NL,sourcenames,verbose=FALSE){

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


  #-----------------------------------------
  # ---- Total Shannon entropy of loci -----
  #-----------------------------------------
  HlT <- vector("numeric",NL)
  for (l in 1:NL){
    factorHlT <- array(0,max(alleles))
    #aux7 <- data.matrix(rawdata[,l])
    aux7 <- rawdata[,l]
    aux8 <- data.matrix(count(subset(aux7,aux7>0)))

    #aux8 <- data.matrix(count(aux7))
    aux9 <- sum(aux8[,2])
    aux10 <- aux8[,2]/aux9

    #pia <- array(0,c(nallele))

    for (ta in 1:length(aux10)){

      factorHlT[aux8[,1][ta]] <- -aux10[ta]*log2(aux10[ta])
      #   pia[aux4[,1][ta]] <- aux10[ta]
    }
    HlT[l]=sum(factorHlT)
    pbcount <- pbcount+1
    if(verbose) setTxtProgressBar(pb,pbcount)
  }

  nosingle<-which(HlT>0) #Selecting loci with positive entropy. Those with zero entropy cannot possibly distinguish different sources.


  #----------------------------------
  # ---- Within-sources entropy -----
  #----------------------------------
  HlW <- vector("numeric",NL)
  icount <- 0
  for (l in 1:NL){
    #pias <- array(0,c(nallele,NR))
    factorHlW <- array(0,c(max(alleles),NR))
    #  pia <- array(0,c(nallele))
    #  factorpia <- array(0,c(nallele,NR))
    icount <- icount+1
    for (k in 1:NR){

      indk <- which(ListNames==Sources[k])
      aux3 <- rawdata[indk,l] #data.matrix(subset(rawdata,popdata==Sources[k])[l])
      aux4 <- data.matrix(count(subset(aux3,aux3>0)))
      #      aux4 <- data.matrix(count(aux3))
      aux5 <- sum(aux4[,2])
      aux6 <- aux4[,2]/aux5
      for (ta in 1:length(aux6)){
        factorHlW[aux4[,1][ta],k] <- -qs[k]*aux6[ta]*log2(aux6[ta])
        #pias[aux4[,1][ta],k] <- aux6[ta]
        #        if(length(aux6)>4) print(length(aux6))
      }
    }
    HlW[l]=sum(as.vector(factorHlW))
    pbcount <- pbcount+1
    if(verbose) setTxtProgressBar(pb,pbcount)
  }

  #-----------------------------------
  # ---- Between-sources entropy -----
  #-----------------------------------
  HlB=HlT-HlW

  end_timeMSD <- Sys.time()
  timeMessage <- end_timeMSD-start_timeMSD

  entropies <- as.data.frame(cbind(HlT,HlW,HlB))
  colnames(entropies) <- c("HlT","HlW","HlB")
  output <- list(NL,Nsamples,cbind(Sources,qs),nallele,alleles,entropies,timeMessage)

  pbcount <- pbcount+1
  if(verbose) setTxtProgressBar(pb,pbcount)

  return(output)
}
