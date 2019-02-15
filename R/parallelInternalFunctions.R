####################################
#### Paralelle inner functions #####
####################################

##### computeChipProfile #####
### Looping over Loci###
.internalChIPLoci <- function(profile,Occup,LocalSet,OccupancyVals,chipMean=chipMean,chipSd=chipSd,
  stepSize=stepSize,norm=norm,chipSmooth=chipSmooth,peakSignificantThreshold=peakSignificantThreshold,ZeroBackground=ZeroBackground,removeBackground=removeBackground,method=method){
  #browser()
    stepIndex<-seq(from=1, to=width(LocalSet), by=stepSize)
    occupancyAbundanceChIPLocal<-rep(ZeroBackground,width(LocalSet))

    occupancyAbundanceChIPLocal[(start(Occup) - start(LocalSet) + 1)] <- OccupancyVals

    occupancyAbundanceChIPLocal[which(occupancyAbundanceChIPLocal <
        removeBackground)] <- 0

    if(any(method=="exact")){
        occupancyAbundanceChIPLocal <- .generateChIPProfile(
        occupancyAbundanceChIPLocal,chipMean,chipSd,
        chipSmooth,norm=norm, quick=FALSE,
        peakSignificantThreshold=peakSignificantThreshold)
    }else if(any(method=="truncated_kernel")){
        occupancyAbundanceChIPLocal <- .generateChIPProfile(
        occupancyAbundanceChIPLocal,chipMean,chipSd,
        chipSmooth,norm=norm, quick=TRUE,
        peakSignificantThreshold=peakSignificantThreshold)
    } else if(any(method=="moving_kernel")){
        occupancyAbundanceChIPLocal <- .generateChIPProfileRcpp(
        occupancyAbundanceChIPLocal, chipMean, chipSd, chipSmooth,
        norm = norm,
        peakSignificantThreshold=peakSignificantThreshold)
    }
    profile$ChIP<-occupancyAbundanceChIPLocal[stepIndex]
    return(profile)
}
##### computeChipProfile #####
### Looping over Parameters ###
.internalChIPParam <- function(SplitGRList,Occup,LocalSet,OccupancyVals,chipMean=chipMean,chipSd=chipSd,
  stepSize=stepSize,norm=norm,chipSmooth=chipSmooth,peakSignificantThreshold=peakSignificantThreshold,ZeroBackground=ZeroBackground,removeBackground=removeBackground,method=method){

  for(j in seq_along(Occup)){
      stepIndex<-seq(from=1, to=width(LocalSet[[j]]), by=stepSize)
      occupancyAbundanceChIPLocal<-rep(
          ZeroBackground,width(LocalSet[[j]]))
      occupancyAbundanceChIPLocal[(start(Occup[[j]]) -
          start(LocalSet[[j]]) + 1)] <- OccupancyVals[[j]]
      occupancyAbundanceChIPLocal[which(occupancyAbundanceChIPLocal <
          removeBackground)] <- 0
      if(any(method=="exact")){
      occupancyAbundanceChIPLocal <- .generateChIPProfile(
          occupancyAbundanceChIPLocal,chipMean,chipSd,
          chipSmooth,norm=norm, quick=FALSE,
          peakSignificantThreshold=peakSignificantThreshold)
      }else if(any(method=="truncated_kernel")){
      occupancyAbundanceChIPLocal <- .generateChIPProfile(
          occupancyAbundanceChIPLocal,chipMean,chipSd,
          chipSmooth,norm=norm, quick=TRUE,
          peakSignificantThreshold=peakSignificantThreshold)
      } else if(any(method=="moving_kernel")){
      occupancyAbundanceChIPLocal <- .generateChIPProfileRcpp(
          occupancyAbundanceChIPLocal, chipMean, chipSd, chipSmooth,
          norm = norm,
          peakSignificantThreshold=peakSignificantThreshold)
      }
      SplitGRList[[j]]$ChIP<-occupancyAbundanceChIPLocal[stepIndex]

  }

  return(SplitGRList)
}

##### computeChipProfile #####
### Split setSequence for parallel processing in computeChipProfile
.internalChIPLociSplit<- function(LocalSet,stepSize){

  stepIndex <- seq(from=1, to=width(LocalSet), by=stepSize)
  SplitSeq <- GRanges(seqnames=S4Vectors::Rle(as.character(
      seqnames(LocalSet)),
      length(stepIndex)),
      ranges = IRanges::IRanges(start = start(LocalSet)+
          stepIndex-1,end =start(LocalSet)+stepIndex-1+(stepSize)),
      strand = "*", ChIP=rep(0,length(stepIndex)))

  names(SplitSeq) <- rep(names(LocalSet), length(SplitSeq))
  return(SplitSeq)
}

##### computeChipProfile #####
### Extracting Occupancy in parallel computeChipProfile if using Parameters
.internalChIPOccupValsParam <- function(Occup){
    OccupVals <- lapply(Occup, function(x){
        x<-x$Occupancy})
    return(OccupVals)
}

##### computeChipProfile #####
### Extracting Occupancy in parallel computeChipProfile if using Loci
.internalChIPOccupValsLoci <- function(Occup){
    OccupVals <- Occup$Occupancy
    return(OccupVals)
}

##### computeGenomeWidePWMScore #####
### Splitting DNA Sequence set Custom function
.splitDNARanges <- function(DNASequenceSet,cores){
    SplitSeq <- floor(seq(1,length(DNASequenceSet),length.out=cores+1))
    splitDNASequenceSet<-vector("list",(cores))
    start <- SplitSeq[1:(length(SplitSeq)-1)]
    end <- c(SplitSeq[2:(length(SplitSeq)-1)]-1,SplitSeq[length(SplitSeq)])

    for(i in seq_len(cores)){
        splitDNASequenceSet[[i]]<-DNASequenceSet[start[i]:end[i]]
    }
    return(splitDNASequenceSet)
}


##### computeGenomeWidePWMScore #####





##### computeGenomeWidePWMScore #####
## DNASequenceSet subset from Access
.internalDNASequenceFromAccess<-function(DNA,Access){

    ## Drop unneccesary levels
    ## surprise surprise it was a factor problem
    drop<-seqlevels(Access)[which(!seqlevels(Access) %in% seqnames(Access))]
    Access<-dropSeqlevels(Access,drop)


    ## Extract DNASequenceSet with respect to Accessibility
    sub<-getSeq(DNA,Access)
    names(sub)<-as.character(seqnames(Access))


    return(sub)
}

##### computePWMScore #####
#### Scores Above threshold

.internalPWMScoreExt<-function(setSequence,DNASequenceSet,DNAAccessibility,PWM,PWMThresholdLocal,
    minPWMScore,maxPWMScore,strand,strandRule){

    ## Building storage lists
    BufferSequence <- vector("list", length(setSequence))
    NoAccess <- c("-")
    ## Naming if names = NULL
    if(is.null(names(setSequence))){
        names(setSequence)<-paste0(as.character(seqnames(setSequence)),":",start(setSequence),"..",end(setSequence))
    }







    for(i in seq_along(setSequence)){

        ## ID matching
        localDNASequenceSet <- DNASequenceSet[unique(match(as.character(seqnames(setSequence))[i],
            names(DNASequenceSet)))]

        ## Extracting Sequence
        #BufferSequence[[i]]<-.internalDNASequenceFromAccess(localDNASequenceSet,setSequence[i],PWM)
        BufferSequence[[i]] <- DNAStringSet(localDNASequenceSet[[
            which(names(localDNASequenceSet)==(as.character(seqnames(
                setSequence[i]))))]],start = start(setSequence[i]),end = end(setSequence[i]))

        ## Scoring Sequence with PWM
        BufferSequence[[i]] <- .scoreDNAStringSet(PWM, BufferSequence[[i]],
        strand=strand,strandRule=strandRule)
        ## Extracting Sites Above Threshold
        indexPWMThresholded <- .getIndexOfPWMThresholded(BufferSequence[[i]][[1]],PWMThresholdLocal)[[1]]

        ## Extracting strand Information for sites above threshold
        if(strand == "+-" | strand == "-+"){
            strandLocal <- rep("*", length(BufferSequence[[i]][[1]][[1]]))
            strandLocal[BufferSequence[[i]][[2]][[1]]] <- "+"
            strandLocal[BufferSequence[[i]][[3]][[1]]] <- "-"
        }
        if(strand == "+"){
            strandLocal <- rep("+",length(BufferSequence[[i]][[1]][[1]]))
        }
        if(strand == "-"){
            strandLocal <- rep("-",length(BufferSequence[[i]][[1]][[1]]))
        }

        ## Building GRanges from indexed PWMScore
        #if(length(indexPWMThresholded)>=1){
        GRLocal <- GRanges(seqnames = S4Vectors::Rle(unique(as.character(
            seqnames(setSequence[i]))),
                length(indexPWMThresholded)),
            ranges = IRanges::IRanges(start=(
                start(setSequence[i]) +
                indexPWMThresholded-1),
                end = (start(setSequence[i])+indexPWMThresholded-
                1 + ncol(PWM)-1)),
            strand = strandLocal[indexPWMThresholded],
            PWMScore=
            BufferSequence[[i]][[1]][[1]][indexPWMThresholded])

        names(GRLocal)<-rep(names(setSequence)[i],length(GRLocal))
        BufferSequence[[i]]<-GRLocal
      #  } else {
          #  GRLocal<-GRanges()

          #  BufferSequence[[i]]<-GRLocal
      #  }

        ## Intersection with DNA Accessibility
        if(!is.null(DNAAccessibility)){
            AccessibleScore <- as.data.frame(findOverlaps(BufferSequence[[i]],DNAAccessibility))


            BufferSequence[[i]]<-BufferSequence[[i]][AccessibleScore[,1]]

            BufferSequence[[i]]$DNAaffinity<-DNAAccessibility$DNAaffinity[AccessibleScore[,2]]
            if(length(BufferSequence[[i]])<1){
                NoAccess <-c(NoAccess,names(setSequence)[i])
            }
        } else {
              BufferSequence[[i]]$DNAaffinity<-rep(0,length(BufferSequence[[i]]))
        }

    }
    return(list(BufferSequence,NoAccess))
}



### Data processing extraction

.internalExtraction<-function(setSequence,profile,maxSignal,backgroundMethod){
    sub<-.extractOccupancyDataAtLoci(profile, setSequence, maxSignal=maxSignal,backgroundMethod=backgroundMethod)
    return(sub)
}
