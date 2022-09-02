####################################
### Paralelle internal functions ###
####################################

### General Functions ###

# splitting by number of cores
.splitRanges <- function(ranges,cores){

    if(cores>1 & cores<length(ranges)){
        SplitRanges <- floor(seq(1,length(ranges),length.out=cores+1))
        rangeSet<-vector("list",(cores))
        start <- SplitRanges[1:(length(SplitRanges)-1)]
        end <- c(SplitRanges[2:(length(SplitRanges)-1)]-1,SplitRanges[length(SplitRanges)])

        for(i in seq_len(cores)){
            rangeSet[[i]]<-ranges[start[i]:end[i]]
        }

  } else if(cores>1 & cores==length(ranges)) {
    idx<-seq_along(ranges)
    rangeSet<-vector("list",cores)
    for(i in seq_len(cores)){
        rangeSet[[i]]<-ranges[idx[i]]
    }
  } else if(cores>1 & cores>length(ranges)){
    idx<-seq_along(ranges)
    rangeSet<-vector("list",length(ranges))
    for(i in seq_along(ranges)){
        rangeSet[[i]]<-ranges[idx[i]]
    }
    cores<-length(ranges)


  }  else{
    rangeSet <-list(ranges)
  }
    return(list(rangeSet=rangeSet,cores=cores))
}


### processingChIP function ###

# extracting score from ChIP profiles at loci

.extractionChIP <- function(profile,loci,noiseFilter="sigmoid"){

    # build dem profiles
    # Extract the score
    idx<-seq_along(loci)
    localProfiles <- lapply(idx,.scoreReplace,loci, profile)

    ## Extracting noise filter threshold and assinging weights
    if(grepl("zero",noiseFilter,ignore.case=TRUE)){
        Background<-0
    }

	  if(grepl("mean",noiseFilter,ignore.case=TRUE)){
		    Background<-mean(sapply(localProfiles,mean),na.rm=T)
		}

		if(grepl("median",noiseFilter,ignore.case=TRUE)){
			  Background<-median(sapply(localProfiles,median),na.rm=T)
		}

		if(grepl("sigmoid",noiseFilter,ignore.case=TRUE)){

      Background<-0
			midpoint<-mean(sapply(localProfiles,quantile,0.95),na.rm=T)
			localProfiles<-lapply(localProfiles,.sigmoidWeight,midpoint=midpoint)

		}
    ## Actually applying the noise filter
    localProfiles<-lapply(localProfiles, function(x){x[x < Background]<-0;return(x)})
    if(is.null(names(loci))){
    names(localProfiles)<-paste0(as.character(seqnames(loci)),
                                 ":",start(loci),"..",end(loci))
    } else {
       names(localProfiles)<-names(loci)
    }
    return(localProfiles)
}
##### computePWMScore #####
#### Scores Above threshold

.internalPWMScoreExt<-function(loci,DNASequenceSet,chromatinState,PWM,PWMThresholdLocal,
    minPWMScore,maxPWMScore,strand,strandRule){

    ## Building storage lists
    BufferSequence <- vector("list", length(loci))

    ## Naming if names = NULL
    if(is.null(names(loci))){
        names(loci)<-paste0(as.character(seqnames(loci)),":",start(loci),"..",end(loci))
    }

    for(i in seq_along(loci)){

        ## ID matching
        localDNASequenceSet <- DNASequenceSet[unique(!is.na(match(as.character(seqnames(loci))[i],
            names(DNASequenceSet))))]

        ## Extracting Sequence
        #BufferSequence[[i]]<-.internalDNASequenceFromAccess(localDNASequenceSet,loci[i],PWM)
        BufferSequence[[i]] <- DNAStringSet(localDNASequenceSet[[
            which(names(localDNASequenceSet)==(as.character(seqnames(
                loci[i]))))]],start = start(loci[i]),end = end(loci[i]))

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
            seqnames(loci[i]))),
                length(indexPWMThresholded)),
            ranges = IRanges::IRanges(start=(
                start(loci[i]) +
                indexPWMThresholded-1),
                end = (start(loci[i])+indexPWMThresholded-
                1 + ncol(PWM)-1)),
            strand = strandLocal[indexPWMThresholded],
            PWMScore=
            BufferSequence[[i]][[1]][[1]][indexPWMThresholded])

        names(GRLocal)<-rep(names(loci)[i],length(GRLocal))
        BufferSequence[[i]]<-GRLocal

        ## Intersection with DNA Accessibility
        if(!is.null(chromatinState)){
            AccessibleScore <- as.data.frame(findOverlaps(BufferSequence[[i]],chromatinState))


            BufferSequence[[i]]<-BufferSequence[[i]][AccessibleScore[,1]]

            BufferSequence[[i]]$DNAaffinity<-chromatinState$DNAaffinity[AccessibleScore[,2]]

        } else {
              BufferSequence[[i]]$DNAaffinity<-rep(1,length(BufferSequence[[i]]))
        }

    }
    names(BufferSequence)<-names(loci)
    return(BufferSequence)
}


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

  ## Honestly you could do this as an lapply as well
  ## No excuse PAtrick!
  ## yeah but time is my excuse


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
          stepIndex-1,end =start(LocalSet)+stepIndex-1+(stepSize-1)),
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



#### Accuracy Estimate #####
.intMSE<-function(gp,chip){
     mse <- sum((1/length(gp))*(gp-chip)^2)
     return(mse)
}

.GoFCheck<-function(gp,chip){
    gpsd<-sd(gp,na.rm=TRUE)==0
    chipsd<-sd(chip,na.rm=TRUE)==0
    gp[which(is.na(gp))]<-0
    chip[which(is.na(chip))]<-0
    return(list("zeroCheck"=c(gpsd,chipsd),"predicted"=gp,"validation"=chip))
}

.GoFPred<-function(GoF,method="all",step=10){
    localGoF<-lapply(GoF,function(x,method,step){
                     return(.methodSwitchGoF(x$prediction,x$validation,method,step=step))
    },method=method,step=step)
    return(localGoF)
}

.GoFLoci <-function(GoF,method="all",step=10,cores=1){
  localGoF<-parallel::mclapply(GoF,function(x,method,step){
                   return(.methodSwitchGoF(x$prediction,x$validation,method,step))
  },method=method,step=step,mc.cores=cores)
  return(localGoF)
}

.methodSwitchGoF <- function(gp,chip,method,step){

     if(any(grepl(method, c("recall","precesion","fscore","MCC","Accuracy","AUC"),ignore.case=T))){
        method<-"all"
     }
     if(grepl("pearson", method)){
         check <- .GoFCheck(gp,chip)
         if(length(which(check$zeroCheck==FALSE))>0){
           local<-c(0,.intMSE(gp,chip))
         }else{
           local<-c("pearson"=cor(check$predicted,check$validation,method="pearson"),"MSE"=.intMSE(check$predicted,check$validation))
         }

     } else if(grepl("spearman",method)){

       check <- .GoFCheck(gp,chip)
       if(length(which(check$zeroCheck==FALSE))>0){
         local<-c(0,.intMSE(gp,chip))
       }else{
         local<-c("spearman"=cor(check$predicted,check$validation,method="spearman"),"MSE"=.intMSE(check$predicted,check$validation))
       }
     } else if(grepl("kendall",method)){
       check <- .GoFCheck(gp,chip)
       if(length(which(check$zeroCheck==FALSE))>0){
         local<-c(0,.intMSE(gp,chip))
       }else{
         local<-c("kendall"=cor(check$predicted,check$validation,method="kendall"),"MSE"=.intMSE(check$predicted,check$validation))
       }
     } else if(grepl("ks",method)){
       check <- .GoFCheck(gp,chip)
         local<-c("K-S"=ks.test(check$predicted,check$validation)[[1]],"MSE"=.intMSE(check$predicted,check$validation))
     } else if(grepl("binary",method)){
       check <- .GoFCheck(gp,chip)
         local <-c(.peakExtractionFScore(check$predicted,check$validation),"MSE"=.intMSE(check$predicted,check$validation))
         ## need to check what this goinmg to return
     } else if(grepl("geometric",method)){
         check <- .GoFCheck(gp,chip)
         local<-c("geometric"=.geometricRatio(check$predicted,check$validation,step),"MSE"=.intMSE(check$predicted,check$validation))
     } else if(grepl("all",method)){
          check <- .GoFCheck(gp,chip)

          local <- c(.allMetrics(check$predicted,check$validation,step),"MSE"=.intMSE(check$predicted,check$validation))
     } else if(grepl("MSE",method)){
          check <- .GoFCheck(gp,chip)
          local<-c("MSE"=.intMSE(check$predicted,check$validation))
     }else{
        stop("Make sure that your goodness of fit method is one the available ones \n")
     }
     return(local)
}
